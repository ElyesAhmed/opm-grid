# What is needed to redistribute a CpGrid (move the grid from one MPI
# partition to another)

Companion to [REDISTRIBUTION-status.md](REDISTRIBUTION-status.md) (which says
*where we are*); this says *what would have to be built*. Grounded in the
`strip-lgr` source as of 2026-06-23.

"Redistribution" here = **rebalancing an already-distributed grid**: take a grid
that is currently split across ranks under partition P0 and re-split it under a
new partition P1, migrating cells/faces/points (and, with LGR, the refinement
hierarchy) rank→rank. This is distinct from the one-time **distribution**
(serial grid on root → parallel) that `loadBalance`/`scatterGrid` already do.

## TL;DR

- Redistribution does **not exist** today — not with LGR, not even without it.
  It is a **CpGrid-level gap**, independent of refinement.
- The blocker is structural, not a small patch: every piece of the current
  parallel path assumes the source is the **serial grid on rank 0**
  (gather-to-root partitioning, root-authored export/import lists, a
  `distributeGlobalGrid` that reads the global serial view). A
  distributed→distributed migration is a different communication pattern.
- With LGR there is a second, larger layer: the level hierarchy and refined ids
  must also migrate, and the refined-grid scatter machinery was **stripped**
  from this fork (`scatterGrid` throws on a refined leaf).
- Strategic call (already argued in REDISTRIBUTION-status.md): real
  redistribution belongs to the **octree/dynamic-AMR** class via root-tree
  migration, not to retrofitting CpGrid. The breakdown below is the cost of
  doing it *in CpGrid* anyway, should that be wanted as an interim step.

## Naive but sufficient: gather to root, then re-scatter

The "big gap" framing above is about doing redistribution *well* (in place,
distributed→distributed, at scale). If performance and memory are set aside, a
**correct** redistribution is much simpler — it is essentially plumbing, not a
research problem.

The leaf is just plain arrays: topology (ints), geometry (doubles), global
ids/tags. `distributeGlobalGrid` (`CpGridData.cpp:1541`) is *already* nothing
more than a sequence of "communicate this array" handle calls. Nothing in the
grid resists serialization.

So the worst-case recipe is:

> **gather the full distributed grid back to rank 0 → rebuild one serial
> `CpGridData` → call the existing `scatterGrid` with the new partition.**

The only *new* code is the gather (the symmetric inverse of the scatter, which
today nobody writes because the partitioners only gather the lightweight graph
and output is rebuilt from the deck). The scatter itself is reused wholesale.

Two things you cannot skip even in this naive version:

1. **Stable global ids.** The point of rebalancing is to carry the solution
   (and well state) across the move, so every cell must keep its identity
   through gather→rescatter. CpGrid's `global_id_set_` already does this
   bookkeeping — just don't lose the ids (order the gathered serial grid
   deterministically, e.g. by global id / Cartesian).
2. **The LGR hierarchy has to ride along.** Either flatten it to per-leaf-cell
   integer tags `(level, father Cartesian, LGR-local index)` and scatter them
   next to `global_cell_` — the `leafHasParentCellIndices()` "flat-leaf" route,
   which also sidesteps the refined-leaf throw at `CpGrid.cpp:243-268` — or
   gather/rescatter each level grid as its own polyhedral grid plus the maps.

Cost honesty: O(global-grid) memory on rank 0 plus a full gather+scatter per
rebalance. Fine for correctness and small/medium cases; unacceptable at scale.
The octree's root-tree migration (see
[DESIGN-parallel-octree.md](DESIGN-parallel-octree.md) §5, and §10's compact
Layer-A representation — migrating tree *bytes*, not leaf vectors) is the
scalable version of this same idea. It is the *performant* path, **not** a
prerequisite for a correct one.

One-line takeaway: redistribution is plumbing (serialize arrays, gather, reuse
the existing scatter); the hard parts in the rest of this doc are about doing it
*without* the gather-to-root bottleneck.

## Where the code blocks redistribution today

All in `opm/grid/cpgrid/CpGrid.cpp` / `CpGridData.cpp`:

1. **Second-call refusal** — `scatterGrid` (`CpGrid.cpp:230`):
   ```cpp
   if (!distributed_data_.empty()) {
       std::cerr << "There is already a distributed version of the grid. ...";
       return {false, {}};
   }
   ```
   A distributed grid cannot be scattered again. There is **no
   `redistribute()`/`repartition()` method** anywhere in the fork.

2. **Partitioners are gather-to-root only** — `zoltan*PartitionGridOnRoot`,
   `metisSerialGraphPartitionGridOnRoot`, `zoltan*PartitioningWithGraphOfGrid`
   (`CpGrid.cpp:351-388`). They build the graph on **rank 0** and emit a single
   global `computedCellPart`. There is no partitioner that operates on the
   *already-distributed* graph (no ParMETIS/parallel-Zoltan path).

3. **Root-authored export/import lists** — `setupSendInterface` /
   `setupRecvInterface` (`CpGrid.cpp:96,115,557`) are fed export/import tuples
   produced by the on-root partitioner. The send/recv interface is a
   serial→parallel scatter, not a rank↔rank migration.

4. **`distributeGlobalGrid` reads the serial source view**
   (`CpGridData.cpp:1541`): its `view_data` argument is the **global** grid; it
   uses `view_data.global_id_set_`, `view_data.cell_to_point_`,
   `view_data.cellIndexSet()` as the authoritative source topology and scatters
   from it. It cannot take a *distributed* source.

5. **LGR guards** — `scatterGrid` (`CpGrid.cpp:243-268`) throws for a refined
   *level* (`level>0`) and for the *leaf of an LGR'd grid*
   (`maxLevel()>0 && level==-1`); only **level 0** of an LGR grid distributes,
   and only via `zoltanGoG`. The distributed-refinement scatter was removed on
   `strip-lgr`.

## Part 1 — generic redistribution (no LGR). Required regardless of refinement.

1. **A `redistribute(const std::vector<int>& newOwner)` entry** on `CpGrid`
   that accepts the *current distributed* grid as the source and produces a new
   distributed layout. Drop the `!distributed_data_.empty()` refusal for this
   path (or build into a fresh `CpGridData` and swap).

2. **A new (or re-gathered) partition.** Minimum viable: gather the distributed
   graph to root, run the existing on-root partitioner, scatter the new part
   vector. Proper version: a **parallel** partitioner (ParMETIS / Zoltan PHG)
   on the distributed graph so rebalancing scales.

3. **Distributed→distributed migration pattern.** Replace the root-authored
   export/import lists with an **all-to-all** plan: each rank computes, for its
   owned cells, the destination rank from `newOwner`, and the symmetric receive
   side. Dune already has the primitive (`RedistributeInterface` /
   `CommunicationGraph`); the work is wiring CpGrid's cell/face/point/geometry
   handles through it instead of `setupSend/RecvInterface`.

4. **Generalize `distributeGlobalGrid` to a distributed source.** The geometry
   and topology gather (`computeCell2Point`, `computeCell2Face`,
   `computeFace2Cell`, `computeGeometry`, the `scatterData` handles for
   `global_cell_`, face tags/normals/bids) must move data *between two
   distributed layouts* keyed by **global id**, not "from the root's view."

5. **Rebuild everything keyed to the old layout** on the new one: cell/point
   `ParallelIndexSet`, `RemoteIndices`, communication interfaces, the **overlap
   layer** (`addOverlapLayer`), partition-type indicators, and the boundary-id
   set.

6. **Migrate simulator state** (the actual reason to rebalance). The
   solution/primary variables, field properties, well connections,
   transmissibilities and NNCs, aquifer cells, the `cartesianIndexMapper`
   compressed map, and the linear solver's `OwnerOverlapCopyCommunication` are
   all built against P0 and must be re-communicated/rebuilt for P1. This is
   simulator-side (`CpGridVanguard` / `FlowProblem`) work layered on the grid
   primitive above.

## Part 2 — redistribution WITH LGR. On top of Part 1.

7. **Migrate the level hierarchy, not just the leaf.** A redistribute must move
   every level grid (`data_[1..maxLevel]`) plus the per-cell hierarchy maps that
   are currently never scattered:
   `child_to_parent_cells_`, `parent_to_children_cells_`,
   `leaf_to_level_cells_`, `level_to_leaf_cells_`, `corner_history_`, and the
   per-level `global_cell_` (`CpGridData.hpp:379-458`).

8. **Keep refinement boxes coherent under P1.** Either preserve the
   **rank-interior** invariant (a box stays whole on one rank — the partitioner
   must respect box→rank grouping, as `applyLgrPartitionCellGroups_` does for
   the initial distribution), or reinstate the **general cross-rank refined
   scatter** that upstream master has and `strip-lgr` removed (refined-id
   prediction, overlap communication of refined cells, winner-selection). The
   latter is the heavy machinery REDISTRIBUTION-status.md and `LGR_GAPS.md` C4
   refer to.

9. **Globally-consistent refined ids across the migration.** Refined cell /
   point / corner global ids must survive the move so `GlobalIdSet`,
   `father()`/`level()` chains, and inter-level NNCs stay valid on P1.

10. **Lift the LGR scatter guards** (`CpGrid.cpp:243-268`) once 7–9 exist.

## Recommended path (consistent with the existing design docs)

Do **not** retrofit generic distributed→distributed migration into CpGrid for
the static leaf first. Per [REDISTRIBUTION-status.md](REDISTRIBUTION-status.md)
§"How this relates to the octree" and
[DESIGN-parallel-octree.md](DESIGN-parallel-octree.md) §5, the dynamic-AMR
class needs a *different* redistribution anyway — **root-tree migration**:
partition level-0 cells by a space-filling curve, rebalance by moving the
compact per-tree refinement structure (bytes) between ranks, then reconstruct
octant geometry locally from deterministic refinement + the global macro input
(no geometry communication). That mechanism *is* the redistribution capability
CpGrid lacks, at a granularity where it is cheap.

So "add redistribution" and "build the octree" are the same project. The static
parallel LGR in CpGrid is correctly scoped as **distribute-once-then-refine,
never rebalance** — adequate for static CARFIN, and the foundation the octree
generalizes.

If an interim CpGrid redistribution is nonetheless required (e.g. to rebalance
a *non-LGR* run), Part 1 alone is a self-contained, useful deliverable and the
smaller half of the work; Part 2 is where the cost concentrates.
