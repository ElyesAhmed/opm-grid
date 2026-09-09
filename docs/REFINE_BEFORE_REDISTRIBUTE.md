# refine-before-redistribute (serial refine, then distribute the leaf)

Branch `refine-before-redist` (off `two_level`). Goal: make
`--refine-before-redistribute=true` work in parallel — refine the full grid on
one processor (serially), then distribute the refined **leaf**, the way the
normal corner-point grid is `process_grdecl`'d on one proc and then
load-balanced. Isolated build: `builds/refined_rbr` (worktrees in
`/Users/hnil/Documents/OPM/rbr_src/`). The `two_level` branch + its compiled
`builds/refined_two_level` are the frozen fallback — do not touch.

## Why it fails today

Two blockers, in order:

1. **classifyBox rejects the box.** With `=true`, `CpGridVanguard::loadBalance`
   calls `addLgrsUpdateLeafView(grid_)` *before* `doLoadBalance_`. `grid_` is the
   replicated global grid (rank 0 owns all cells, other ranks hold them as
   *overlap*). `ConformingBlockBuilder::build` sees `comm().size() > 1`, runs in
   distributed mode, and `classifyBox` throws "box 'CENTER' is split across MPI
   ranks or touches the overlap" because on non-owner ranks the box cells are
   overlap. The refinement must instead run **serially** (refine every cell).

2. **scatterGrid refuses a refined grid.** Even after a correct serial refine,
   `CpGrid::scatterGrid` throws at `CpGrid.cpp:252` ("Loadbalancing a leaf grid
   view with local refinement is not supported, yet") — the distributed-
   refinement machinery was stripped from this fork.

## Key insight — it is NOT the level hierarchy

The leaf *is* a flat corner-point grid (cells/faces/points/geometry) that could
be scattered like level zero. The only real obstacle is that the leaf's
`global_cell_` is **not unique** for refined cells: `LeafGridAssembler.cpp:762`
gives every refined cell its **parent's** Cartesian index (so all children of a
coarse cell share one id), deliberately, so refined cells inherit the parent's
field properties (poro/perm are keyed by Cartesian index). The scatter graph,
`ParallelIndexSet`, and `cartesian→compressed` map all assume `global_cell_` is
a unique per-cell id, so siblings can't be told apart.

**The unique id already exists:** `CpGridData::stableCellId()`
(`CpGridData.cpp:92`) returns a unique int64 per leaf cell — coarse = plain
Cartesian, refined = tagged `(parentCart << 20) | childIdx`. So distribution
needs a *second* id (stableCellId for identity) while keeping `global_cell_`
(parent Cartesian) for field-prop lookup.

## Plan (3 steps)

### Step 1 — serial refinement of the full grid
In `CpGridVanguard::loadBalance` `refineBeforeRedistribute` branch, build the
full refined grid the way `outputGrid_` already does
(`GenericCpGridVanguard.cpp:579`): a self-communicator
(`Dune::MPIHelper::getLocalCommunicator()`) CpGrid copy of the global grid +
`addLgrsUpdateLeafView`. With a self-comm, `ConformingBlockBuilder` runs serial
(no `classifyBox`, all cells refined) and `assembleLeafGrid` skips its
`cc.size() > 1` parallel block (no collectives) — exactly the full refined leaf
a serial run produces. Make this refined grid `grid_` (replacing the coarse one)
so the subsequent `doLoadBalance_` distributes *it*.

### Step 2 — distribute the leaf with stableCellId
Generalise `CpGrid::scatterGrid` (`CpGrid.cpp:208`) + `distributeGlobalGrid`
(`CpGridData.cpp:1541`) to distribute the **leaf** (`selectedLevel = leaf`)
instead of throwing at `CpGrid.cpp:252`:
- partition the leaf cells (graph keyed on `stableCellId()`, not `global_cell_`);
- build `distributed_data_[0]->cellIndexSet()` with `stableCellId()` as the
  global id (the rank-interior leaf already does the equivalent in
  `assembleLeafGrid` lines ~815-833 — fresh per-rank-disjoint ids for refined
  cells; reuse that scheme);
- keep `global_cell_` = parent Cartesian on the scattered cells so
  `distributeFieldProps_` still inherits the parent's properties.

### Step 3 — wire-up + guards
- Drop / condition the `maxLevel()>0` leaf-distribution guard for this path.
- `addLgrs()` already early-returns when `refineBeforeRedistribute && maxLevel>0`.
- Wells: COMPDATL/coarse connections resolve against the distributed leaf
  (`compressedIndexForInteriorLGR` / `compressedIndexForInterior`) — verify.

## Verification
- Default (`=false`) byte-identical (rank-interior path untouched).
- `=true` np=2/4: completes; solution matches `=false` / serial to 1e-6 on
  `SIMPLE_2PH_W_FAULT_LGR` and `SPE1CASE1_CARFIN1`.
- Output (EGRID/INIT/UNRST per-LGR) is a separate track; the solve is the goal.

## Discovered blocker ladder (the flow has more interlocking parts than 3 steps)

Each refine-before blocker found and cleared, in order:
1. **classifyBox** rejected the box — cleared by Step 1 (serial refine).
2. **rank-0-full / rank-1-empty layout.** Before load balancing the global grid
   is on rank 0 (all cells) and *empty* on the others (NOT replicated). Forcing
   every box Owned made empty ranks throw "no active parent cell." Cleared by
   classifying boxes by cell *presence* (`boxCellPresence`): rank 0 refines,
   empty ranks get placeholder level grids (commit afbc48e3).
3. **EQUIL / output-grid setup (CURRENT blocker).** `GenericCpGridVanguard`
   builds `equilGrid_ = CpGrid(*grid_)` and the I/O output grid *before*
   scatterGrid, and `equilGrid()` asserts `mpiRank==0` / assumes a *coarse*
   `grid_`. The pre-refined `grid_` trips this (SIGABRT, line 681) before the
   scatter checkpoint is reached. Needs: build `equilGrid_` from level zero
   (not the refined leaf), and make the output-grid/eclOutputGrid path tolerate
   the already-refined `grid_`.
4. **The leaf scatter itself** (partition propagate -> stableCellId index set +
   interfaces -> distribute) — still ahead, behind (3).

## Status
- [x] Isolated branch + worktrees + `builds/refined_rbr`; fallback preserved.
- [x] Root cause + design verified (this doc).
- [x] **Step 1 — serial refine via self-comm (commit f438eb93).** Added
      `CpGrid::isDistributed()` (`!distributed_data_.empty()`) and gated the
      builder's `distributed`/comm on it. `=true` np=2 now refines the full grid
      (past classifyBox) and reaches the leaf-distribution step; default
      (`=false`) unchanged.
- [ ] **Step 2 — leaf scatter (the deep core; next focused chunk).** Refined,
      lower-risk plan that avoids rewriting the partitioners:
      1. partition **level 0** with the existing machinery (untouched), then
         *propagate* to leaf cells via the parent:
         `leafPart[c] = level0Part[ compressed0(global_cell_[c]) ]` (a leaf
         cell's `global_cell_` is its parent Cartesian). The LGR cell groups
         keep each box's coarse cells (hence its refined leaf cells) on one
         rank, so the propagated partition is rank-interior by construction.
      2. take the distribution communicator from level zero
         (`data_[0]->ccobj_` = world), NOT `data_[leaf]->ccobj_` (which is the
         Stage-1 self-comm).
      3. `distributeGlobalGrid(leaf, leafPart)` scatters the leaf and builds
         `cellIndexSet()` with **`stableCellId()`** as the global id (unique;
         siblings no longer collide), keeping `global_cell_` = parent Cartesian
         for `distributeFieldProps_`.
- [ ] Step 3 — wire-up, guards, well check, verification.

## Stage 2 — concrete code-level plan (from reading scatterGrid + distributeGlobalGrid)

`scatterGrid` (`CpGrid.cpp` ~270-566) and `distributeGlobalGrid`
(`CpGridData.cpp:1541`), for the refine-before leaf, need:

1. **Comm:** `cc = data_[selectedLevel]->ccobj_` (line 270) is the Stage-1 leaf
   *self*-comm (size 1) → the whole `if (cc.size()>1)` block is skipped. Use
   `data_[0]->ccobj_` (world) as the distribution comm for the leaf path. The
   leaf `distributed_data_[0]` is created `make_shared<CpGridData>(worldCc, ...)`.
2. **Partition (propagate, don't re-partition the leaf):** run the existing
   level-0 partitioner to get `level0Part`, then
   `leafPart[c] = level0Part[ compressed0( leaf.global_cell_[c] ) ]`
   (leaf `global_cell_` = parent Cartesian; `compressed0` inverts level-0
   `globalCell()`).
3. **Lists with unique ids:** the index set is built from `importList`
   (`std::get<0>` = global id) at lines 547-553. Build `exportList`/`importList`
   for the leaf **manually** from `leafPart` using `leaf.stableCellId()` as the
   global id (NOT `createListsFromParts`, which keys on the colliding
   `global_cell_`). Owner = `leafPart[c]==rank`; overlap = leaf cells face-
   adjacent (`cell_to_face_`/`face_to_cell_`) to an owned cell but owned by
   another rank (attribute `copy`/`overlap`). This mirrors the rank-interior
   index-set build in `assembleLeafGrid` lines ~795-847, but partition-driven.
4. **distributeGlobalGrid** then scatters cells/faces/points/geometry as-is
   (it reads `view_data.cell_to_point_/cell_to_face_/face_to_point_/geometry_/
   global_cell_` and the index set); `global_cell_` stays parent-Cartesian so
   `distributeFieldProps_` inherits parent props. Verify `global_id_set_` for
   the leaf returns `stableCellId` (else swap it in via the `map2GlobalCellId`
   path at line 1566-1579).
5. **Guards:** replace the `CpGrid.cpp:252` throw with the leaf path for the
   refine-before case; keep it for unsupported combinations.

Test ladder: (a) np=2 distributes without throw; (b) cell counts per rank sane;
(c) solve completes; (d) solution == `=false`/serial to 1e-6 on
`SIMPLE_2PH_W_FAULT_LGR`. Commit each rung.
