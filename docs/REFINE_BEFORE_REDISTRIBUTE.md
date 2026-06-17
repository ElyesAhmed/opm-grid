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

## Status
- [x] Isolated branch + worktrees + `builds/refined_rbr`; fallback preserved.
- [x] Root cause + design verified (this doc).
- [ ] Step 1 — serial refine via self-comm.
- [ ] Step 2 — leaf scatter on stableCellId.
- [ ] Step 3 — wire-up, guards, well check, verification.
