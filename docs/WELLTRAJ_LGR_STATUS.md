# WELTRAJ/COMPTRAJ wells on a static-LGR grid — status

Status of making well-trajectory (WELTRAJ/COMPTRAJ) wells work on a
**statically-refined (CARFIN/LGR) CpGrid**, by re-intersecting the trajectory
against the refined leaf as a **post-process** (not inside the well-model solve).

Last updated: 2026-06-25.

## Problem

WELTRAJ/COMPTRAJ define perforations by a 3-D trajectory. opm-common turns the
trajectory into connections by intersecting it with the **coarse `EclipseGrid`**
(`WellConnections::loadCOMPTRAJ` → `RigEclipseWellLogExtractor` → the cvf AABB
`BoundingBoxTree`). When a CARFIN box refines the cells the trajectory crosses,
those coarse connections do not correspond to leaf cells, so the well perforates
the wrong/absent cells. The goal is for welltraj to work on the static-LGR grid.

## Approach (chosen, distinct from PR #5914)

- **Post-process** the trajectory → connections against the **refined CpGrid
  leaf**, reusing the existing cvf bbox machinery, fed from explicit leaf-cell
  corners instead of the EclipseGrid.
- Run it once after the grid is refined (start-of-schedule assumption), written as
  a **reusable per-well routine** so a later report-step (re)build can call it too.
- **No recompute in the Newton/assembly loop** — that is PR #5914's rejected
  approach.

Pipeline:
1. opm-common `RigEclipseWellLogExtractorGrid` — geometry-fed sibling of
   `RigEclipseWellLogExtractor`; builds the same cvf `BoundingBoxTree` from
   supplied leaf-cell corners (OPM→ResInsight hex permutation `{0,1,3,2,4,5,7,6}`).
   Keeps opm-common free of any opm-grid dependency.
2. opm-common `WellConnections::recomputeTrajectoryConnections(...)` +
   `Schedule::recomputeTrajectoryConnections(...)` — rewrite a well's
   `WellConnections` with leaf-indexed connections, snapshot get/update write-back.
3. opm-simulators `CpGridVanguard::recomputeWellTrajectoriesInLgr_()` — after
   `addLgrs()`: build leaf corners + per-cell LGR-local ijk, intersect, compute a
   perm-based Peaceman CTF on the refined cell, tag single-LGR wells, and encode
   each connection in opm-common LGR indexing (LGR-local ijk + LGR grid number +
   LGR-cell global index, COMPDATL convention).

See `LGR_GAPS.md` for the surrounding static-LGR gap list.

## Current commits to compile

The feature spans three OPM repos plus this one. Build order:
opm-common → opm-gridrefined → opm-simulators.

| Repo | Branch | HEAD commit | Remote (fork) | Upstream |
|---|---|---|---|---|
**Consolidated (2026-06-25):** welltraj now lives on the integrated `new_lgr`
branch alongside the nested / refine-before / field-props / solver fixes (it used
to be on a separate `welltraj-lgr-postproc` line). Build from `new_lgr` in all
repos; see `STATUS.md` for the full picture.

| Repo | Branch | code HEAD | upstream |
|---|---|---|---|
| opm-common      | `new_lgr` | `6291b8165` | `OPM/opm-common` |
| opm-gridrefined | `new_lgr` | `5030fe10`  | `OPM/opm-grid` |
| opm-simulators  | `new_lgr` | `8970d1b92` | `OPM/opm-simulators` |
| opm-tests       | `new_lgr` | `3b86809c`  | `OPM/opm-tests` |

opm-gridrefined is this workspace's fork of **opm-grid** (`OPM/opm-grid`); the
`adaptive-cpgrid` branch carries the static-LGR + nested-LGR + parallel-LGR work
that the welltraj post-process refines on top of. The welltraj feature itself
lives in opm-common + opm-simulators; it only needs an opm-grid that supports
static CARFIN refinement of CpGrid.

### Key commits (opm-common, `welltraj-lgr-postproc`)
- `96f25e253` geometry-fed extractor `RigEclipseWellLogExtractorGrid`
- `ea76bb1ef` `WellConnections::recompute…` against a supplied grid
- `378c974cb` `Schedule::recompute…` post-process (snapshot write-back)
- `cafedd66f` flag trajectory-LGR wells via `flag_lgr_well()`
- `7e9daa1f7` give trajectory-LGR wells an LGR-local head for output
- `316aebb9c` encode trajectory-LGR connections in opm-common LGR indexing
- `bb42965c2` `RegionCache`: skip LGR connections in the global-grid region cache

### Key commits (opm-simulators, `welltraj-lgr-postproc`)
- `adf1581a5` post-process WELTRAJ connections against the refined leaf
- `c29947231` perm-based CTF for welltraj-LGR connections
- `121e6fefb` give welltraj-LGR connections the opm-common LGR index
- `5952fe789` use the COMPDATL LGR grid-number convention (off-by-one fix)

## Alternative codes / branches

- **OPM/opm-simulators PR #5914** — recomputes trajectory connections **inside**
  `BlackoilWellModel` (in the solve). Rejected here in favour of the post-process;
  kept as the reference for the intersection logic.
- **OPM/opm-common PR #4446** — adds `RigEclipseWellLogExtractorGrid` +
  `recomputeConnections`. Right shape but incomplete upstream (stubbed
  perm/poro/depth, `grid.get_cell` unimplemented, LGR index wiring absent). This
  branch ports and completes it (real refined-leaf perm/geometry, LGR encoding).
- **opm-grid `master`** — has static CARFIN LGR but not the nested/parallel/
  conformity work on `adaptive-cpgrid`; the serial single-LGR welltraj path builds
  on plain master-style static LGR.

## Publicly available test cases (opm-tests)

On opm-tests branch `welltraj-lgr-postproc` (fork `hnil/opm-tests`), reproducible
from public grids:

| Deck | Purpose |
|---|---|
| `weltraj/WELTRAJ-01_CARFIN.DATA` | SPE1 welltraj injector + producer, two CARFIN boxes; validated against `weltraj/WELTRAJ-01.DATA` (no LGR) to <1% |
| `model5/MODEL5_WELTRAJ_CARFIN.DATA` | corner-point grid, producer B-1H welltraj + CARFIN box |
| `lgr/TLGR_VSTACK_42.DATA` | touching LGRs, 4/2 vertical refinement (compatible, builds) |
| `lgr/TLGR_VSTACK_HCOMPAT.DATA` | touching LGRs, different-but-compatible lateral refinement |
| `lgr/TLGR_SIDE_BAD.DATA` | touching LGRs, incompatible refinement → expected clear error |

A second, **local-only** inj+prod verification pair lives in this workspace under
`data/model2_lgr/lgr/` (`TEST_LGR_WT_INJ2_PROD1.DATA` welltraj vs
`TEST_LGR_REF_INJ2_PROD1.DATA` COMPDATL); it depends on local model2 grid
includes and so is not in opm-tests. See `RUNNING_WELLTRAJ_LGR_TESTS.md` in the
workspace root for run/compare instructions.

## What works

- Serial, **single-LGR-per-well** WELTRAJ/COMPTRAJ: perforates the refined leaf
  cells, perm-based CTF, full ECL output (EGRID/INIT/UNRST/SMSPEC).
- Verified to reproduce the explicitly-completed equivalent:
  - SPE1 welltraj+CARFIN vs non-LGR welltraj: well cumulatives/rates <1%.
  - model2 welltraj (PROD1, INJ2 in their own LGRs) vs COMPDATL centre-column on
    the same boxes: WBHP/WOPR/WWPR/WGPR within ~0.2%; an untouched COMPDATL
    control well matches to 0.000%.
- Correct ResInsight placement (after fixing the LGR grid-number off-by-one).
- No regression to COMPDATL / static-LGR / non-LGR welltraj paths (the recompute
  only triggers for refined grids with trajectory wells).

## Known limitations / follow-ups (not done)

- **One LGR per well** — a single well's trajectory must stay inside one CARFIN
  box. Multi-LGR-per-well hits an LGR grid-number convention conflict that would
  break COMPDATL with ≥2 LGRs (prototype reverted).
- **Keep the perf inside the box** — COMPTRAJ perf endpoints on exact cell
  boundaries bleed into the coarse neighbour outside the box → mixed LGR+coarse →
  `Cells … not found in grid`. Pull perf endpoints inward of the box K-range.
- **Parallel (np>1) untested** — field-props local/global indexing on the refined
  leaf not yet validated under domain decomposition.
- Recompute runs per report-step snapshot (could be cached).
- PRT production report shows the LGR-local head (cosmetic; restart IWEL is
  father-converted correctly).
- Cosmetic error-message noise from the vendored `RigWellLogExtractor.cpp`.
- Touching-LGR compatible-but-different mosaic (gap A2 in `LGR_GAPS.md`) — error
  messages in place; the sub-face mosaic build is not implemented (effort: medium).
