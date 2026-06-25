# LGR on CpGrid — status

Status of the local-grid-refinement (CARFIN/LGR) rebuild across the OPM modules.
Companion docs: `lgr_review.md` (review), `PLAN.md` (roadmap), `LGR_GAPS.md`
(gap list), `REDISTRIBUTION-status.md`/`-requirements.md`, `WELLTRAJ_LGR_STATUS.md`,
`NESTED_LGR_PLAN.md`/`NESTED_LGR_TESTING.md`, `REFINE_BEFORE_REDISTRIBUTE.md`,
`DESIGN-parallel-octree.md`/`DESIGN-builder.md` (AdaptiveCpGrid design), and
`RUNNING.md` (workspace root, full build/run/test guide).

Last updated: 2026-06-25.

## Approach

Route A: **strip the LGR layer out of CpGrid and rebuild static refinement as a
separate builder** under `opm/grid/cpgrid/refinement/` (`GrdeclRefinement`
resamples COORD/ZCORN onto sub-pillars; `LevelGridAssembler` runs the
preprocessor per block so faults inside a block are matched by `findconnections`;
`LeafGridAssembler` / `assembleNestedLeafGrid` merge into the leaf;
`ConformingBlockBuilder` is the registered backend). opm-simulators gets only
minimal additive hooks. The parallel model is **rank-interior** (distribute level
0 first, then each rank refines the boxes it owns), with
**refine-before-redistribute** as an optional path for better load balance.

### Geometry note (precise)
Refinement is the **corner-point-native** resampling of COORD/ZCORN (refined
corners lie on straight sub-pillars — ECLIPSE LGR semantics), which conserves the
parent cell's volume and volume-weighted centre of mass. The per-cell **trilinear
map** of a hex's 8 corners refines that single cell's interior correctly for any
hexahedron, and across a **matching, non-faulted** shared face the two neighbours'
refined sub-faces are automatically conformal (both restrict to the same bilinear
surface) — i.e. trilinear alone suffices for **hexes not touching a fault**. Where
a cell touches a fault the two sides are *different* surfaces, so the interface is
matched by intersection via the preprocessor (`findconnections`), not assumed from
the trilinear map. See `LESSONS.md` item 12 and `DESIGN-builder.md`.

## ✅ Now consolidated on one branch: `new_lgr`

Previously the LGR work was split across feature branches (welltraj vs
nested/refine-before) that were never merged, so no single branch had everything.
**As of 2026-06-25 all of it is consolidated onto `new_lgr` in every repo** —
`nested-lgr-parallel` (nested + refine-before + field-props/solver fixes) with the
welltraj + conformity + design-doc work cherry-picked on top. Built and verified
end-to-end (see Verification). Develop on `new_lgr` going forward.

### Branches / commits to compile

| Repo | Branch | code HEAD | fork remote | upstream |
|---|---|---|---|---|
| opm-common      | `new_lgr` | `6291b8165` | `git@github.com:hnil/opm-common.git` | `OPM/opm-common` |
| opm-gridrefined | `new_lgr` | `5030fe10`* | `git@github.com:hnil/opm-gridrefined.git` | `OPM/opm-grid` |
| opm-simulators  | `new_lgr` | `8970d1b92` | `git@github.com:hnil/opm-simulators.git` | `OPM/opm-simulators` |
| opm-tests       | `new_lgr` | `3b86809c`  | `git@github.com:hnil/opm-tests.git` | `OPM/opm-tests` |

\* opm-gridrefined `new_lgr` tip advances with this docs commit; the code tip is
`5030fe10` (conformity). Build order: opm-common → opm-gridrefined(as opm-grid) →
opm-simulators. opm-gridrefined is the fork of `OPM/opm-grid`.

## Capability status (all on `new_lgr`)

| Capability | Status | Notes |
|---|---|---|
| Serial static LGR, bit-identical to master | ✓ | CARFIN1 = 52 Newton |
| Volume/CoM-conserving corner-point refinement | ✓ | more correct than master on skewed/faulted cells |
| Faults inside a CARFIN block | ✓ | preprocessor matches per block |
| Edge/face-sharing (touching) boxes | ✓ | `CARFIN`, `CARFIN_FLEX` |
| Touching-box conformity check | ✓ | equal in-face subdivisions required; compatible→mosaic-not-implemented error; incompatible→clear error |
| Parallel rank-interior single-level LGR | ✓ | CARFIN1 np2; mass balance matches serial |
| LGR ECL output (serial + parallel cell/restart) | ✓ | EGRID/INIT/UNRST consistent, ResInsight-ready |
| **Field-props correct on refined cells** | ✓ | EQLNUM/PVTNUM/SWATINIT (equil), FIPNUM/FPR, datum-region pressure — via LookUpData leaf-mapping |
| THPRES + LGR | ✓ | global restart vectors copied per level |
| Nested LGR (serial, fully contained): build+solve+output | ✓ | `*_NESTED_CONTAINED` 54 Newton |
| Nested LGR (parallel rank-interior) | ✓ | np2 writes EGRID/INIT/UNRST |
| refine-before-redistribute (grid+props+solve+output) | ✓ | `--refine-before-redistribute=true` |
| WELTRAJ/COMPTRAJ on LGR (single-LGR/well, serial) | ✓ | reproduces COMPDATL <0.3%; see `WELLTRAJ_LGR_STATUS.md` |
| Non-black-oil (solvent/polymer/biofilm/MICP) + LGR | ✗→err | clear error; black-oil is the supported scope |
| Whole-grid (box==grid) LGR in parallel | ✗→err | needs distributed refinement; rank-interior refuses clearly |
| Nested boundary-touching child | ✗→err | rejected with containment message |
| Redistribution / rebalancing a distributed grid | ✗ | CpGrid-level gap; belongs to dynamic AMR (`REDISTRIBUTION-status.md`) |
| Dynamic AMR | ✗ | design only (`DESIGN-parallel-octree.md`, AdaptiveCpGrid) |

## New / notable run options

- `--refine-before-redistribute=true` — refine the full grid on rank 0 then
  distribute the leaf (better balance for big LGRs). Default off = rank-interior
  (each rank refines the boxes it owns; a box must fit one rank).
- **Parallel-LGR linear solver** — the rank-interior partition is often very
  unbalanced (whole refined region on one rank), which can empty an AMG coarse
  level. Use one of: `--linear-solver=ilu0` (simplest); a cprw JSON with
  `coarsenTarget` ≥ the first coarse level size (~2000); or a cprw JSON with a
  HYPRE coarse solver (build with `-DUSE_HYPRE=ON`). An empty-partition now throws
  a clear actionable error instead of segfaulting.
- **WELTRAJ/COMPTRAJ** wells are post-processed against the refined leaf so a
  trajectory perforates the correct LGR cells (one LGR per well; keep the
  perforation inside the box). No `--enable-lgr` needed — CARFIN auto-enables.
- **Nested CARFIN** — a box may name a parent LGR; its I/J/K are parent-LGR-local.
- `OPM_LGR_POISON_REFINED=1` — opt-in diagnostic canary (poisons refined cells'
  global Cartesian index to catch any late re-derivation; default off, zero
  impact).
- All CARFIN decks need `--parsing-strictness=low` on any flow build.

## How to run a case

```
# build (dependency order); enable HYPRE if you want the HYPRE coarse solver
ninja -C builds/refined flow_blackoil          # opm-common + opm-grid + flow

FLOW=builds/refined/opm-simulators/bin/flow_blackoil

# serial LGR
$FLOW opm-tests/lgr/SPE1CASE1_CARFIN1.DATA --parsing-strictness=low --output-dir=/tmp/c1

# parallel LGR (rank-interior)
mpirun -np 2 $FLOW opm-tests/lgr/SPE1CASE1_CARFIN1.DATA \
    --parsing-strictness=low --linear-solver=ilu0 --output-dir=/tmp/c1p

# nested LGR
$FLOW opm-tests/lgr/SPE1CASE1_CARFIN1_NESTED_CONTAINED.DATA --parsing-strictness=low --output-dir=/tmp/nest

# refine-before-redistribute
mpirun -np 2 $FLOW opm-tests/lgr/SPE1CASE1_CARFIN1.DATA \
    --parsing-strictness=low --linear-solver=ilu0 --refine-before-redistribute=true --output-dir=/tmp/rbr

# welltraj on LGR
$FLOW opm-tests/weltraj/WELTRAJ-01_CARFIN.DATA --parsing-strictness=low --output-dir=/tmp/wt
```
The workspace `RUNNING.md` has the full build-tree map, the local model2 decks,
the regression harness, and comparison tools.

## Verification (2026-06-25, on `new_lgr` / `builds/refined`)

Grid tests `conforming_builder_test` (16) + `grdecl_refinement_test` (5) green;
serial CARFIN1 = 52 Newton; parallel CARFIN1 np2 = 71 Newton; nested contained =
54 Newton; refine-before np2 = End of simulation; model2 welltraj vs COMPDATL
≤0.3% (with the field-props fix now active); WELTRAJ-01_CARFIN = 88 / model5 = 251
Newton; touching-box decks build (compatible vertical) or reject with the correct
compatible/incompatible messages.

## Remaining gaps (details in `LGR_GAPS.md`)

- **Redistribution** of an already-distributed grid — unsupported (CpGrid-level);
  belongs to the dynamic-AMR class, not a static add-on.
- **Whole-grid LGR in parallel** — needs distributed refinement (box across ranks).
- **refine-before deep reconstruction** for highly irregular leaves — rank-interior
  is the working path for big/faulted cases.
- **Parallel solver empty-partition** — handled by workarounds (ilu0 / coarsenTarget
  / HYPRE); a true min-per-rank AMG redistribution is upstream Dune work.
- **MINPV in the parallel output grid** — not yet handled.
- **Compatible-but-different touching subdivisions (A2)** — error in place; the
  sub-face mosaic build is not implemented (effort medium).
- **Parallel (np>1) welltraj-LGR** — untested.
