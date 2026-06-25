# LGR (corner-point local grid refinement) — known gaps & test decks

Status of the `opm-gridrefined` LGR rebuild plus its `opm-simulators` /
`opm-common` output integration, as of 2026-06-16. "Works" = exercised and
verified this far; "gap" = throws, is skipped, or produces wrong output.

What already works (for reference):
- Serial LGR: grid build, solve, summary/PRT, EGRID/INIT/UNRST cell output,
  fault crossing an LGR box boundary (`SIMPLE_2PH_W_FAULT_LGR`).
- Parallel LGR (rank-interior boxes): solve (bit-identical global mass
  balance vs serial), summary, **and cell/restart output** (EGRID/INIT/UNRST
  per-grid arrays match serial — added 2026-06-16).
- Wells completed inside an LGR (`COMPDATL`), serial and parallel.

---

## A. Refinement topology (opm-gridrefined `ConformingBlockBuilder`)

These fail at grid build with a clear message. Each already has a *grid-unit*
test in `tests/cpgrid/conforming_builder_test.cpp`; a *flow-level* deck is only
useful as an end-to-end smoke test.

| # | Gap | Where it throws | Deck |
|---|-----|-----------------|------|
| A1 | **Nested LGR** (an LGR refined inside another LGR, `parentGridName != GLOBAL`) — **partially landed 2026-06-16**, see `docs/NESTED_LGR_PLAN.md`. Phase A (simulator parent pass-through + topological ordering) and Phase B (build each nested level grid over its parent LGR) are done and the GLOBAL path is byte-identical; the recursive **leaf** assembly (Phase C) and parallel nested (Phase D) remain, so the builder now throws a precise *"Nested LGR leaf assembly is not implemented"* after building the nested level grids (was an immediate throw at `ConformingBlockBuilder.cpp`). | `ConformingBlockBuilder.cpp` (post-level-grid Phase-C guard) | `opm-tests/lgr/SPE1CASE1_CARFIN1_NESTED.DATA` (CARFIN inside `LGR1`); unit tests `nestedRefinementReachesLeafBoundary`, `nestedChildBeforeParentThrows`. Master supports it. |
| A2 | **Touching boxes, different in-face subdivisions** (different `cellsPerDim` in a *shared/overlap* direction). The builder requires the in-face subdivisions to be **equal**; a different-but-*compatible* (multiple) factor is **also** rejected, even though it is geometrically conformalizable — see the note below. | `ConformingBlockBuilder.cpp` (the per-direction `cellsPerDim` check) | unit test `faceSharingNonMatchingSubdivisionsThrow`; flow decks `opm-tests/lgr/TLGR_VSTACK_HCOMPAT.DATA` (compatible 4-vs-2 → rejected), `TLGR_SIDE_BAD.DATA` (incompatible 3-vs-2). |
| A3 | **Box that cannot be kept on one rank** (spans whole grid, or would empty another rank) — parallel only | `ConformingBlockBuilder.cpp:97` (classifyBox) / `CpGridVanguard::applyLgrPartitionCellGroups_` | `opm-tests/lgr/SPE1CASE1_CARFIN_GR.DATA` (`LGR1` spans the entire grid) — runs serial, throws an actionable error in parallel. Master supports it (distribute-then-refine). |
| A4 | **A refinement box crossing a coarse fault that is *not* on the box boundary in some configs** | — | covered by `SIMPLE_2PH_W_FAULT_LGR` (boundary-crossing fault works); deep fault/pinch interactions inside a box are not separately tested. |

Note: `CARFIN.DATA`/`CARFIN_FLEX.DATA` (diagonal, edge/corner-touching) now build
and are bit-identical to master. `CARFIN_FAULTS.DATA` / `*XYZ-NON.DATA` fail in
opm-common `FAULTS` parsing **on master too** — not a refinement gap.

### A2 detail — compatible touching boxes via a sub-face mosaic (not implemented)

Two refinement boxes that share a 2-D face are conformal only if their
subdivisions match in the **in-face** directions (the two directions *parallel*
to the shared face); the **touch direction** (perpendicular to it) is free. So:

- *Stacked* boxes (touch at a horizontal face): the in-face directions are I,J.
  A different **vertical** (K) factor is fine (e.g. K = 4 vs 2 builds); a
  different **I or J** factor is the constrained case.
- *Side-by-side* boxes (touch at a vertical face): the in-face directions are
  the other two; the along-touch factor is free.

For the constrained (in-face) directions the builder currently demands the two
factors be **equal**. But there are two distinct sub-cases:

- **Compatible** — the finer factor is an integer multiple of the coarser
  (4 vs 2). Geometrically this *is* conformalizable: subdivide the coarser
  side's interface faces into a matching `hi/lo × hi/lo` mosaic so each fine
  cell meets exactly one sub-face. The interface cells on the coarse side become
  hexes with a **split top/bottom (or side) face** — i.e. >6-face cells, exactly
  the multi-face-hex the leaf already uses at a box↔*coarse-neighbour* boundary
  (see [DESIGN-parallel-octree.md](DESIGN-parallel-octree.md) §2/§10). **This is
  a real, supportable case the builder rejects today.**
- **Incompatible** — neither factor divides the other (3 vs 2). The sub-faces
  cannot be aligned by any straight subdivision; the interface is genuinely
  non-conformal and **must** be rejected.

The error message now distinguishes the two (commit on `adaptive-cpgrid`):
*"…different but compatible subdivisions … a sub-face mosaic on the coarser side
… not implemented yet"* vs *"…incompatible subdivisions … neither is a multiple
of the other …"*.

**Why it's not done yet:** the existing mosaic machinery in the leaf assembler
splits a refined box's boundary face against the **parent-level (unrefined)**
neighbour — a single one-level jump. A box↔box interface with a *compatible*
factor needs a mosaic between **two refined levels** (match the finer to the
coarser), which the assembler doesn't construct.

**Implementation effort — medium (a few days), self-contained.** The leaf
representation already permits >6-face cells (no container change). Concretely:
(1) relax the A2 guard from "equal" to "equal **or** the in-face factor of one is
a multiple of the other" in the shared directions; (2) at the shared face,
generate the sub-face mosaic on the coarser box's interface cells by reusing the
parameter-space pillar clip used for the box↔coarse boundary, now matching the
finer box's sub-pillars instead of the parent pillars; (3) wire the resulting
sub-faces into `face_to_cell_`/`cell_to_face_` (variable-length rows already
support it) and dedup the shared corners via the existing exact-coordinate pool.
Risks: only the **compatible** ratios; degenerate/fault cells at the interface
are out of scope (Restriction A); parallel (rank-interior boxes never touch
across ranks) is unaffected. A good first cut is **2-D-in-face mosaics for a
single compatible ratio** (e.g. 2×), validated that a 4-vs-2 stacked pair builds
a conformal leaf and matches a uniformly-4-refined reference where they overlap.

---

## B. Output (opm-simulators / opm-common) — has flow-level test decks

New decks live in `opm-tests/lgr/`, all derived from `SPE1CASE1_CARFIN1.DATA`.

| # | Gap | Where | Test deck | Observed |
|---|-----|-------|-----------|----------|
| B1 | **MINPV with parallel LGR cell output.** The I/O-rank output grid is built with `processEclipseFormat(input_grid, nullptr)` (no `EclipseState`), so the post-MINPV pruning is not applied; for a deck that prunes cells the output grid (full) mismatches the simulation grid (pruned). | `opm-simulators/.../GenericCpGridVanguard.cpp` (outputGrid_ construction) | **`SPE1CASE1_CARFIN1_MINPV.DATA`** (one coarse cell PORO=0 + `MINPV 1.0`, prunes 300→299) | serial fine; **parallel restart PRESSURE differs from serial** (compareECL: 280 errors). Fix: build the output grid through the post-MINPV path. |
| B2 | **RFT output for LGR is skipped.** | `opm-common/.../EclipseIO.cpp:1223` ("RFT file is currently skipped for LGR grids") | **`SPE1CASE1_CARFIN1_RFT.DATA`** (adds `WRFTPLT`) | run completes, **no `.RFT` file written** for the refined grid. |
| B3 | **Restarting (reading) a refined-grid restart is rejected.** | `opm-simulators/.../FlowProblemBlackoil.hpp:1295` (`readEclRestartSolution_`) | **`SPE1CASE1_CARFIN1_RESTART.DATA`** (run `SPE1CASE1_CARFIN1` first, then this — needs `UNIFIN`) | throws *"Refined grids are not yet supported for restart"*. |
| ~~B0~~ | **MPI_RANK in parallel LGR INIT — FIXED (2026-06-16).** Two parts: (a) *sizing* — `MPI_RANK` was written at the full refined-leaf size in the *main-grid* INIT slot (e.g. 924 where the main grid has 300); fixed in `EclGenericWriter_impl.hpp` writeInit (commit a2176e09d) by reducing the leaf `globalRanks_` to level-0 via `getOrigin()`. (b) *per-grid* — the simulator integer maps (incl. `MPI_RANK`) were written only for the main grid via `writeIntegerMaps()`, so the array was absent on refined LGR cells (ResInsight showed it only on the coarse region). Fixed in opm-common `WriteInit.cpp` (commit 2c2e859a5): `writeLGRLocalProperties` now mirrors each integer map onto every LGR via `filterArray(value, global_fathers)`, so each refined cell inherits its father's rank (= the box's owning rank). Verified CARFIN1 np=2: INIT `MPI_RANK` present in all 3 grid sections (main 300 with ranks 0/1, each LGR 324 uniformly its box rank); serial unchanged; parallel UNRST matches serial. |
| B4 | **Parallel LGR trans/NNC in INIT.** The 2026-06-16 fix routes only the *cell-data gather* through the refined output grid; transmissibility/NNC still use the coarse `equilGrid_`. Inter-level `TRANNNC` in the parallel INIT is therefore not yet verified to match serial. | `EclGenericWriter` (equilGrid_ used by `computeTrans_`/`exportNncStructure_`) | reuse any LGR deck; **compare INIT `TRAN*`/`TRANNNC` serial vs np=2** | not yet checked — likely gap. |
| B5 | **Block summary vectors at refined cells** (e.g. `BPR` inside an LGR, ECLIPSE `LGR`-qualified block syntax) | summary config | none yet | unverified. |

---

## C. Parallel correctness / infrastructure

| # | Gap | Notes |
|---|-----|-------|
| C1 | **Distributed wells across an LGR boundary** | rank-interior model keeps each box on one rank; a well perforating cells on both sides of a rank cut, or `--enable-distributed-wells` with LGR, is untested. |
| C1b | **Well completed inside an LGR (`COMPDATL`) at high rank counts** | works at np≤4, **hangs at np≥6** (`SIMPLE_2PH_W_FAULT_LGR`, injector I1 in the WELLI1 box). The well's assigned rank and the box's refining rank diverge, so `compressedIndexForInteriorLGR` finds the connection cells on *no* rank → `ParallelWellInfo.cpp:812` "cells not found" → `checkAllConnectionsFound` throws asymmetrically and deadlocks. Coarse-well LGR decks (e.g. `CARFIN1`) scale to np=8 fine. Partial fix landed (compressedIndexForInteriorLGR returns -1 instead of throwing out_of_range); the remaining fix is to anchor the LGR well to its box's parent coarse cell so its rank matches the box's, and/or make the connection-check collective-safe. **Reproducer:** `mpirun -np 8 flow SIMPLE_2PH_W_FAULT_LGR.DATA ...`. A run should finish in <10 s; a longer one is this hang. |
| C2 | **Timestep-path stability serial vs parallel** | small numbering/geometry differences can move the adaptive controller to a different ministep count → different-length SMSPEC (values agree). Stabilising the controller is a separate task; the regression harness reports these as `Lnz` not FAIL. |
| C3 | **Acceptance harness** | port upstream `tests/cpgrid/lgr/LgrChecks.hpp` invariants (equal cell/intersection geometry, father/siblings, id-consistency) and revive `global_refine` / `lgr_cartesian_idx` / `lgrIJK` tests once nested/parallel features catch up. |
| C4 | **Refine-before-distribute (refine then load-balance) — kept as opt-in `--refine-before-redistribute`, default off (2026-06-16).** The fork uses the *rank-interior* model: load-balance the coarse grid, then refine each box on its owning rank. The alternative (refine globally on rank 0, then distribute the already-refined grid — what upstream master does) is wired in as a selectable option: `Parameters::RefineBeforeRedistribute` (commit `e988172dc`); when on, `CpGridVanguard::loadBalance()` refines `grid_` before `doLoadBalance_` and skips the cell-group partitioning, and `addLgrs()` no-ops if already refined. The corner-point broadcast guard in `CpGrid.cpp` (commit `7c503cf2`) now enters on `comm size > 1` unconditionally (equivalent for the default path, fixes the refine-before asymmetry). **Default stays false / not viable in parallel:** `CpGrid`'s scatter only distributes **level 0** and discards the refinement (the distributed-refinement / refined-grid-scatter machinery was *stripped* from this fork on `strip-lgr`); with the option on in parallel, `addLgrs()` re-refines after scatter and the run hangs in a second `assembleLeafGrid`. Serial / single-rank is correct (verified identical to default). The option is preserved so it becomes the on-switch once refined-grid distribution is reinstated — a substantial feature, separate track. |

---

## How to run the new decks

```sh
FLOW=builds/refined/opm-simulators/bin/flow
CMP=builds/refined/opm-common/bin/compareECL
cd opm-tests/lgr

# B1 MINPV: serial ok, parallel cell output wrong (gap)
$FLOW SPE1CASE1_CARFIN1_MINPV.DATA --parsing-strictness=low --output-dir=/tmp/mp_s
mpirun -np 2 $FLOW SPE1CASE1_CARFIN1_MINPV.DATA --parsing-strictness=low --output-dir=/tmp/mp_p
$CMP -t UNRST -k PRESSURE -n /tmp/mp_s/SPE1CASE1_CARFIN1_MINPV /tmp/mp_p/SPE1CASE1_CARFIN1_MINPV 0.01 1e-3   # expect errors -> gap

# B2 RFT: completes, no .RFT produced
$FLOW SPE1CASE1_CARFIN1_RFT.DATA --parsing-strictness=low --output-dir=/tmp/rft
ls /tmp/rft/*.RFT 2>/dev/null || echo "no RFT (gap)"

# B3 restart: base run then restart -> expected throw
$FLOW SPE1CASE1_CARFIN1.DATA --parsing-strictness=low --output-dir=/tmp/rb
( cd /tmp/rb && $FLOW <repo>/opm-tests/lgr/SPE1CASE1_CARFIN1_RESTART.DATA --parsing-strictness=low --output-dir=. )
# -> "Refined grids are not yet supported for restart"
```

When a gap is fixed, the corresponding deck becomes a positive regression
(compareECL serial==parallel for B1; `.RFT` present + correct for B2; restart
continues for B3).
