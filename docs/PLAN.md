# opm-gridrefined — development plan

Goal: **static LGR for CARFIN decks on general corner-point grids** (faults, pinch-outs), correct and ECL-compatible, eventually upstreamable to OPM/opm-grid as a replacement-with-evidence. Full analysis in [lgr_review.md](lgr_review.md) (Parts I–IV); distilled lessons from the current implementation for the refactoring in [LESSONS.md](LESSONS.md).

## Strategy (Part IV, "Route A")

This repo carries the full opm-grid history (`upstream` = OPM/opm-grid). The plan is *not* a new grid and *not* an in-place refactor of upstream:

- Strip the current LGR/adaptivity layer from CpGrid (keep `mark`/`adapt` as inert stubs; DUNE facade and class name unchanged).
- Rebuild static refinement as a separate **builder component** behind a composition-style interface (§6): LevelHierarchy state + RefinementBuilder.
- Preferred builder: **preprocessor-level (§8.2-S2)** — refinement expressed at the `processed_grid` level, reusing `findconnections()` so faults and pinch-outs are handled by the same code that handles them on level 0; intersections by 2D clipping in pillar parameter space (§8.1, no CGAL).
- `CpGridData` already represents what is needed (variable-length `face_to_point_`, collapsed hexes, >6-face leaf cells) — the work is the build algorithm, not the container.

## Repository structure

1. **opm-gridrefined** (this repo, module name `opm-grid`): grid work, grid-level tests/tools. No opm-simulators dependency.
2. **opm-simulators: untouched** — rebuilt against this module; plain `flow` compiles and runs (without LGR while the builder matures).
3. **Thin app repo later (e.g. opm-flowrefined)**: experimental flow executable (`main()` + TypeTag + custom Vanguard wiring the new builder API) + deck-matrix integration tests. opm-geomech pattern. Deferred until the builder API diverges from `addLgrsUpdateLeafView`.

Fork discipline: rebase regularly on `upstream/master`; the diff stays mostly deletions + the new layer. This is a development vehicle targeting replacement, not a permanent fork. Final upstream PRs staged through hnil/opm-grid.

## Track 0 — save the current implementation (short term, Part II)

Evaluate applying the Part II short-term items to the *existing* LGR code so it stays useful for real decks while the rebuild matures. Per item, decide: small targeted PR to upstream (mergeable, unlike refactors) or carry in this fork.

- [ ] Relax global MULTZ / TRANX-Y-Z checks to "only when intersecting an LGR" (currently any MULTZ + any LGR throws).
- [ ] Pinched/degenerate parents: auto-deactivate inside patches (warn; optionally compensate pore volume).
- [ ] Fault-adjacent patch boundaries, stage 1: NNC-style coupling from parameter-space overlap areas (flow-only).
- [ ] Diagnostics / partial degradation: report offending cells + reason instead of aborting the deck.
- [ ] Overlap-disjointness check for CARFIN boxes (currently unvalidated → silent corruption; Part III.3).
- [ ] Guard: global COMPDAT inside a CARFIN box must throw (currently silently attaches to an arbitrary refined child; Part III.4).

## Track 1 — the rebuilt static refinement (main track)

Milestones, gated on the deck matrix (below):

1. **Strip** — **done** (branch `strip-lgr`, 2026-06-12): removed `LgrHelpers.*`, `refineAndUpdateGrid` + private helpers, `NestedRefinementUtilities`, `CpGridUtilities`, `refineSingleCell`/`refineCellifiedPatch`, mark/adapt plumbing and the LGR test suite (~17.8k lines). Kept the public DUNE facade as inert/throwing stubs, the multilevel data model, and `LgrOutputHelpers`. Verified: all 39 remaining opm-grid tests pass (serial+parallel); *unmodified* opm-simulators builds (superbuild variant in `superbuild-refined/`, build dir `builds/refined`); SPE9 runs serially and with `mpirun -np 2`; a CARFIN deck aborts with a clear "refinement removed" message.
   Discovery: mainline flow lists CARFIN/ENDFIN as *critical unsupported keywords* (`UnsupportedFlowKeywords.cpp`) — CARFIN decks only reach the grid machinery with `--parsing-strictness=low`, even on upstream master with full LGR support. Note for Track 0/upstreaming: enabling CARFIN by default upstream is itself a pending step.
2. **Seam** — **in progress** (2026-06-12): `Opm::Refinement::Builder` interface + `BlockRefinement` request type landed (`opm/grid/cpgrid/refinement/`); `CpGrid::addLgrsUpdateLeafView` now validates (per-box consistency **and pairwise disjointness** — the Part III.3 silent-corruption fix, a hard error here) and dispatches to the registered builder, throwing a clear message when none is registered. Covered by `refinement_seam_test` (40/40 tests green). The milestone-2 remainder (state model and construction-stable ids) is now specified in [DESIGN-builder.md](DESIGN-builder.md): shared corner pool across levels and leaf (D1, kills the corner-identification layer and `corner_history_`), materialized leaf topology with shared corner geometry (D2), packed 64-bit construction-stable ids (D3), and the staged S2 pipeline with per-stage property tests (D4–D5).
3. **Builder, conforming core** — **in progress** (2026-06-12): pipeline stages 2+3 landed as `refinement/GrdeclRefinement.{hpp,cpp}` — pure-function resampling of COORD/ZCORN/ACTNUM onto sub-pillars (corner-point-native semantics, see DESIGN-builder D4 note). Property-tested in `grdecl_refinement_test`: uniform-lattice exactness, volume conservation on distorted vertical-pillar grids, **faults inside the block preserved**, inactive inheritance. Also landed: stages 4+6a as `refinement/LevelGridAssembler.{hpp,cpp}` + `GridStateWriter` (the single befriended writer of multilevel state) — one block request becomes a fully processed refined *level grid* (preprocessor run on the resampled block, so faults inside the block are matched by `findconnections`), with parent relations driving `father()`/`geometryInFather()` and CARFIN-local Cartesian indices for free via `global_cell_`. Property-tested in `level_grid_assembler_test` (volume partition per parent incl. faulted blocks, inactive parents childless, reference volumes in father). **Leaf assembly landed** (stages 5+6b, `refinement/LeafGridAssembler.{hpp,cpp}`): leaf merges level 0 and the refined grids — children replace parents in the cell ordering, block-boundary parent faces are replaced by the refined mosaic (paired with the coarse neighbor), corners follow the D1 shared pool (level-0 prefix preserved; identification through parent cells so faults inside blocks stay consistent), and ids work through the existing delegation (`leaf_to_level_cells_` + `corner_history_` populated — no insertIdSet machinery beyond registering the new views). Boundary conformity is checked hard: fault-split or partial parent faces at the block boundary throw (incl. the diagonal-neighbor case where the cartesian neighbor count alone would pass). **First backend registered**: `refinement/ConformingBlockBuilder` (serial, GLOBAL parents, separated boxes; faults *inside* blocks supported) — `CpGrid::addLgrsUpdateLeafView` works end-to-end again. `conforming_builder_test`: leaf volume conservation, parent/child structure via the facade, two-sided intersection symmetry over the whole leaf (incl. across a refined interior fault), mosaic neighbor sees 9 intersections, unique global ids for cells and points, guards for touching boxes and faulted boundaries. **Structural-equivalence oracle tests added** (`refined_structure_comparison_test`, mirroring upstream's LgrChecks methodology): (a) a whole-grid box on a *faulted* grid produces a leaf structurally identical to the directly processed refined description — cell/point counts, interior face pairs, per-cell volume and centroid at 1e-10, including along the refined fault; (b) for a *local* box containing a fault, the refined region matches the direct block oracle cell-by-cell and in interior connectivity. 44/44 tests green.
   **Flow integration landed (2026-06-12): CARFIN decks run end-to-end through unmodified flow.** Two pieces closed the gap: `getParentIntersectionFromLgrBoundaryFace` restored (upstream implementation, only surviving APIs — transmissibility needs it at LGR boundaries), and the post-MINPV corner-point description is retained in `CpGridData` when the deck has LGRs (`RetainedCornerPointInput`, captured in `processEclipseFormat`), so `addLgrsUpdateLeafView` self-constructs the `ConformingBlockBuilder` — no opm-simulators change, no Vanguard hook needed.
   **Deck matrix, opm-tests/lgr, serial, `--parsing-strictness=low`** (rebuilt flow in `builds/refined`): CARFIN1 ✅ (53 Newton its vs upstream's 51 — quantitative comparison TODO), CARFIN_GR ✅, 2DCORNERPOINT_XY ✅, 3DCORNERPOINT_XYZ ✅; CARFIN and CARFIN_FLEX ❌ with the loud touching-boxes restriction (**upstream runs these** — known feature gap, needs mosaic-mosaic identification); CARFIN_FAULTS and XYZ-NON ❌ at EclipseState parsing **on upstream master too** (not a refinement issue).
   **Edge/corner-sharing boxes landed (2026-06-12):** the leaf assembler merges refined corners shared between touching boxes through an exact-coordinate pool — sound because equal subdivisions on the shared boundary make the resampling arithmetic identical, so shared corners coincide *bitwise* (not a floating-point heuristic). The builder classifies box pairs (separated / edge-or-corner / face) and requires equal factors for touching boxes; **face-sharing is still rejected** (mosaic-mosaic pairing pending). The real `CARFIN` deck (diagonal LGR1/LGR2) now runs end-to-end and is **bitwise-identical to upstream**: PORV and TRANX/Y/Z in all three sections, and all 2x216 boundary TRANGL connections (connection-keyed). `edgeSharingBoxesMergeCorners` test asserts no coincident-distinct vertices; `faceSharingBoxesThrow` guards the unsupported case. 44/44 tests green.
   **Face-sharing boxes landed (2026-06-12):** boundary faces of two refined blocks on a shared plane pair into interior faces, matched by their merged corner set (corners already merge, so paired faces have identical sets); the existing mosaic-outside mechanism supplies the second cell (refined or coarse, uniformly). A box's boundary faces may be *mixed* — some toward another refined box, some toward coarse — handled per parent. The builder enforces matching subdivisions only in the *overlap* (shared) dimensions and only for pairs that actually interact (no gap in any dimension); out-of-face factors are free. **All CARFIN decks in opm-tests/lgr now run** (CARFIN, CARFIN1, CARFIN_GR, CARFIN_FLEX, 2D/3D cornerpoint) and are **bitwise-identical to upstream**, including CARFIN_FLEX (6 LGRs, mixed separated/edge/face adjacencies): PORV+TRANX/Y/Z across all 7 sections and all 580 boundary TRANGL connections, 0.0 diff. Tests `faceSharingBoxes` (A/B blocks connected, vertex merge, ids) and `faceSharingNonMatchingSubdivisionsThrow`. 44/44 green. The only opm-tests/lgr decks that still fail (CARFIN_FAULTS, XYZ-NON) fail at EclipseState parsing on upstream master too.
   **Faults AT the box boundary — NEXT (planned, 2026-06-15).** Faults *inside* a box work; a fault *on* the box boundary throws (`LeafGridAssembler::outsideNeighborOf`: "a fault-split / partial face at the block boundary"). Repro: `opm-tests/flow_diagnostic_test/SIMPLE_2PH_W_FAULT_LGR.DATA` (box `CENTER`, geometric corner-point fault through parent cell 1347). Upstream master also fails (`LgrHelpers.cpp:507`, ">6 faces", hexahedral-only) — the fork's polyhedral leaf is the right substrate (`face_to_point_` is a variable-length `SparseTable`; coarse box-neighbours already carry >6 faces). **Approach (decided, conformal split-faces):** at a faulted boundary, build the connecting faces with the existing corner-point face processing rather than new clipping. Integration **confirmed: use `process_grdecl` on a "mini-grdecl"** (the box-boundary region + a coarse-neighbour shell, refined to the box factors), then read split connections from its `processed_grid` (`face_neighbors`/`face_nodes`/`face_tag`/`local_cell_index`) and map shell sub-cells back to the coarse neighbour leaf cells — NOT `findconnections` standalone (it is wired into preprocess.c's unique-point/accumulator internals, impractical to call directly). Edge-conformal via the `edge_conformal` flag + the existing `edgeConformalizeLeaf(*leaf)` pass. Full design + phases in the approved plan file. Scope: rank-interior (serial / box-on-one-rank); independent of the D3/output track.
   **Build phasing (decided 2026-06-15):** (1) get the **topology and true polyhedral geometry** right — the split faces from `process_grdecl` (with the boundary corners inserted on the shared coarse pillars, so the result is edge-conformal by construction), their new vertices added to the leaf corner pool, and *real* polygon geometry (centroid / area / unit normal) for the synthetic faces; verify structurally (builds, volume conserved, two-sided intersections, no regression on clean boundaries). (2) **Separately, on top**, apply the **ECLIPSE-convention geometry** that transmissibility needs — the face centres / distance vectors used by `faceCenterEcl` / `distanceVector` in `Transmissibility` are ECLIPSE-style (pillar-projected), not true centroids, so correct flow behaviour at these faces is a second layer added after the real geometry is in place. Do **not** conflate the two: phase 1 is "a valid polyhedral grid", phase 2 is "ECLIPSE transmissibility behaviour".
   Implementation status: **phase 1 landed.** `faultedBoundaryConnections` (commit ea26c45c) computes the split connections; `LeafGridAssembler` now detects faulted box sides (pre-pass), suppresses the box's faces there, and rebuilds them from the processor — split faces to coarse neighbours plus domain (fault-scarp) boundary faces — inserting the new split vertices into the corner pool with real polygon geometry. Faulted-boundary decks now **build** a structurally valid leaf (volume conserved, every cell closed, edge-conformal): `conforming_builder_test` `faultAtBoxBoundary*Builds`, `faulted_boundary_test`, `edge_conformal_refinement_test::refineBoxBoundaryOnFaultBuilds`. `SIMPLE_2PH_W_FAULT_LGR.DATA` builds the grid (the "partial face" throw is gone). **Phase 2 — trans fix landed (commit 91128d7f):** the solver NaN was `faceCenterEcl` on the non-planar (skewed) fault faces — its fixed four-vertex average collapsed onto a cell centre (zero distance). Now planar quads keep the ECLIPSE four-vertex average (clean LGR unchanged, CARFIN1 bit-for-bit at 52 Newton) and non-planar / many-vertex faces use the true centroid. **The case now simulates end-to-end (36/36 report steps, exit 0).** Two follow-ups remain: (a) the `RPTSCHED FIP` fluid-in-place *report* crashes for the LGR leaf (`Inplace::get(FIELD)` in `LogOutputHelper::fip`) — an output/reporting issue, not the solve (the run completes with FIP suppressed); (b) a detailed comparison against the ECLIPSE reference in `opm-tests/flow_diagnostic_test/eclipse-simulation/` (a tracer deck with a non-standard summary set, so the keyword sets must be reconciled first).
   Remaining for milestone 3: quantitative output comparison vs upstream (done for trans/PORV — bitwise; restart fields at Newton tol), and **porting the applicable upstream test invariants** — upstream `tests/cpgrid/lgr/LgrChecks.hpp` (checkEqualCellGeometrySet, checkEqualIntersectionsGeometry, checkFatherAndSiblings, id-consistency checks, …) is the systematic acceptance harness to adapt; several upstream tests (global_refine, lgr_cartesian_idx, lgrIJK semantics) apply conceptually to the new implementation and should be revived once nested/parallel features catch up. *Related work to review before finalizing*: JutulDarcy (SINTEF's Julia simulator) has implemented refinement on corner-point grids — compare its approach (and its handling of fault/pinch-out cases).
4. **General corner-point**: pinched parents; fault-adjacent boundaries via `findconnections` + lateral splitting at sub-pillar positions (§8.2-S2); edge-conformal post-pass option for VEM (§3).
5. **Wells-in-LGR** (WELSPECL/COMPDATL through the existing opm-common path) + ECL output.
6. **Parallel** — design settled (2026-06-12); LGR-aware partitioning API landed, with a key finding about the underlying scatter:

   **Landed:** `CpGrid::setPartitionCellGroups(vector<set<int>>)` — groups of Cartesian cell ids the partitioner must keep on one rank. The GraphOfGrid path (`addPartitionCellGroups`, the `zoltanGoG` method LGR requires) contracts each group into one vertex, the same `addWell` mechanism used for well cells. No-op when unset, so non-LGR runs are unaffected (verified: SPE9 np=2 unchanged at 233 Newton its; 46/46 grid tests incl. a 4-proc test). The intent: keep each refinement box on one rank so the builder never splits a box.

   **First-implementation model (accepted): LGRs are rank-interior.** Each box (with its halo) lives entirely inside one rank's interior — enforced by the cell-group contraction above. Then parallel refinement is embarrassingly simple: each rank refines only its own boxes, the rest of the grid stays coarse, refined cells are interior-only (no cross-rank refined entities), and the existing coarse-cell overlap/index sets are untouched. No refined-entity communication, no id prediction.

   **Overlap finding (acted on the "add corners to overlap" hint):** contracting a *fully-interior* region with overlap layer **1** trips the `computeFace2Cell` "oneValid" scatter assertion — CpGrid's overlap-layer-1 gap ("cells sharing only corners/edges, not faces, with interior cells are not seen"). **Overlap layer 2 fixes it** (it captures the corner/edge neighbors): with overlap 2, a contracted interior box stays whole on one rank on **2 ranks** (verified, `partition_cell_groups_test::boxStaysWholeOnOneRank`). So parallel LGR should use overlap ≥2. **Still open:** on **≥4 ranks** the contracted-region scatter still asserts at multi-rank junctions regardless of box placement — a deeper CpGrid distribution limitation (the face distribution around contracted regions) that needs dedicated work before >2-rank LGR runs. **Next:** the distributed builder (refine local rank-interior boxes; refined cells appended as interior with per-rank-disjoint ids; coarse overlap unchanged), tested first on 2 ranks via the working overlap-2 partitioning; then the ≥4-rank scatter fix.

   Remaining mechanics after that (unchanged from the design):
   - **LGR-aware partitioning first, general redistribution later**: restrict the level-zero partition so no refinement box is split across ranks — implementable by merging each box's cells into one vertex in the partitioning graph, exactly the mechanism `GraphOfGrid` already uses for well cells. Requires the boxes at loadbalance time (flow calls loadBalance before addLgrs): an additive, backward-compatible `loadBalance` parameter (cell-group list), wired from the app-repo Vanguard or a small upstream PR. Redistribution *after* refinement stays unsupported (same as upstream) until a general solution.
   - **Communication-free refined geometry**: the retained corner-point input is global on every rank and refinement is deterministic, so any rank can construct the children of any parent it sees (owned or ghost) locally — no geometry/topology communication. Each rank refines every parent in its owned+overlap region.
   - children inherit the parent's partition type (assembler extension); leaf/level `cell_index_set_` built from the construction-stable ids (D3 — rank-independent, so no sync step); communication interfaces rebuilt from ids by the existing machinery; 2-rank CARFIN deck as the acceptance test.

   **Distributed builder infrastructure landed (2026-06-13):** `ConformingBlockBuilder::build` now has a distributed path. Per rank it classifies each box as *owned* (all cells local + interior) or *absent* (no cell local), throwing if a box is split or touches the overlap (the rank-interior enforcement — `distributed_builder_test::boxTouchingOverlapThrows`, passes np=2/4). The owning rank refines the box; other ranks get an empty placeholder level grid (`assembleEmptyLevelGrid`) so the level count matches across ranks. Refined level grids are built with a **self-communicator** (`MPIHelper::getLocalCommunicator()`), which was essential: building them with the world communicator made `processEclipseFormat` collective and deadlocked when only the owning rank refined. Serial is unaffected (CARFIN decks unchanged; 48/48 grid tests).

   **Distributed leaf assembly now works (2026-06-13).** The fault was the leaf assembler indexing `boxOfCell[]`/`leafIdxOfCell0[]` with `std::numeric_limits<int>::max()` — the off-rank sentinel that a distributed `face_to_cell` uses for a neighbor on another rank. Fixes: skip the sentinel in the face-drop and `outsideNeighborOf` scans; preserve it in the leaf `face_to_cell` rows; build the leaf `cell_to_face` manually (skipping the sentinel) since `makeInverseRelation` sizes its table from the max cell index and the sentinel blows it up. The leaf `cell_indicator_` partition types are now set (coarse cells from level zero, refined cells from their parent — via a new befriended `GridStateWriter` accessor on `PartitionTypeIndicator`), so overlap cells are no longer mis-counted as interior. `rankInteriorBoxRefinedInParallel` is **enabled and passes on np=2 and np=4**: exactly one rank refines, 144 children, and the interior volume is conserved across the distributed leaf. 48/48 grid tests; serial unchanged (CARFIN1 at 53 Newton its).
   **Flow wiring done (2026-06-13, opm-simulators branch `lgr-partition-cell-groups`).** `CpGridVanguard::loadBalance` now keeps each CARFIN region on one rank: it builds one cell group per connected set of boxes (union-find; touching/overlapping boxes merged) **plus a 2-cell halo** (so the box never lands in a neighbour's overlap), calls `setPartitionCellGroups`, and forces overlap-2 + `zoltanGoG`. Guarded by a compile-time trait so it compiles against upstream opm-grid too. The fork also skips the original's global-view id-sync (the rank-interior builder doesn't use it). Grid-side, `addLgrsUpdateLeafView` now **broadcasts the retained corner-point input** from rank 0 (it lived only there) so every rank's builder can resample the global box geometry.
   With these, **parallel CARFIN refinement runs on all ranks** (both ranks refine their rank-interior boxes, no "split across ranks").

   **Parallel leaf `cell_index_set_` + interfaces done (2026-06-14, grid commit `Build parallel cell index set + interfaces for the refined leaf`).** `assembleLeafGrid` now builds the distributed leaf's parallel structures (coarse cells inherit level-zero global id + attribute; refined cells are rank-interior owners with fresh per-rank-disjoint global ids; then `cellRemoteIndices().rebuild`, `computeCellPartitionType/PointPartitionType/CommunicationInterfaces`). This fixed the `dofVolume>0` assertion (overlap cells were stuck at zero volume because the leaf had no communication interfaces for opm-models' border sync). Also: a `CollectDataOnIORank::collect` guard in opm-simulators so parallel-LGR output degrades instead of segfaulting.

   **State: `SPE1CASE1_CARFIN1` runs to completion on 2 ranks and converges.** Correctness probe (partition-invariant sum over interior cells): identical initial condition and **bit-identical global mass-balance residual** vs serial — global conservation is correct. **Two remaining items for fully-correct parallel LGR:**
   1. A **localized refined-region transmissibility discrepancy**: local CNV residual ~27× larger at iteration 0 in parallel while MB is identical — a *mass-conserving* flux misrouting between adjacent refined cells (tightening tolerances does not remove the ~1e-4 solution difference; the non-LGR control matches serial↔parallel to ~1e-7). Next: dump per-face transmissibilities serial vs parallel over the refined region, trace the differing face to geometry / partition-type / NNC handling.
   2. **Parallel-LGR output collection** (`CollectDataOnIORank`): upstream-unimplemented; needs the global (equil/IO) grid refined consistently and the gather keyed on **construction-stable global ids (D3)** rather than cartesian indices (which collide for refined siblings) — the same id work the octree makes load-bearing.

   Serial and the static path are unaffected (CARFIN1 53 Newton its; 49/49 grid tests, incl. the now-enabled `rankInteriorBoxRefinedInParallel` np=2/4).

7. **Verification vs upstream — full-matrix check (`scripts/compare_lgr_output.py`)**: two comparisons. (a) **Full element-wise, every INIT array and restart field, per LGR section** — valid *because the cell numbering is identical*, since INIT/restart are written in per-section Cartesian (active-cell) order; identical numbering makes a direct element-wise comparison meaningful and bitwise-identical INIT then proves grid/trans/PORV/depth are identical. (b) **Connection-keyed TRANGL**, keyed by (NNCG, NNCL) cell pairs from the EGRID, immune to face-numbering differences.
   Result across the passing deck matrix (CARFIN, CARFIN1, CARFIN_GR, CARFIN_FLEX): **every INIT array bitwise-identical** (worst-rel 0.0; e.g. CARFIN_FLEX 152/152 arrays over 7 sections) and **every boundary TRANGL connection bitwise-identical** (e.g. CARFIN_FLEX 580/580, no only-A/only-B). Restart differs only at solver tolerance: on CARFIN_FLEX, peak |dP| 15.5 bar at step 1 (injection startup, 0.26% of ~5870 bar) decaying monotonically to ~0.7 bar (~0.01%) — the signature of linear-solver-path roundoff from internal leaf cell ordering, not discretization (a real difference would persist and INIT would not be bitwise-identical). SGAS's large *relative* diffs are near-zero-saturation artifacts (peak abs 0.04 at the moving gas front, also decaying). Newton-count differences (52-vs-49, 100-vs-99) are the same ordering effect.
   *Optional future check*: reproduce upstream's exact internal leaf ordering → bitwise-identical restart and matching Newton counts. Not needed for correctness (the above already proves identical discretization), but would make the equivalence trivially demonstrable.

## Edge-conformal grids (serial) — works; no fixes needed

`edge_conformal_refinement_test` (2026-06-13): refinement on grids built
with `edge_conformal=true` works serially — a box away from a fault, and a
box *containing* a fault, both conserve volume; a box boundary on a fault
throws (the general faulted-boundary restriction, not edge-conformal
specific). The **four-columns-faulted-around-a-pillar** torture grid is
covered: `edge_conformal` demonstrably inserts the neighbour columns' nodes
(total face-node count grows vs the plain grid), and refining the flat
corner away from the complex pillar conserves volume — the edge-conformal
coarse faces survive refinement. **Serial needs no fixes** for these cases.

**Edge-conformal leaf — implemented as a toggle (2026-06-13).**
`refinement/EdgeConformal.{hpp,cpp}::edgeConformalizeLeaf` makes the refined
leaf edge-conformal: for every face edge it inserts the leaf nodes lying on
that edge's interior into the face's node list (a spatial-hash sweep — the
"know all nodes on a pillar/edge first" step). This closes the hanging-node
gap that refinement opens on the LGR-boundary pillars and lateral edges for
diagonal coarse neighbours. Only `face_to_point_` changes — no new corners,
so cell geometry, face area and normal are untouched (verified: identical
volume on/off). The toggle is a `ConformingBlockBuilder` constructor flag,
**default off (current face-conformal behaviour)**; the deck path turns it
on automatically when level zero was built `edge_conformal=true` (carried on
`RetainedCornerPointInput`). Covered by
`edge_conformal_refinement_test::edgeConformalToggleInsertsBoundaryNodes`.
Not needed for flow/TPFA (face centroid/area/normal suffice); intended for
VEM. Note: the refined level grids themselves remain `edge_conformal=false`
internally — the leaf post-pass makes the *assembled* leaf conformal, which
is what consumers see.

## Track 2 — construction-stable ids (D3) and parallel-LGR output

Status of parallel LGR (2026-06): the rank-interior model works and is correct
where it overlaps upstream — Cartesian and non-faulted corner-point decks agree
with master to 1e-6, serial and parallel (see `scripts/lgr_regression.sh`). Two
gaps remain, with a shared root, and one explicit non-goal.

**Key correctness note — do not "fix" the solve path.** In the rank-interior
model a refined cell is interior-only and never appears in another rank's
overlap, so the per-rank-disjoint parallel `cell_index_set_` ids
(`LeafGridAssembler.cpp`, the `refinedNext` block) only need *local* uniqueness
and are correct for the solve. D3 is **not** about that `int32` parallel index.

**D3 — construction-stable global id (the keystone).** Give each refined leaf
cell an id derivable from the input alone, so it is identical regardless of
which rank (or run) built it: `encode(parent stable id, child index)`. The
inputs already exist at leaf assembly time — `leafChildToParent` (parent
level-0 cell → stable `globalCell()` Cartesian) and `leafIdxInParent` (child
index within the parent). The Dune entity `IdType` is `int64_t`, so the packed
id has room. This is what the output gather keys on; it is also the prerequisite
for any future distributed refinement (matching a box's overlap cells across
ranks) and for a forest/octree core.

**Output collection (the immediate payoff).** `CollectDataOnIORank` currently
bails for parallel + LGR (index maps left empty; we degrade to no gathered
cell output — VTK is the interim cell output, which works in parallel). To
restore summary-correct, cell-correct output: refine the IO/equil grid
consistently (today it stays coarse, 300 cells) and gather distributed refined
cells into the per-LGR-section output order keyed on the D3 id (Cartesian keys
collide for refined siblings). Well **summary** rates are already correct in
parallel; this is about the gridded restart/INIT cell arrays. Acceptance: the
`compare_lgr_output.py` INIT/restart element-wise check passes parallel-vs-serial.

**Distributed refinement — only if cheap; otherwise out of scope.** Letting a
box span ranks (each rank refines its part) would close the whole-grid-LGR gap
(`CARFIN_GR` in parallel, which master does and we now reject cleanly). Decision
(2026-06): keep the rank-interior model; the partition is adjustable
(`setPartitionCellGroups`), which is sufficient for most decks. Pursue
distributed refinement only if it falls out of D3 without much added
complexity. D3 + output do **not** depend on it.

## Track 3 — AdaptiveCpGrid (dynamic, parallel local refinement)

Promoted from "out of scope" to an **active design track** (2026-06-23). This is
the dynamic successor to the static builder: a parallel, locally/dynamically
refinable corner-point grid. Full design: [DESIGN-parallel-octree.md](DESIGN-parallel-octree.md)
(§10–§14 added this session); representation/id details in
[DESIGN-builder.md](DESIGN-builder.md) (D1–D3); redistribution context in
[REDISTRIBUTION-requirements.md](REDISTRIBUTION-requirements.md) /
[REDISTRIBUTION-status.md](REDISTRIBUTION-status.md). Work lives on branch
**`adaptive-cpgrid`**.

**Name & constraint.** The class is **`AdaptiveCpGrid`** — new and **additive**:
it must not modify or regress `CpGrid` or the static LGR path (those tests stay
green); it reuses the shared kernels (`GrdeclRefinement`, `assembleLeafGrid`,
`edgeConformalizeLeaf`) unchanged.

**Decisions fixed 2026-06-23:**
- **Scope** (octree §2): always refined from a corner-point grid, **never
  general polyhedral**. Refinable parent = clean hex (8 corners); leaf cells are
  hex *geometry* but may carry **>6 faces** (conformal sub-face mosaics).
  Interfaces stay cheap — a 2-D clip in pillar parameter space (review §8.1).
- **Representation** (octree §10): two layers. *Layer A* (persistent / migrated /
  serialized, compact) = the macro corner-point input (`RetainedCornerPointInput`)
  + a forest of per-root trees. *Layer B* (derived, per-rank, owned+ghost only) =
  the materialized leaf `CpGridData`; multi-face cells need **no new container**
  (`cell_to_face_` is already a variable-length `SparseTable`). No full global
  vectors.
- **Fast incremental refinement** (octree §12): mutate the leaf vectors locally —
  append + free-list + D3-id identity + stencil-only edits + lazy compaction.
  Full-leaf rebuild is only the correctness oracle. Serial-first.
- **Split policy** (octree §13): arbitrary anisotropic **first** split per root
  (CARFIN-exact), then **factor-2 anisotropic** per dynamic level (Morton-clean
  ids, textbook 2:1 balance, bounded ≤2×2 interface mosaic).
- **AMR backend** (octree §14): **deferred, gated on a spike.** p4est/t8code are
  ruled out under the anisotropic policy (they hard-code isotropic 1→8); in-house
  forest is the working assumption. Reuse Zoltan/ParMETIS (`zoltanGoG`) + a small
  SFC + Dune comm + the existing leaf kernels.

**Incremental path** (octree §9, gating the backend decision):
1. **Serial in-house forest + leaf rebuild** — per-cell trees over level zero;
   mark/balance/apply; rebuild leaf via the generalised assembler; coarsening.
   *Gate:* a uniform-depth forest reproduces a CARFIN refinement **bitwise**.
2. **Packed ids (D3)** replacing delegation — serial first (must stay
   bitwise-equal), then they enable parallel and fast incremental adapt.
3. **Parallel** — SFC root partition + ghost-tree layer; owned+ghost leaf build;
   ids-based index sets. 2-rank then 4-rank (the cases CpGrid's own scatter
   failed).
4. **Dynamic loop in flow** — `adaptGrid` hook + conservative state transfer.

**Redistribution note.** AdaptiveCpGrid's root-tree migration *is* the
redistribution capability CpGrid lacks (move tree bytes, rebuild leaf locally).
A correct-but-slow interim alternative for plain CpGrid (gather to root,
re-scatter) is documented in REDISTRIBUTION-requirements.md but is not on this
track's critical path.

## Acceptance gate — deck matrix (Part II)

CI matrix from day one: unfaulted Cartesian / faulted / pinched / NNC / aquifer × serial / parallel × wells in/out of LGR; ECLIPSE reference output where available. Each milestone and Track 0 item is "done" when its rows pass.

## Out of scope (recorded, not planned)

A polyhedral-first new grid core (Route B, §IV.2) — revisit only if
VEM-on-general-grids becomes a firm goal.

(Dynamic refinement, previously listed here, was **promoted to Track 3**
(`AdaptiveCpGrid`) on 2026-06-23 — design active, implementation pending the
serial spike.)
