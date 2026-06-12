# Review of the LGR work (opm-grid / opm-common / opm-simulators)

State reviewed: master checkouts as of 2026-06-12 in `~/Documents/OPM/opm_vscode_clean`.
Last updated: 2026-06-12 (added §8, alternatives for full general refinement; Part II, way forward for static LGR on current input; Part III, follow-ups; Part IV, repository strategy after lgr_refactor).
Canonical copy: `docs/lgr_review.md` in github.com/hnil/opm-gridrefined. Development roadmap: `docs/PLAN.md`.

## Executive summary

- The LGR implementation (mainly Antonella Ritorto, ~240 commits since 2022) is functionally broad and well-tested for its intended use: **static** CARFIN/AUTOREF refinement at initialization, including nested LGRs, parallel runs, wells-in-LGR and ECL-compatible output. The mark-based kernel is general; CARFIN is a wrapper (§1, §5).
- It is restricted to clean hexahedral parents: no fault-split cells (>6 faces), no pinched cells, no NNCs, no coarsening, no growth of refined regions, no repartitioning after refinement (§1).
- **General corner-point support** fails at the LGR *boundary*, not the interior: entity identification is index arithmetic that cannot express partial face overlap. Proper support needs geometric matching — months of work, not a patch (§2). The most robust route is to do refinement at the *preprocessor* level, reusing `findconnections()` which already handles non-matching/fault stacks (§8.2-S2); all required intersections are 2D polygon clipping in pillar-pair parameter space — no CGAL needed (§8.1, §8.3).
- **VEM**: the refined leaf is face-conformal but not edge-conformal; a post-pass analogous to `make_edge_conformal` makes it VEM-usable. LGR on `edge_conformal`-constructed grids is currently untested/likely broken (§3).
- The implementation is fundamentally **static** (full leaf rebuild per call, O(global) cost); dynamic refinement should reuse the level/parent-child data model but replace the build algorithm — best as a forest-of-octrees state over the level-0 macro-grid with 2:1 balance, Morton ids and refinement-invariant cached fault overlaps (§4, §7, §8.4).
- The `CpGridLGR` subclass split (refactor branches) moves declarations, not coupling. Separate the code by **composition** instead: a LevelHierarchy state component plus a RefinementBuilder interface, keeping one `CpGrid` type and creating the seam for alternative backends — the current static builder, a general corner-point builder, a dynamic tree-based builder (§6, §8.5).
- **Way forward after the rejected lgr_refactor** (Part IV): for the static-LGR goal, no complete new grid is needed — `CpGridData` can already represent the result; the blockers live in the build algorithm. Develop in **opm-gridrefined** (opm-grid history, LGR layer stripped, new builder component), keep opm-simulators untouched, and host the experimental executable/Vanguard in a separate thin app repo later (opm-geomech pattern).

## 1. What is there — architecture and scope

Main author of the core work is Antonella Ritorto (~240 commits since 2022 on the central files).

**Code layout.** `CpGrid.hpp`/`cpgrid/CpGrid.cpp` hold the public entry points: `addLgrsUpdateLeafView()` (CARFIN-style block patches, incl. nested LGRs via parent-grid names), `autoRefine()` (AUTOREF), `globalRefine()`, and the DUNE adaptivity interface `mark()/preAdapt()/adapt()/postAdapt()`. All funnel into one ~540-line orchestrator, `CpGrid::refineAndUpdateGrid()`. The bulk of the algorithm lives in `cpgrid/LgrHelpers.cpp` (~2300 lines of free functions in `Opm::Lgr`), plus `NestedRefinementUtilities`, `LgrOutputHelpers`, `CpGridData::refineSingleCell()` and `Geometry<3,3>::refineCellifiedPatch()`. 33 dedicated test files under `tests/cpgrid/lgr`, with substantial parallel coverage.

**Data model.** `CpGrid::data_` is a vector of `CpGridData`: `[0]` level zero, `[1..maxLevel]` refined level grids, `back()` the leaf view. The leaf is fully materialized — own geometry, topology and index sets; coarse-cell geometry is *duplicated* into the leaf. Relations maintained: `child_to_parent_cells_`, `cell_to_idxInParentCell_` (powers `geometryInFather()`), `level_to_leaf_cells_`/`leaf_to_level_cells_`, `corner_history_`, per-level global ids with cross-rank synchronization, `LevelCartesianIndexMapper`.

**Mechanism.** Every marked cell is refined *individually* into a throwaway single-cell `CpGridData` LGR using the trilinear map of its 8 corners (`refineSingleCell` → `refineCellifiedPatch` with patch `{1,1,1}`). Shared refined corners/faces between neighboring marked cells are identified by closed-form index arithmetic (`replaceLgr1CornerIdxByLgr2CornerIdx`, `getParentFaceWhereNewRefinedCornerLiesOn`, …). The leaf is rebuilt from scratch on every call; parent faces on the LGR boundary are replaced by the mosaic of refined faces, so the coarse face-neighbor becomes a polyhedral cell with >6 faces in `cell_to_face_` (its `cell_to_point_` stays 8 corners).

**Simulator integration** is real and fairly deep: `CpGridVanguard::addLgrs()` adds LGRs after load balancing (and again on the global view, then `syncDistributedGlobalCellIds()`); transmissibility handles LGR-boundary intersections via `getParentIntersectionFromLgrBoundaryFace()`; wells inside LGRs work through the stack (`is_lgr_well`, per-connection grid ids); output writes EGRID LGR sections (`EclipseGridLGR`), LGR restart headers (`LgrHEAD*`), per-level restart values. opm-common parses CARFIN/ENDFIN only (`LgrCollection`).

**Strengths.** Functionally broad (nested LGRs, LGRs sharing faces with equal subdivision, inactive parent cells, distributed grids, wells-in-LGR, ECL-compatible I/O incl. restart); honest DUNE adaptivity interface; unsupported combinations throw instead of silently producing wrong results; good test discipline.

**Weaknesses (code level).** `refineAndUpdateGrid` is a monolith feeding 15+ helpers with 10–20 parameters each through ~20 `std::map<std::array<int,2>,…>` relation maps — O(N_global·log N) with large constants on *every* call, plus one heap-allocated grid object per marked cell. The index-arithmetic identification hard-codes tensor-product/hexahedral assumptions in dozens of places — the single biggest obstacle to generalization. The leaf duplicates all geometry (≈2× memory plus levels).

**Restrictions enforced today.**
- Parent cells must have exactly 8 distinct corners (no pinched/degenerate cells) and at most 6 faces (no fault-split cells — "more than six faces … not supported yet").
- No NNCs in refined cells; aquifer cells filtered or throw; `MULTZ` throws with LGRs; `TRANX/Y/Z` edits throw.
- No coarsening (`mark(-1)` not supported).
- A level-0 cell touching an existing LGR boundary cannot be marked → refined regions cannot grow incrementally.
- Neighboring LGRs sharing faces must have matching subdivisions (conforming interfaces only).
- Load balancing: distribute level 0 then add LGRs works (zoltanGoG); redistributing an already-refined grid is unsupported.
- `autoRefine` factors must be odd.

## 2. General corner-point grids — how much work?

**What already works:** arbitrary (distorted) hexahedra anywhere in a corner-point grid — there is an explicit test for non-rectangular faces — as long as the patch *and its boundary* avoid fault-split cells, pinch-outs and NNCs. Interior geometry is the trilinear image of each parent; mosaics on shared parent faces are consistent between the two sides because both restrict to the same bilinear surface.

**The gap** is (a) cells with more than 6 faces (fault neighbors), (b) degenerate cells (<8 distinct corners), (c) NNC connections. The experimental branch `hnil/copilot/add-failed-test-lgr-case` attacks (a) by skipping duplicate (tag, orientation) faces in `refineSingleCell` — but that is a band-aid: the child faces span the whole logical face while each fault-split parent face covers only part of it, so the leaf-face replacement is geometrically wrong across the fault.

**Effort assessment: substantial, not incremental.** The identification layer assumes (i) each parent face is one of six logical faces, (ii) the refined mosaic covers it exactly, (iii) neighbor matching is solvable by index arithmetic. Fault-adjacent LGR boundaries need *geometric* intersection of refined boundary faces against arbitrary neighbor faces — essentially the fault-matching machinery of `processEclipseFormat`/preprocess applied at the LGR boundary, producing partially-overlapping (non-matching) leaf faces. Degenerate parents additionally need collapse-aware refinement templates. Realistic estimate: a few months of focused work for fault-adjacent boundaries via geometric matching; degenerate cells on top of that. A cheaper engineering shortcut — represent fault-side connections of refined cells as NNC-like transmissibilities without proper face geometry — would unblock flow (TPFA only needs trans) but breaks DUNE intersection semantics and is useless for VEM.

## 3. VEM compatibility / edge conformity

Two distinct questions:

**(a) Is the refined leaf edge-conformal? No.** It is *face*-conformal: coarse faces on the LGR boundary are replaced by the refined mosaic, so face-neighbors become proper polyhedra containing all hanging vertices through their faces. But coarse cells touching the LGR only along an edge — and the lateral faces of face-neighbors — keep their original 4-node faces while new refined corners lie in the interior of those faces' edges. Hanging nodes on edges → a nodal VEM space on the leaf is inconsistent there.

**(b) Can it be fixed? Yes, with a contained post-pass, and there is direct precedent in-tree.** `cpgpreprocess/make_edge_conformal.cpp` (added 2025 for the VEM/geomechanics work) does exactly this operation at preprocessing for fault-induced hanging nodes: insert nodes lying on a face's edges into that face's node list. The same pass applied to the leaf after refinement (insert refined boundary corners into edges of unrefined faces containing them) is well-defined and moderate work — the candidate corners are already classified during refinement (`newRefinedCornerLiesOnEdge`, `isRefinedNewBornCornerOnLgrBoundary`). Caveats:
- `cell_to_point_` is hard-wired `std::array<int,8>`; VEM must derive cell vertex sets from `face_to_point_` (the existing VEM code does this for edge-conformal processed grids already).
- LGR on a grid *constructed* with `edge_conformal=true` is untested and will likely break: faces then carry >4 nodes and the mosaic replacement drops the inserted nodes. Conformity must be (re)established **after** refinement, not only before.

Bottom line: LGR + VEM is feasible via an edge-conformalization post-pass on the leaf, provided LGRs stay away from faults until §2 is solved.

## 4. Static vs dynamic in light of the fundamental implementation

The implementation is a static, build-once design:
- every `adapt()`/`addLgrsUpdateLeafView()` call rebuilds the **entire** leaf (cost ∝ global grid size, map-heavy, allocation-heavy);
- levels accumulate per call and persist (the old leaf is overwritten by the first new level grid, then a new leaf is appended);
- there is no coarsening;
- refined regions cannot grow (marking next to an LGR boundary throws) → no moving refinement window;
- no repartitioning after refinement.

For the intended use — CARFIN/AUTOREF at initialization — this is fine; correctness, not speed, is the criterion, and that is what is wired into flow. As a foundation for *dynamic* refinement the cost model is wrong: per-step adaptation would pay O(global) rebuilds, unbounded level growth and no load rebalancing. A dynamic capability would reuse the level/parent-child data model and the DUNE interface, but replace the build algorithm.

## 5. General refinement vs the DATA-file (CARFIN) representation

The kernel is general: `refineAndUpdateGrid()` takes per-cell marks plus per-level subdivision factors; CARFIN is a thin wrapper converting IJK boxes into marks, and `adapt()` exercises the general path with arbitrary marked sets. But CARFIN assumptions leak into the surrounding layers: `logical_cartesian_size_` and local Cartesian indices are only defined for block LGRs (explicit TODO in the code for the general case); the LGR naming / output / well plumbing assumes named block LGRs; `adapt()` hard-codes `{2,2,2}`; subdivision-compatibility checks exist only for block patches. So: machinery general, productization CARFIN-shaped — and the two are not separated in the code, which is precisely the refactoring question of §6.

## 6. The separation attempt — is the split the way to go?

Master already contains a half-split: the `Opm::Lgr` free-function layer (LgrHelpers, ~2300 lines), `LgrOutputHelpers`, `NestedRefinementUtilities` — while `CpGrid.cpp` keeps the orchestrator and `CpGridData` keeps all multilevel state. The branches (`copilot/refactor-lgr-code-to-class`, compilation-fixed in `hnil/lgr_refactor`) move the public LGR/adaptivity API into `CpGridLGR : public CpGrid` (~1400 lines moved; the branch history shows non-LGR methods were accidentally moved and had to be moved back — the operation is mechanical but error-prone).

**Assessment: separation yes, inheritance no.** The subclass moves declarations, not coupling — `CpGridData` still owns all the LGR members, helpers still reach into both classes, and a second grid type would propagate through every simulator template (the Vanguards instantiate `Dune::CpGrid`), forcing templating or casts. It also enables no alternative implementation: the algorithm remains one fixed pipeline.

**Recommended shape: composition.** Keep a single `CpGrid` type and extract
1. a *LevelHierarchy* state component out of `CpGridData` (level vector, parent-child maps, corner history, LGR names), and
2. a *RefinementBuilder* interface, with the current index-arithmetic implementation as the first backend.

`CpGrid` keeps the thin DUNE adaptivity API and forwards. This preserves current behavior, keeps simulator code untouched, and creates the seam where an alternative backend can live — a geometric-matching builder for faulted corner-point grids (§2), or an incremental builder for dynamic refinement (§7). The existing `LgrHelpers` layer is already ~70% of this move.

## 7. What is needed for dynamic refinement

**Pieces that already exist in-tree:**
- Time-loop hook: `FvBaseDiscretization::advanceTimeLevel()` calls `adaptGrid()` under `EnableGridAdaptation`; the dune-fem implementation `fvbasediscretizationfemadapt.hh` (split out by hnil's commit 37b87aedb, 2023-08; machinery originating in Robert Klöfkorn's 2015 ewoms/dune-fem work) performs mark → adapt → `AdaptationManager` restrict/prolong. It requires a dune-fem-adaptive grid (ALUGrid; `flow_blackoil_alugrid` exists behind `BUILD_FLOW_ALU_GRID`). The earlier unmerged demonstrations lived outside the OPM mainline (ewoms-era forks / dr-robertk's repositories) and were never wired to flow's blackoil+wells production path.
- CpGrid side: `mark/adapt` interface; parent-child maps + `cell_to_idxInParentCell_` + `geometryInFather()` provide everything an FV restrict/prolong needs; a `PersistentContainer` specialization exists (though it is leaf-index-backed, i.e. not truly persistent across adapt — transfer must go through ids or the parent maps).

**Missing at grid level:** coarsening; incremental/local leaf update instead of global rebuild; lifting the "no marking at LGR boundaries" restriction — or, pragmatically, a *rebuild-from-level-0* strategy per adaptation event (re-derive the full mark set on the original grid and rebuild), which sidesteps both coarsening and the boundary restriction and fits the current static builder; level compaction so repeated adaptation does not grow `data_` unboundedly; repartitioning of refined grids.

**Missing at simulator level (the larger half):** blackoil state transfer (conservative restriction; prolongation with per-phase/composition consistency); rebuild of transmissibility, element mappers and field properties on the new leaf (`LookUpData` makes property re-evaluation straightforward); well perforation re-indexing when cells split/merge; aquifers and NNCs; ECL output for a time-varying grid (most plausibly: always report on level 0 by restriction — the parent maps support this).

**Pragmatic path:** (1) keep the static builder and implement "re-adapt = rebuild from level 0 with a new mark set + state transfer via parent maps" — acceptable if adaptation events are infrequent (front-tracking with hysteresis on the marks); (2) only if per-event rebuild cost proves prohibitive, invest in an incremental builder behind the §6 interface. The dune-fem/ALUGrid path remains the fastest route to a research demonstrator of the physics under dynamic refinement, but it bypasses what makes flow production-grade (CpGrid parallel machinery, ECL I/O, wells) — which is why it never merged.

## 8. Alternatives for a full, proper CpGrid refinement

This section sketches concrete routes to a refinement capability that handles *general* corner-point grids (faults, pinch-outs), separated into the static and the dynamic part, plus the two cross-cutting questions: how to compute the intersections at refinement boundaries, and whether to exploit the logical Cartesian structure (octrees).

### 8.1 The key geometric observation

Almost every intersection problem that LGR on a corner-point grid produces is **two-dimensional, in the parameter space of a pillar-pair bilinear surface** — not a general 3D problem:

- A refined LGR-boundary face is an axis-aligned rectangle `[i/nx,(i+1)/nx] × [k/nz,(k+1)/nz]` in the reference coordinates of its parent face (rational coordinates with denominator `cells_per_dim`).
- A coarse or fault-split neighbor face lies on the *same* bilinear surface between the same pillar pair; in that parameter space it is a quadrilateral with straight edges.
- Therefore matching refined faces against arbitrary (faulted) neighbor faces is convex-polygon clipping of rectangles against quads in 2D — exactly solvable with snapped rational/integer coordinates, no exact-arithmetic library required.

The only truly 3D cases are NNCs and grids post-processed by `make_edge_conformal` — both better kept as NNC-like connections than solved by general 3D intersection.

### 8.2 Static refinement — three options

**S1 — Evolve the current per-cell machinery (incremental).** Keep `refineSingleCell` and the leaf rebuild, but (a) replace the pairwise index-translation layer (`replaceLgr1CornerIdxByLgr2CornerIdx`, `markedElemAndEquivRefinedCorn_to_corner`, …) with *global refined-entity keys* — hash of (parent cell ijk, local corner/face ijk) — which makes deduplication O(1) and deletes most of the fragile arithmetic; and (b) add the §8.1 parameter-space clipping only at LGR-boundary faces whose parent face does not exactly match the neighbor (fault-split, partial overlap). Degenerate parents need collapse-aware templates on top (subdivide the reference cube, then collapse refined corners that land on collapsed parent edges — the trilinear map already sends them to the right points; the work is in the topology bookkeeping). Moderate-to-large effort, preserves all existing behavior and tests.

**S2 — Do refinement at the preprocessor level (most robust to generality).** The corner-point preprocessor already solves the hard problem: `process_vertical_faces`/`findconnections` ([facetopology.c](opm-grid/opm/grid/cpgpreprocess/facetopology.c)) match arbitrary non-matching z-interval stacks per pillar pair — faults and pinched cells included — and `make_edge_conformal` and `geometry.c` run on its output. Refinement fits this framework naturally:
- *Inside a patch*: insert sub-pillars (interpolated along COORD lines and laterally on the bilinear surfaces) and process refined columns exactly like ordinary columns.
- *Lateral LGR boundary*: a pillar-pair band where one side carries the coarse column and the other the refined columns. In z, `findconnections` handles the non-matching intervals *with the same code path that handles faults today*; laterally, the coarse face is split at sub-pillar parameter positions — exact, since both sides lie on the same bilinear surface.
- *Fault at the LGR boundary*: nothing special — it is just another non-matching stack, which is the preprocessor's home turf.

Afterwards, the level hierarchy and parent-child maps are derived from the known column structure and attached to CpGridData. Global refinement (AUTOREF) is the trivial special case: resample COORD/ZCORN and reprocess the whole deck — faults and degeneracies handled "for free". Caveats: `preprocess.c` is C with dense per-pillar point lists, so an LGR-local application needs either a windowed run (patch + one-cell halo, then stitching) or extending `process_vertical_faces` to laterally-split bands. This option composes directly with `make_edge_conformal` → edge-conformal refined leaf → VEM (§3).

**S3 — Nonconforming leaf (mortar/NNC).** Compute only overlap *areas* at LGR boundaries (parameter-space clipping or dune-grid-glue) and represent the couplings as NNC-like transmissibilities without constructing matched faces. Cheapest route for TPFA flow across faulted LGR boundaries; but it breaks DUNE intersection semantics and is useless for VEM — acceptable only as an interim flow-only feature.

### 8.3 Computing the intersections — preprocessor parts vs CGAL vs lightweight

- **Reuse of the preprocessor** is the natural first choice: `findconnections` *is* the 1D (z-interval) special case of the needed overlay and already deals with faults, pinch-outs and point uniquification; its per-pillar-pair structure is the right decomposition. The extension to lateral subdivision is a per-band 2D overlay rather than a 1D sweep.
- **Hand-rolled 2D clipping** in face parameter space covers everything §8.1 identifies: Sutherland–Hodgman on convex quads/rectangles, with coordinates snapped to the rational lattice defined by `cells_per_dim` and the z-point indices (making the arithmetic exact in integers). This is a few hundred lines, dependency-free, and robust by construction — recommended.
- **CGAL** would only be justified for genuinely 3D Boolean operations (Nef polyhedra, mesh booleans) which the pillar structure makes unnecessary; it is a heavy dependency with exact-arithmetic cost and real maintenance burden (license is GPL-compatible with OPM, so that is not the blocker). If an external library is wanted at all, **Clipper2** (Boost license, robust integer 2D clipping) matches the actual need far better; **dune-grid-glue** is the in-ecosystem option for surface-mesh intersection/projection and fits S3 (mortar) specifically.
- Keep NNCs and `edge_conformal`-processed inputs out of the intersection machinery entirely — represent them as NNC-like leaf connections.

### 8.4 Exploiting the logical Cartesian structure — octrees for the dynamic part

The corner-point grid hands us a perfect macro-structure for a **forest of octrees**: every active level-0 cell is a tree root, the forest connectivity is the level-0 adjacency (including fault faces, already computed by the preprocessor), and an octant's geometry is the trilinear image of its reference box in the root cell — exactly the arithmetic `refineSingleCell` already implements. With anisotropic factors the trees are k-ary rather than octal, which changes nothing structurally.

What this buys over the current design:

- **Refinement state becomes data, not grids**: a per-cell tree (or per-column, if pillar alignment is kept) replaces the accumulating `CpGridData` levels. Coarsening = delete a subtree. Incremental adaptation touches only the affected trees and their neighbors — no global leaf rebuild.
- **Hanging-node management by rule**: adopting a 2:1 balance constraint (adjacent leaves differ by at most one level) bounds the interface cases to a small enumerable set, replacing the open-ended index-translation layer. The current code has no grading rule at all (only "equal subdivision on shared faces"); 2:1 balance is worth adopting in any redesign.
- **Stable ids for free**: global ids as (root cell index, Morton code of octant) — hash-based, persistent under adaptation, no synchronization maps.
- **Faults stay cheap under dynamic adaptation**: the fault overlap polygons computed once at level 0 (by the preprocessor, in pillar parameter space) are *refinement-invariant*. A refined boundary face's overlap is the clip of its parameter rectangle against the stored polygon — local, exact, never re-run globally. This removes the main cost concern of dynamic AMR on faulted grids.
- **Parallelism**: space-filling-curve partitioning over (tree, octant), ghost layers, and load rebalancing after adaptation are exactly what **p4est** provides, battle-tested at scale (deal.II, ForestClaw). Alternatively, a modest in-house tree implementation is viable because the macro connectivity is already known; p4est's value is the parallel machinery, not the trees.

Architecturally this lands precisely on the §6 composition split: CpGrid = static corner-point **macro-grid** (current preprocessor, untouched and already general) + **tree state** + a **leaf-view builder** that assembles the DUNE view by forest traversal. The current static CARFIN path then becomes "initialize trees from IJK boxes", and the CARFIN-specific layers (names, local Cartesian indices, ECL output) sit on top unchanged.

### 8.5 Recommendation

- **Static, general corner-point**: S2 (preprocessor-level refinement) is the most robust target because it inherits fault/pinch-out handling and edge-conformality from code that already works; S1 with global keys + parameter-space clipping is the fallback if preprocessor surgery is judged too invasive.
- **Intersections**: 2D clipping in pillar parameter space with snapped integer coordinates, written in-house or with Clipper2; reuse `findconnections`' decomposition; no CGAL as a core dependency; dune-grid-glue only if a mortar (S3) interim is wanted.
- **Dynamic**: forest-of-octrees state over the level-0 macro-grid with 2:1 balance, Morton-code ids, cached level-0 fault overlaps, and SFC-based rebalancing (p4est or in-house) — implemented as a second `RefinementBuilder` backend behind the §6 composition seam, so the current static implementation remains the production path while the dynamic one matures.

---

# Part II — Recommendation: working static LGR for current input (CARFIN decks)

Goal: flow runs *real* decks containing CARFIN correctly. This part orders the work by what actually blocks such decks, separated into short term (unblock users on the existing implementation) and long term (structural fix).

## What actually blocks real decks today

Real corner-point models almost always contain pinched/degenerate cells, faults, NNCs, aquifers and transmissibility edits. Against that, the current implementation throws on: parents with <8 distinct corners or >6 faces, NNCs in patches, aquifer cells in patches, MULTZ with LGRs, and TRANX/Y/Z edits with LGRs. The refinement *interior* math is not the problem — any distorted hexahedron refines fine. "Working for current input" is therefore mostly about (a) patch **contents and boundaries** tolerating real-model features and (b) a handful of **simulator-side** gaps.

## Short term — unblock decks on the existing implementation

Ordered by ratio of unblocked decks to effort:

1. **Simulator-side gaps first** (independent of any grid work, well-defined): implement MULTZ and TRANX/Y/Z edits on the LGR leaf in `Transmissibility_impl.hpp`, and aquifer-adjacent handling. These throw today even when the grid machinery would be fine. Order weeks each.

2. **Pinched/degenerate cells inside patches.** Two options, in increasing correctness and cost:
   - *(a) Auto-deactivate degenerate parents inside the patch* (children never created; warn, and optionally add the lost pore volume to a vertical neighbor). Pragmatic, days-to-weeks, mirrors how inactive parents are already supported — but changes active-cell counts vs. the unrefined run.
   - *(b) Collapse-aware deduplication*: refine in reference space as today — the trilinear map already places refined corners of collapsed edges coincidently — then merge coincident refined corners and drop zero-area faces. The existing corner-dedup infrastructure can be extended to do this; 1–2 months. This is the correct end state and also what general pinch-outs need.

3. **Fault-adjacent patch boundaries**, in two stages:
   - *Stage 1 (fast, flow-only)*: keep proper matched faces everywhere except fault-split boundary faces, and build those connections as **NNC-like transmissibilities** from overlap areas computed by 2D clipping in pillar parameter space (§8.1). TPFA only needs the trans; wells and output are unaffected. Roughly 1–2 months, and it unblocks the large class of decks where CARFIN boxes touch faults.
   - *Stage 2*: upgrade the same boundaries to true matched leaf faces using the same clipping output (S1(b)), restoring proper DUNE intersection semantics (and the VEM path, §3). Additional 2–3 months.

4. **Diagnostics and partial degradation.** When a patch still hits an unsupported feature, report the specific cells/faces and the reason, and continue where semantics allow, instead of aborting the whole run. Cheap, and the single biggest usability improvement while gaps remain.

5. **Do not spend short-term effort on**: the `CpGridLGR` subclass split (moves declarations, unblocks nothing), octree/dynamic work (§8.4 — different goal), or performance of the static builder (one-time setup cost is acceptable).

## Long term — structural fix

1. **Composition refactor first** (§6): extract the LevelHierarchy state and a RefinementBuilder interface, with the current pipeline as backend #1. Mechanical (the `LgrHelpers` layer is ~70% of the move), 1–2 months, and it protects the production path while the replacement matures.

2. **Preprocessor-level builder (S2) as backend #2** (~6 months): build refined grids at the `processed_grid` level, reusing `findconnections()` for all non-matching stacks. Faults and pinch-outs are then handled *by construction* with the same code that handles them on level 0, `make_edge_conformal` composes on top (VEM-ready leaf), and the short-term boundary patches from items 2–3 above can be deleted rather than maintained. Run both backends on the same deck matrix in CI and require equivalence on the currently-supported subset before switching the default.

3. **Keep dynamics decoupled.** The only static decisions that prejudge dynamic refinement favorably are the composition seam and parameter-space overlap caching — do both; defer trees until dynamic refinement is actually scheduled.

## Test gate for every step

A deck matrix run in CI from day one: unfaulted Cartesian, faulted, pinched, NNC, aquifer; serial and parallel; with and without wells inside the LGR; ECLIPSE reference output where licenses allow. Each short-term item above is "done" when its row of the matrix passes — this is what keeps the incremental fixes from regressing each other, given how interconnected the identification logic is.

---

# Part III — Follow-up: scope reductions, touching/overlapping LGRs, deck conversion and wells

## III.1 Does anything change if NNC and MULTZ/TRAN edits are *not* needed inside LGRs?

The structure of the recommendation survives; the short-term ordering changes:

- Part II item 1 (MULTZ, TRANX/Y/Z on the LGR leaf) drops almost entirely. The top short-term priorities become pinched cells, then fault-adjacent boundaries, then diagnostics.
- One residue remains even under the reduced scope: the current checks are *global* — e.g. MULTZ anywhere in the deck plus an LGR anywhere in the grid throws (`Transmissibility_impl.hpp`), even when the MULTZ region never touches a patch. Relaxing the checks from "deck contains X and grid has LGRs" to "X intersects an LGR" is cheap and should be done regardless, otherwise the restriction blocks decks it does not need to block.
- "No NNCs in refined cells" is worth keeping as a *permanent* restriction (refine around them): explicit deck NNCs inside a refinement box are rare, and the restriction costs little. Note this concerns explicit NNC keywords — geometric fault connections are a different mechanism and are exactly what Part II item 3 addresses.
- The long-term recommendation (composition seam + preprocessor-level builder) is unchanged: it is motivated by faults and pinch-outs, not by NNC/MULTZ.

## III.2 Restrictions that simplify the code a lot at little cost in value

In rough order of (simplification gained)/(value lost):

1. **Disallow LGR patches that touch each other** (require ≥1 coarse cell of separation, or merge adjacent boxes into one). This deletes the LGR-to-LGR identification layer — `replaceLgr1CornerIdxByLgr2CornerIdx`/`...FaceIdx...`, `sharedFaceTag`, the last-appearance bookkeeping — which is among the most complex and fragile code in LgrHelpers. Value lost is small: equal subdivisions are *already* required on shared faces, so two touching boxes can usually be expressed as one larger box; and ECLIPSE itself does not allow adjacent independent LGRs (that is what AMALGAM exists for), so the restriction matches the input format's semantics anyway.
2. **Degenerate parents auto-deactivated** instead of collapse-aware refinement (Part II item 2a as the end state, not a stopgap): avoids coincident-corner merging and zero-area-face logic entirely; loses a small amount of pore volume unless compensated.
3. **Single refinement level (no nested LGRs)**: removes the shifted-level arithmetic, `corner_history_` chains and the per-parent-grid repeated leaf rebuilds. Nested CARFIN is rare in practice — but this feature was just built and tested, so this is a "if starting fresh" observation rather than a proposal to remove it.
4. **Conforming-only interfaces** (equal factors on shared faces — already enforced): keep forever; it is what makes the mortar machinery unnecessary.
5. **Not worth restricting**: faults near patches. Wells sit near faults in real models; excluding fault-adjacent LGRs would remove much of the practical value (it is also exactly the gap Part II prioritizes).

## III.3 Current status: LGRs that touch or overlap

- **Touching through a shared face**: supported if and only if the subdivisions match on the shared face — `compatibleSubdivisions` is checked (synchronized across ranks) and throws otherwise; covered by `lgrs_sharing_faces_test`.
- **Touching through an edge or corner only**: *no compatibility check exists* — `patchesShareFace` detects face sharing only. The corner-identification machinery handles shared edges between equal-factor refinements; for *different* factors meeting along an edge nothing validates the situation, and duplicate coincident leaf corners are the likely outcome. Untested corner case; should either be checked or covered by a test.
- **Overlapping boxes: not validated at all.** `validStartEndIJKs` checks only `start < end` and size consistency; the "disjoint patches" promise in the `CpGrid.hpp` documentation is not enforced anywhere. A cell inside two boxes gets assigned to whichever LGR comes last in the marking loop (`markElemAssignLevelDetectActiveLgrs` overwrites `assignRefinedLevel`), leaving the first LGR with holes, and the CARFIN local Cartesian indexing (`computeGlobalCellLgr` assumes a full box) then produces wrong local indices — i.e. **silent corruption, not an error**. ECLIPSE requires disjoint CARFINs, so a disjointness check that throws is the right fix and is cheap. This is a concrete addition to the §1 restriction list: it is the one unsupported case that fails silently rather than loudly.
- **Effect on the review**: this strengthens III.2 item 1 — a large share of the identification complexity exists precisely to serve touching refinements, while the touching cases that users actually need are largely expressible by merging boxes; and it adds one real bug-risk finding (unchecked overlap) to an implementation that otherwise fails safely.

## III.4 Putting an LGR into an existing deck — especially wells

**Grid section.** Add `CARFIN ... ENDFIN` (name, IJK box, nx ny nz) — that is the supported input surface: `LgrCollection` reads only CARFIN records (RADFIN is an explicit TODO), and the LGR control family (`LGR`, `LGRCOPY`, `LGRFREE`, `LGRLOCK`, `LGRON/OFF`) is flagged as unsupported-error in flow.

**Properties.** Refined cells inherit all static properties from their parent cell (`LookUpData`). Property keywords placed *inside* a CARFIN…ENDFIN block are parsed into `CarfinManager`, but I found no consumer of it — per-LGR property edits are silently **not applied**. If the deck relies on, e.g., PERMX edits inside the LGR block, the run will be wrong without warning. (Worth a parse-time warning until implemented.)

**Wells — the main conversion task.** Any well perforating cells inside a refinement box must be rewritten to the local-grid form: `WELSPECL` (declares the well in a named LGR) and `COMPDATL` (connections in LGR-local IJK). Both are implemented end-to-end (Schedule handlers → per-LGR `CompletedCells` → `BlackoilWellModel` `is_lgr_well`/per-connection grid ids → LGR-aware output). Defaulted connection factors are then recomputed from the *refined* cell geometry — which is the accuracy benefit LGR exists for.

**Pitfall (found in this review): leaving a global well inside the box is not caught.** `cartesianToCompressed_` (FlowBaseVanguard) is built by looping leaf cells, and all refined children share their parent's level-0 Cartesian index — so the map is last-writer-wins and a global `COMPDAT` inside a CARFIN box silently attaches the connection to one arbitrary refined child, with a connection factor computed from level-0 data. Deck conversion must therefore be complete (every affected well converted), and flow should gain a guard that throws when a global connection resolves to a refined cell.

**Placement guidance.**
- Keep the box at least one cell away from anything the implementation refuses: fault-split cells, pinched cells, explicit NNCs, aquifer cells — and from other boxes (III.3).
- Choose odd nx/ny factors so the well column stays centered in refined cells (`AUTOREF` enforces odd factors for exactly this reason; CARFIN does not, so it is on the deck author).
- Size the box so the well is interior, with a few refined cells of buffer before the LGR boundary — the boundary cells carry the coarse-fine transmissibility transition.
- A well cannot straddle the LGR lateral boundary: it must be entirely global or entirely within one LGR (standard ECLIPSE limitation; amalgamation, which lifts it, is not supported).
- Parallel runs: wells inside LGRs are supported with zoltanGoG partitioning (well cells kept together; covered by the distributed-LGR-with-wells tests).

---

# Part IV — After lgr_refactor: copy, new grid, or in-place refactor — and repository strategy

Background: the refactor attempt github.com/hnil/opm-grid/tree/lgr_refactor (the `CpGridLGR : public CpGrid` subclass, §6) was not received positively for merging. This part answers: does further development need a complete new CpGrid, or a copy of the old one to start from — and what would one change in CpGrid starting from scratch? Scope: the goal is **static LGR** (Part II); dynamic refinement is recorded but deferred.

## IV.1 New grid, copy, or in-place refactor — for static LGR

Three options, with a clear recommendation:

- **In-place refactor upstream: not viable as the leading move.** The lgr_refactor reception is direct evidence that invasive restructuring of the production grid will not pass review when presented as a refactor. Structural change has to arrive as a *replacement with evidence*, not as a rearrangement PR.

- **A complete new grid implementation ("Route B"): not needed for static LGR.** The decisive technical observation is that `CpGridData` can already *represent* everything static general LGR produces: `face_to_point_` is variable-length (`SparseTable`), `cell_to_point_` already holds collapsed hexahedra on level 0, and leaf cells at LGR boundaries already carry more than 6 faces in `cell_to_face_`. The blockers catalogued in §1–§2 live in the refinement *build algorithm* (index-arithmetic identification, 6-logical-face assumptions, the monolithic rebuild) — not in the container. A new grid core only pays for itself if VEM-on-general-grids and *dynamic* refinement become firm goals; its design notes are recorded in IV.2 and §8.4.

- **Recommended ("Route A"): a copy of CpGrid with the current LGR layer removed, plus a new builder component — hosted in `opm-gridrefined`.** Keep the mature and expensive-to-reproduce parts untouched (preprocessor including `make_edge_conformal`, parallel machinery, ECL I/O, the DUNE facade); strip the LGR/adaptivity layer (`LgrHelpers.*`, `refineAndUpdateGrid` and helpers, `NestedRefinementUtilities`, corner-history and level bookkeeping — `mark`/`adapt` stay as inert stubs so the DUNE interface remains intact); rebuild static refinement as a separate builder component behind a §6-style interface, preferably the §8.2-S2 preprocessor-level builder (faults and pinch-outs by construction, composes with `make_edge_conformal` for VEM).

**Fork discipline.** `opm-gridrefined` carries the full opm-grid history (seeded from `OPM/opm-grid`, `upstream` remote configured) — that is what makes rebasing on upstream and eventual upstream PRs possible. It is *not* a GitHub-fork (the fork slot is taken by `hnil/opm-grid`, which stages final PRs); a standalone repo with shared history is equivalent. Rebase regularly; the diff is mostly deletions plus the new layer, so rebases stay manageable. The copy is a development vehicle targeting replacement, not a permanent fork.

**Repository structure (keeping opm-simulators conflict-free).** Three layers:
1. **opm-gridrefined** — the grid module (module name stays `opm-grid`): grid work and grid-level tests/tools only. It cannot host simulator executables (opm-simulators depends on opm-grid, not vice versa).
2. **opm-simulators — untouched.** Because Route A preserves the `Dune::CpGrid` class name and public facade, unmodified opm-simulators rebuilds against this module; plain `flow` works (without LGR while the builder matures).
3. **A separate thin app repo (e.g. `opm-flowrefined`), created when needed**: the experimental flow executable — `main()` + TypeTag + custom Vanguard wiring the new builder API — plus the deck-matrix integration tests. The pattern is proven by opm-geomech (custom executables on top of installed opm-simulators). The custom Vanguard only becomes necessary once the builder API diverges from `addLgrsUpdateLeafView`; until then the app repo is deferred.

**Bring-up order** (acceptance gate = the Part II deck matrix): serial unfaulted decks → pinched/faulted decks → wells-in-LGR → parallel.

## IV.2 What to change when rebuilding — and from-scratch notes

Design rules for the rebuilt refinement layer (Route A):

1. **Levels as views, not copies.** One corner/geometry pool; level grids and the leaf are index views into it. This eliminates the leaf geometry duplication, `corner_history_`, and most of the identification maps in one stroke.
2. **Global refined-entity keys** — (parent cell, local ijk) hashes — instead of pairwise index translation (`replaceLgr1CornerIdxByLgr2CornerIdx` and friends); geometric matching at boundaries by 2D clipping in pillar parameter space (§8.1).
3. **A single distribution story.** Do not extend the `data_`/`distributed_data_`/`current_data_` dual-state into the new layer — it is the source of the add-LGRs-twice Vanguard dance and the global-id synchronization complexity. Refinement state is defined on the distributed grid only, with global ids stable by construction (root cell id, child index).
4. **First-class non-matching connections.** Fault-side LGR boundaries, NNCs and (interim) nonconforming couplings share one explicit connection representation with overlap geometry, feeding the same assembly.

Recorded for a future Route B (complete new core — only if VEM on general grids and dynamic refinement become firm goals): polyhedral-first topology (drop the `array<int,8>`/6-face assumptions entirely); no friend web — narrow public topology/geometry accessors, entities as index handles, core usable without DUNE (which the VEM/geomech side wants anyway); forest-of-trees refinement state with Morton ids and 2:1 balance (§8.4); geometry computed/cached rather than stored per level; distribution as a layer over the same data. The honest cost statement stands: CpGrid's parallel machinery and ECL integration are the moat — year-scale to parity — which is exactly why Route A keeps them.
