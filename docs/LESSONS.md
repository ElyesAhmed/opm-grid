# Lessons learned from the current CpGrid LGR implementation

Distilled from the 2026-06 review ([lgr_review.md](lgr_review.md)) for whoever does the next refactoring — what to keep, what to avoid, and the pitfalls that were found the hard way. Section references point into the review.

## Architecture lessons

1. **The container was never the problem — the build algorithm is.** `CpGridData` can already represent everything static general LGR produces: variable-length `face_to_point_` (SparseTable), collapsed hexahedra in `cell_to_point_` on level 0, leaf cells with >6 faces at LGR boundaries. A refactoring that starts by redesigning storage is solving the wrong problem (§IV.1).

2. **Index arithmetic does not scale to geometry.** The closed-form index translation between neighboring refinements (`replaceLgr1CornerIdxByLgr2CornerIdx` and friends) hard-codes tensor-product/6-logical-face assumptions in dozens of places, and it is the single reason faults, pinch-outs and partial face overlap are unsupported. Any boundary generality requires *geometric* matching — and the key simplification is that all required matching is **2D in pillar-pair parameter space** (rectangles vs straight-edged quads; exact with snapped integer coordinates) (§2, §8.1).

3. **Subclassing the grid fails — proven, not just argued.** The `CpGridLGR : public CpGrid` attempt moved ~1400 lines of declarations but no coupling (`CpGridData` kept all multilevel state, the friend web stayed, a second grid type would propagate through every simulator template) and was rejected upstream. Separation must be by **composition**: a state component (LevelHierarchy) plus a builder interface, behind an unchanged `CpGrid` facade (§6).

4. **A materialized leaf is expensive in code, not just memory.** Duplicating geometry into the leaf forces the `corner_history_` machinery and most of the ~20 `std::map<std::array<int,2>,…>` relation maps that make the pipeline hard to follow. Levels-as-views into one corner/geometry pool eliminates whole categories of bookkeeping at once (§IV.2-1).

5. **Dual serial/distributed state is a complexity multiplier.** The `data_`/`distributed_data_`/`current_data_` pointer-swap plus `switchToGlobalView()/switchToDistributedView()` is what forces adding LGRs *twice* in the Vanguard, `syncDistributedGlobalCellIds()`, and the global-id prediction/winner-selection machinery. Define refinement on one (the distributed) view only, with ids stable by construction — (root cell id, child index) needs no synchronization (§IV.2-3, §7).

6. **Big-bang orchestration resists testing and change.** One ~540-line function feeding 15+ helpers with 10–20 parameters each cannot be tested in pieces; every fix risks the whole. Per-marked-cell throwaway `CpGridData` objects add allocation cost for no architectural gain. The replacement should be a pipeline of small stages exchanging plain value-type artifacts, each testable alone (§1).

## Correctness pitfalls found (do not rediscover these)

7. **The two silent failures were both missing input validation, not algorithm bugs.** (a) Overlapping CARFIN boxes are never checked — overlapped cells go to whichever LGR is marked last, corrupting the local Cartesian indexing without any error (§III.3). (b) A global COMPDAT inside a CARFIN box silently attaches to an arbitrary refined child, because `cartesianToCompressed_` is built last-writer-wins over leaf cells sharing the parent's Cartesian index (§III.4). Lessons: validate deck-level invariants eagerly (disjointness), and never build index maps whose keys can silently collide — with multilevel grids, Cartesian→compressed must be level-aware or refuse ambiguity.

8. **Scope restrictions to intersections, not existence.** MULTZ anywhere in the deck plus an LGR anywhere in the grid throws, even when they never touch (§III.1). Global feature-exclusion checks block decks they do not need to block; check "does X intersect a refined region" instead.

9. **Edge conformity must be (re-)established *after* refinement.** The refined leaf is face-conformal but not edge-conformal (hanging nodes on edges of unrefined faces), and LGR on a grid *built* with `edge_conformal=true` would drop the inserted nodes during face replacement. The `make_edge_conformal` post-pass pattern is correct — but it has to run on the leaf, not only at preprocessing (§3).

10. **Deck-format semantics leak if not layered.** `logical_cartesian_size_` is undefined for non-block refinements (open TODO in the code), `adapt()` hard-codes `{2,2,2}`, and naming/output/well plumbing assume named block LGRs. Keep CARFIN concepts (names, local Cartesian indices, ECL output) in a layer above the refinement kernel (§5).

## Design rules that paid off (keep these)

11. **Conformity-by-invariant beats coupling-by-computation.** Requiring equal subdivisions on shared LGR faces made mortar machinery unnecessary. The generalization for any future dynamic work is a 2:1-balance-style grading rule: design the invariants that bound interface cases first, then code to them (§8.4).

12. **Reference-space refinement with the trilinear map is the right kernel — with one precise boundary qualifier.** The trilinear map of a hex's 8 corners refines that *single* cell's interior geometry correctly for any hexahedron, including degenerate ones (corners of collapsed edges map to coincident points — only topology dedup is missing). Across a **matching, non-faulted** shared face the refined sub-faces of the two neighbours are automatically consistent because both sides restrict to the *same* bilinear surface — so for hexes **not touching a fault** the trilinear map alone gives a conformal mosaic, no mortar/clipping needed. This consistency does **not** hold where a cell touches a fault: the two sides of a fault face are *different* bilinear surfaces (the fault offsets the corners), so the interface is not a simple sub-face mosaic and must be matched by **intersection** — the preprocessor's `findconnections` / §8.1 parameter-space clipping — rather than assumed from the trilinear map. (The current builder handles faults *inside* a refined block precisely because it runs the preprocessor per block; the trilinear-only mosaic assumption is what holds at the non-faulted box boundaries.) Fault overlap polygons in parameter space are *refinement-invariant*, hence computable once and cached (§8.1, §8.4).

13. **The preprocessor is the crown jewel — reuse, don't reimplement.** `findconnections()` already solves non-matching stacks (faults, pinch-outs) per pillar pair; `make_edge_conformal` shows the hanging-node insertion pattern; `geometry.c` computes the derived quantities. Refinement expressed at the `processed_grid` level inherits all of it (§8.2-S2).

14. **Failing loudly was the right call.** Nearly all unsupported combinations throw instead of producing wrong results — which is why the review could map the support boundary precisely. Preserve this property in the rebuild; the two counterexamples (item 7) show the cost of the alternative.

15. **The simulator-side integration model works.** Per-connection grid ids and `is_lgr_well` through Schedule → well model → output is sound; `LookUpData` makes property re-evaluation on a new leaf straightforward; the Vanguard mechanism lets a new grid or builder come up without touching opm-simulators (opm-geomech pattern). Keep all of it (§1, §IV.1).

## Testing lessons

16. **Tests that mirror the implementation validate the implementation, not the behavior.** Much of the 33-file test suite re-derives the same index arithmetic it checks. Prefer property-based checks that survive a rebuild: child volumes sum to parent volume, refined face areas partition the parent face, global ids stable under repartitioning, leaf/level consistency — tied to the Part II deck matrix (unfaulted/faulted/pinched/NNC/aquifer × serial/parallel × wells in/out) as the acceptance gate. The `LgrChecks.hpp` harness is a good starting point and worth porting.
