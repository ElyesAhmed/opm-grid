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

1. **Strip**: remove `LgrHelpers.*`, `refineAndUpdateGrid` + helpers, `NestedRefinementUtilities`, corner-history/level bookkeeping; keep DUNE interface stubs; unmodified opm-simulators still builds and runs the non-LGR test decks. (Branch: `strip-lgr`.)
2. **Seam**: LevelHierarchy state component + RefinementBuilder interface (§6); ids by construction (root cell, child index) — no insertIdSet/global-id sync machinery (§IV.2).
3. **Builder, conforming core**: S2 preprocessor-level builder for unfaulted patches; levels as views into one corner/geometry pool (§IV.2-1).
4. **General corner-point**: pinched parents; fault-adjacent boundaries via `findconnections` + lateral splitting at sub-pillar positions (§8.2-S2); edge-conformal post-pass option for VEM (§3).
5. **Wells-in-LGR** (WELSPECL/COMPDATL through the existing opm-common path) + ECL output.
6. **Parallel**: refinement defined on the distributed grid only (single distribution story, §IV.2-3).

## Acceptance gate — deck matrix (Part II)

CI matrix from day one: unfaulted Cartesian / faulted / pinched / NNC / aquifer × serial / parallel × wells in/out of LGR; ECLIPSE reference output where available. Each milestone and Track 0 item is "done" when its rows pass.

## Out of scope (recorded, not planned)

Dynamic refinement (§7, §8.4: forest-of-trees state, 2:1 balance, Morton ids, cached fault overlaps) and a polyhedral-first new grid core (Route B, §IV.2) — revisit only if VEM-on-general-grids/dynamic refinement become firm goals.
