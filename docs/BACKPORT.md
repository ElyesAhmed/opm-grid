# Learnings from the rebuild that can extend the *original* upstream LGR code

The rebuild's purpose is to handle what the original (`Opm::Lgr` +
`refineAndUpdateGrid` in upstream opm-grid) structurally cannot. But several
findings made along the way are discrete improvements to the **original**
code, independent of the rebuild — small, mergeable upstream PRs. This
catalogues them by value/effort, and is explicit about what cannot be
back-ported. Claims below were checked against the unmodified opm-grid /
opm-simulators checkouts.

## A. Directly back-portable bug fixes (small, high value, low risk)

These are defects in the original, fixable in isolation.

1. **Overlapping CARFIN boxes are not validated — silent corruption.**
   Upstream `validStartEndIJKs` (LgrHelpers.cpp) checks only `start < end`
   per box; it never checks that boxes are pairwise disjoint. Overlapping
   CARFIN boxes assign a cell to whichever LGR is marked last, leaving the
   other LGR with holes and wrong local Cartesian indices — no error
   (review Part III.3). Fix: a pairwise-disjointness check per parent grid
   (the rebuild's `validateBlockRefinements`). ~20 lines, pure validation.

2. **A global COMPDAT inside a CARFIN box is silently mis-mapped.**
   `cartesianToCompressed_` (FlowBaseVanguard) is built last-writer-wins
   over refined children sharing the parent's Cartesian index, so a global
   well connection landing in a refined region attaches to an arbitrary
   refined child with a level-0 connection factor (review Part III.4). Fix:
   when `grid.maxLevel()>0`, throw if a global connection resolves to a
   refined cell, telling the user to convert to WELSPECL/COMPDATL. (This is
   the task already flagged in this session.)

3. **MULTZ / TRANX-Y-Z are rejected globally, not where they matter.**
   `Transmissibility_impl.hpp` throws on `grid_.maxLevel() > 0` — i.e. *any*
   LGR anywhere disables MULTZ/TRAN edits *everywhere*, even when the edit
   region never touches an LGR (review Part III.1). Fix: scope the check to
   "the edit intersects a refined region". Unblocks decks that need not be
   blocked.

## B. Generic CpGrid parallel/partitioning improvements

These are CpGrid-level (not tied to the refinement implementation), so they
help the original's *distributed* LGR path too.

4. **`setPartitionCellGroups` — keep cell groups (LGR boxes) on one rank.**
   A small API on CpGrid plus a GraphOfGrid contraction
   (`addPartitionCellGroups`, reusing the well-cell `addWell` mechanism).
   The original distributes level zero first and then handles LGRs that may
   be split across ranks with substantial communication; keeping each box on
   one rank removes that need for the common case. Directly portable — it
   does not depend on the rebuild.

5. **Overlap layer 2 fixes the contracted-region scatter.**
   The original itself documents the overlap-layer-1 limitation (CpGrid.cpp:
   "due to the overlap layer size (equal to 1) cells that share corners or
   edges (not faces) with interior cells are not seen by the process"). The
   rebuild confirmed that contracting an interior region under overlap 1
   trips `computeFace2Cell`'s assertion, and that **overlap 2 resolves it**
   (2-rank). Insight, not a drop-in: it suggests the original's parallel LGR
   could be simplified by combining box-on-one-rank partitioning with
   overlap≥2, replacing much of the id-prediction / winner-selection /
   overlap-communication machinery. (≥4-rank contracted scatter still needs
   a deeper CpGrid fix — relevant to both implementations.)

6. **Rank-local grids must use a self-communicator.**
   Building a rank-local grid with the world communicator makes any
   collective inside it (e.g. `processEclipseFormat`'s broadcast) deadlock
   when only one rank participates. The original sidesteps this (its level
   grids are built by index arithmetic, no `processEclipseFormat`), but the
   lesson is general for anyone extending the parallel paths.

## C. Verification tooling (implementation-agnostic — works on the original)

7. **Connection-keyed output comparison (`scripts/compare_lgr_output.py`).**
   Full element-wise INIT comparison per LGR section, plus boundary
   transmissibilities (TRANGL) keyed by the EGRID `(NNCG, NNCL)` cell pairs
   so it is immune to face-numbering differences. Nothing in it is specific
   to the rebuild; it is a ready regression/CI tool for the original
   (before/after refactors, or OPM-vs-reference).

8. **Structural-equivalence oracle testing.**
   Comparing a refined leaf against a grid built *directly* from the refined
   corner-point description (the upstream LgrChecks spirit, sharpened). It
   applies to the original on its supported (unfaulted) cases and is a
   stronger check than the property tests.

## D. What cannot be back-ported (this is why the rebuild exists)

- **Faults inside refinement blocks.** The original throws ("more than six
  faces … not supported"). Its entity identification is closed-form index
  arithmetic that assumes each parent face is one of six logical faces with
  a matching refined mosaic — fundamentally incompatible with fault-split
  faces. The rebuild handles faults only because it refines at the
  *preprocessor* level (`findconnections` matches the non-matching stacks).
  There is no small patch that gives the original this; it is the central
  motivation for the rebuild.
- **Pinch-outs / degenerate (<8-corner) parents** — same root cause.

## E. The one architectural lever worth stating plainly

If the original were ever to gain fault/pinch-out support without the full
rebuild, the only viable route is the same one the rebuild took: express
refinement at the `processed_grid` level and let the existing preprocessor
(`findconnections`, `make_edge_conformal`) do the geometric matching. That
is effectively the rebuild, which is why "extend the original for faults" and
"adopt the rebuild" converge.

## Suggested upstreaming order (independent of the rebuild)

1. A.1 disjointness check and A.2 well-mapping guard — pure safety, no
   behaviour change for valid decks.
2. C.7 comparison tool — immediately useful for any LGR CI.
3. A.3 MULTZ/TRAN scoping — unblocks real decks.
4. B.4 `setPartitionCellGroups` — generic partitioning improvement.

Items in B.5/B.6 and D are insights/limitations to record rather than PRs.
