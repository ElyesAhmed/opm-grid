# Nested LGR — how to test

This document describes how to test nested LGR support (a `CARFIN` whose parent
is another LGR, `parentGridName != "GLOBAL"`). See `NESTED_LGR_PLAN.md` for the
design and phase status. Work happens on branch `nested-lgr-serial` (off
`refine-before-redist`).

## 0. Use a fresh build tree (important)

Switching a build directory between branches repeatedly leaves stale objects and
produced misleading parallel results once already. For nested work use a
dedicated build dir against the `nested-lgr-serial` checkout, and do **not**
branch-switch it:

```sh
# rbr_src/{opm-gridrefined,opm-simulators} on nested-lgr-serial, opm-common on refine-before-redist
cmake -S rbr_src/superbuild-refined -B builds/refined_nested -G Ninja \
  -DCMAKE_PREFIX_PATH=<dune-install> -DUSE_MPI=ON -DWITH_NATIVE=ON -DBUILD_TESTING=ON
ninja -C builds/refined_nested flow_blackoil conforming_builder_test
```

If you must reuse a churned build dir, force a clean recompile first
(`find builds/<dir>/opm-grid builds/<dir>/opm-simulators -name '*.o' -delete`)
before trusting any fine-grained (esp. parallel) comparison.

## 1. Grid unit tests (`conforming_builder_test`)

```sh
ninja -C builds/refined_nested conforming_builder_test
builds/refined_nested/opm-grid/bin/conforming_builder_test          # all cases
builds/refined_nested/opm-grid/bin/conforming_builder_test --run_test=nestedFullyContainedBuilds
```

Nested cases (in `tests/cpgrid/conforming_builder_test.cpp`):

| test | what it checks |
|------|----------------|
| `nestedChildBeforeParentThrows` | a child ordered before its parent is rejected |
| `nestedRefinementReachesLeafBoundary` | a nested box reaching the parent LGR's boundary builds its level grids (Phases A/B) then throws the precise Phase-C "not implemented" message |
| `nestedFullyContainedBuilds` | a child **strictly interior** to its parent LGR (the simplest case); same Phase-A/B-then-throw until Phase C lands |

**Touching (sibling) LGR tests** (same parent, adjacent) already pass and must
stay green: `twoSeparatedBoxesAndGuards`, `edgeSharingBoxesMergeCorners`,
`faceSharingBoxes`, `faceSharingNonMatchingSubdivisionsThrow`.

## 2. Phase-C acceptance (flip the contained test)

When the leaf assembler handles nesting, change `nestedFullyContainedBuilds`
from `BOOST_CHECK_EXCEPTION(...)` to the build call plus the asserts listed in
that test's header comment:
- `grid.maxLevel() == 2`, `getLgrNameToLevel() == {LGR1:1, NEST1:2}`;
- total leaf volume == volume before refinement (conservation);
- every `NEST1` leaf cell's father is an `LGR1` cell whose father is `GLOBAL`
  (child → parent-LGR → GLOBAL chain), `geometryInFather().volume() == 1/8`;
- leaf intersections are conformal (each interior face seen from both sides with
  matching area) — the same check `endToEndSingleBox` runs.

## 3. Flow deck

`opm-tests/lgr/SPE1CASE1_CARFIN1_NESTED.DATA` is a `CARFIN` inside `CARFIN1`'s
box. Serial only for now (parallel nested is Phase D).

**Contained variant (works, with ECL output):**
`opm-tests/lgr/SPE1CASE1_CARFIN1_NESTED_CONTAINED.DATA` — NEST1 is strictly
interior to LGR1 (and reaches LGR1-local k=8, beyond the global NZ=3, exercising
the parent-relative CARFIN validation, opm-common `LgrCollection::addLgr`). It
parses, builds the 3-level leaf, runs the SOLVE, and writes EGRID/INIT/UNRST:

```sh
builds/refined_nested/opm-simulators/bin/flow_blackoil \
  opm-tests/lgr/SPE1CASE1_CARFIN1_NESTED_CONTAINED.DATA \
  --parsing-strictness=low --output-dir=/tmp/nested
```
Expect `End of simulation` and `.EGRID/.INIT/.UNRST` (inspect in ResInsight).

The per-LGR ECL output requires the simulator's LGR level order to match
opm-common's deck order; Phase A (`GenericCpGridVanguard`) uses a stable
topological order (deck order, child after parent) for this. A depth-grouping
sort would mispair the per-LGR data (e.g. NEST1's trans with LGR2's grid ->
`EclipseGrid::compressedVector` "Input vector must have full size").

**Bundled `SPE1CASE1_CARFIN1_NESTED.DATA`** has a NEST that *touches* LGR1's
boundary in k, so it is (correctly) rejected with the containment message;
boundary-touching nesting is a future extension.

## 4. Regression — current behavior must not change

Nested support must keep every GLOBAL-parent (non-nested) deck byte-identical.
After any nested change, confirm:

```sh
# single-level LGR, serial, unchanged
builds/refined_nested/opm-simulators/bin/flow_blackoil \
  opm-tests/lgr/SPE1CASE1_CARFIN1.DATA --parsing-strictness=low --output-dir=/tmp/c1_nested
# vs a reference build of refine-before-redist (or master): compareECL must match
builds/refined_nested/opm-common/bin/compareECL -t UNRST \
  /tmp/c1_ref/SPE1CASE1_CARFIN1 /tmp/c1_nested/SPE1CASE1_CARFIN1 1e-9 1e-12
```

and the full `conforming_builder_test` stays green.
