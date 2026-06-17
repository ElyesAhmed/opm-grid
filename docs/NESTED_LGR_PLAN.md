# Nested LGR (CARFIN inside an LGR) — plan & status

A *nested* LGR is a `CARFIN` whose parent is another LGR, not `GLOBAL`
(gap A1 in `LGR_GAPS.md`). This document is the implementation plan plus the
current landing status. **Overriding constraint: the single-level
(`parentGridName == "GLOBAL"`) path must stay byte-identical.** Every
nested-specific branch is guarded so the existing decks are unaffected; this
is verified by the regression (CARFIN1 serial + `conforming_builder_test`)
after every step.

## Cross-layer state (investigation, 2026-06-16)

| Layer | Nesting support | Where |
|-------|-----------------|-------|
| Deck parsing (opm-common) | **already works** | `Carfin::PARENT_NAME()` (default `"GLOBAL"`), `LgrCollection`, `readKeywordCarfin` |
| EclipseState grid + ECL output (opm-common) | **already works** | `EclipseGridLGR` recursive father chain `get_*_child_to_top_father` (walks `father_label != "GLOBAL"`); `perform_refinement(parent_coord, parent_zcorn, parent_nxyz)` is parameterized on the parent geometry |
| Simulator vanguard (opm-simulators) | **gap — dropped** | `GenericCpGridVanguard::addLgrsUpdateLeafView` calls the 4-arg `grid.addLgrsUpdateLeafView(...)`, never passing parent names; no parent-before-child ordering |
| Grid builder (opm-gridrefined) | **gap — throws** | `ConformingBlockBuilder::build` throws on `parentGridName != "GLOBAL"`; `assembleBlockLevelGrid` pins parent grid to level 0 (`childToParent = {0, …}`); `assembleLeafGrid` assumes a flat 2-level hierarchy |

Already-nesting-aware infrastructure: `CpGrid::getLgrNameToLevel()` /
`lgr_names_` (name→level); `validateBlockRefinements` only enforces
disjointness *within the same parent*; the leaf assembler already routes
sources through `(src.grid == 0 ? level0 : boxes[src.grid-1])`.

## Phases

### Phase A — simulator: pass the parent through  *(small)*
`GenericCpGridVanguard::addLgrsUpdateLeafView`:
- build `lgr_parent_grid_name_vec` from `lgrCarfin.PARENT_NAME()` and call the
  5-arg `grid.addLgrsUpdateLeafView(...)`;
- **topologically order** boxes so a parent is always added before its child
  (stable sort by depth in the parent chain), since the builder appends level
  grids in request order and a child must resolve an already-built parent;
- a nested `CARFIN`'s `I1..I2/J1..J2/K1..K2` are already expressed in the host
  LGR's local Cartesian space, so they map directly to `startIJK/endIJK`;
  `cellsPerDim = N{X,Y,Z}/(extent)` is likewise parent-local.

GLOBAL-only decks: parent vector is all `"GLOBAL"`, the sort is a no-op
(depth 0), indices unchanged → identical.

### Phase B — builder: build a level grid from an arbitrary parent level  *(medium)*
`ConformingBlockBuilder::build` + `assembleBlockLevelGrid`:
- resolve `parentGridName → parentLevel` (0 for `"GLOBAL"`);
- keep, per built level, its `RefinedBlockGrdecl` (dims/coord/zcorn/actnum);
  for a nested box, refine from the **parent level's** refined geometry
  instead of the global `coord_/zcorn_/dims_`;
- record `childToParent = {parentLevel, parentIdx}` (was hardcoded `0`) and
  invert the *parent level's* `globalCell()` rather than level 0's;
- drop the blanket early throw; the unsupported boundary now lives in the
  leaf assembler (Phase C) with a precise message.

GLOBAL boxes (parentLevel 0): same geometry, same `{0, parentIdx}` → identical.

### Phase C — leaf: recursive stitching  *(large — the remaining work)*
`assembleLeafGrid` currently does one coarse→fine pass: iterate level-0 cells,
replace a refined parent by its children, identify box corners/faces with
level 0. Nesting needs this to become a walk over the **refinement tree**:
- a leaf cell is the finest descendant; a nested box's cells replace its parent
  *LGR's* cells (themselves refined), not level-0 cells;
- generalize `boxOfCell`, `childrenOfParent`, corner identification
  (`cornerEquiv` via `cellToPoint0[parent]`) and the faulted/mosaic boundary
  logic from "parent = level 0" to "parent = any level grid";
- keep composite global ids (`levelOffset[level] + level-Cartesian idx`)
  consistent across >2 levels.

Until Phase C lands, a nested deck builds its level grids (Phase A/B) and then
throws a precise *"nested leaf assembly not yet implemented"* at the leaf — no
wrong output, GLOBAL path untouched.

### Phase D — parallel nested (rank-interior)  *(large — defer)*
`classifyBox` must classify a nested box against its parent LGR's owner rank
and keep the whole nested subtree on one rank. Serial first.

### Phase E — tests
- grid unit test: a nested `CARFIN` builds, volume conserved, father chain
  (child→parent-LGR→GLOBAL) correct, ids consistent, edge-conformal across the
  nested boundary;
- flow deck: `opm-tests/lgr/SPE1CASE1_CARFIN1_NESTED.DATA` (a CARFIN inside
  CARFIN1's box) runs on the fork once Phase C lands.

## Index convention (confirmed 2026-06-16)

A nested `CARFIN`'s `I1..I2/J1..J2/K1..K2` are **parent-LGR-local** (1-based in
the parent LGR's own refined Cartesian space), *not* global. opm-common resolves
them via `getActiveIndex` on the parent LGR grid object
(`EclipseGrid.cpp:2561/2567`), and the simulator/builder treat them as
parent-local too — so Phase A/B are consistent end to end, no index translation
is needed.

Caveat (separate opm-common gap): `Carfin.cpp:assert_dims` validates *every*
CARFIN's indices against the **global** grid dims (the `Carfin` is always built
with `m_globalGridDims_`). That wrongly rejects legitimate parent-local indices
that exceed the global dims (e.g. a nested box with `K2 > global NZ` even though
`K2 <= parent NZ`). Nested decks must currently keep indices within global dims;
fixing the validation to use the parent LGR's dims is an opm-common task.

## Status
- [x] Investigation / plan (this file)
- [x] Phase A — simulator parent pass-through + topological ordering
      (`GenericCpGridVanguard::addLgrsUpdateLeafView`); GLOBAL-only decks
      byte-identical (CARFIN1 serial UNRST+INIT verified)
- [x] Phase B — level-grid build from an arbitrary parent level
      (`assembleBlockLevelGrid` takes `parentGrid`/`parentLevel` + exposes its
      `RefinedBlockGrdecl`; `ConformingBlockBuilder::build` threads each level's
      geometry and resolves `parentGridName -> level`). GLOBAL path identical;
      nested builds its level grids then stops at the Phase-C leaf boundary.
- [ ] Phase C — recursive leaf assembly (the remaining large piece; builder
      currently throws a precise staged message after building the nested
      level grids)
- [ ] Phase D — parallel nested (guarded: throws in parallel for now)
- [x] Phase E (partial) — grid unit tests `nestedRefinementReachesLeafBoundary`
      and `nestedChildBeforeParentThrows`; flow deck
      `opm-tests/lgr/SPE1CASE1_CARFIN1_NESTED.DATA` (builds level grids, then
      the documented Phase-C throw). Full leaf/volume asserts pending Phase C.
