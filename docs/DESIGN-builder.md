# Static refinement builder — design

Design for Track 1 milestones 2 (remainder) and 3 of [PLAN.md](PLAN.md). Grounded in the post-strip code (`strip-lgr` branch): every claim about "what must be populated" was checked against the surviving readers, not against what the old implementation used to fill.

## 1. Scope

A backend implementing `Opm::Refinement::Builder` (the seam) that creates conforming block refinements (CARFIN semantics) on general corner-point level-0 grids. First increment: unfaulted, non-degenerate patches (milestone 3); faults/pinch-outs follow (milestone 4) without changing this design — they extend the matching stage only.

## 2. The contract: what a builder must populate

Inventory of multilevel state with *surviving* readers (everything else the old machinery maintained is dead weight and must not be resurrected):

| State | Surviving readers | Notes |
|---|---|---|
| `data_` levels + leaf, `level_data_ptr_`, `level_` | CpGrid facade, Entity, iterators | leaf is `data_.back()` |
| per-level topology (`cell_to_face_`, `face_to_cell_`, `face_to_point_`, `cell_to_point_`, `face_tag_`, `face_normals_`) | everything | leaf needs its own (DUNE facade) |
| per-level geometry (`geometry_`) | everything | **corners shareable, see D1** |
| `child_to_parent_cells_`, `cell_to_idxInParentCell_` | `Entity.hpp` (`father()`, `geometryInFather()`) | |
| `leaf_to_level_cells_` | `Entity.hpp`, `Indexsets.hpp` | |
| `cells_per_dim_`, `refinement_max_level_` | `Entity.hpp` (`isLeaf`, `geometryInFather`) | |
| `global_cell_`, `logical_cartesian_size_` | CartesianIndexMapper, LookUpData, output | CARFIN local indices for level grids |
| `corner_history_` | `Indexsets.hpp` only (corner ids) | **eliminated by D1+D3** |
| `lgr_names_` | CpGrid facade, output | |
| index sets / id sets per level + leaf | Indexsets, parallel, output | **replaced by D3** |
| `level_to_leaf_cells_` | none survive | drop — do not populate |

## 3. Design decisions

### D1 — One shared corner pool across all levels and the leaf

`DefaultGeometryPolicy` stores its geometry vectors as `shared_ptr<EntityVariable<...>>`, so all level grids and the leaf can point at a *single* corner vector: level-0 corners first, then refined corners appended per level. Cells and faces reference corners by index into the pool.

Consequences, all simplifying:
- "Is this refined corner the same as that one" disappears — same pool index *is* the same corner. The entire corner-identification layer of the old implementation (its largest and most fragile part) has no counterpart here.
- `corner_history_` becomes unnecessary: a leaf corner's identity is its pool index. Its only surviving reader (`Indexsets.hpp` corner ids) is rewritten to use the pool index (D3).
- Corner geometry is stored once — the old leaf duplicated every corner.
- Shared-face corners between adjacent boxes coincide by construction *if* the builder generates boundary corners deterministically from (parent face, lattice position) — which the conforming case guarantees.

Faces and cells keep per-level/leaf entries (their geometry is small: centroid + area/volume); only codim-3 is pooled.

### D2 — Leaf: materialized topology, shared geometry

A fully "leaf as view" design is not reachable without an index-indirection layer in Entity/Iterators (they assume each `CpGridData` owns dense tables). Decision: keep the leaf's *topology* materialized (int tables — cheap), share the *corner geometry* through D1, and accept duplicated face/cell centroid records. This preserves the DUNE facade unchanged and still removes the expensive duplication.

### D3 — Construction-stable ids (no insertIdSet, no sync)

`IdType` is `std::int64_t` (confirmed in `Indexsets.hpp`). Scheme:
- Level-0 entity → its existing level-0 id (unchanged, so unrefined runs are bit-identical).
- Refined cell → `pack(parent level-0 cell id, child index within parent)`; child index is the (i,j,k) lattice index in the parent's subdivision — deterministic, partition-independent.
- Refined corner/face → `pack(owning parent entity id, lattice position)`, owner chosen by a fixed rule (lowest parent cell id touching the entity) so both sides of a shared boundary compute the same id without communication.
- Packing: high bits parent id, low bits child code; bounds checked at build time (a 2^31 parent id with 2^20 children fits comfortably in 63 bits). If bounds ever fail, fall back to a deterministic hash map built identically on every rank.

This removes `insertIdSet`, id prediction/winner-selection, and `syncDistributedGlobalCellIds` for the new path: ids agree across ranks *by construction* because they are functions of level-0 ids (already synchronized) and deterministic local codes.

### D4 — Build pipeline (S2, staged)

Small stages exchanging value types (LESSONS #6), each unit-testable:

1. **Plan**: requests → per-column subdivision map; reject unsupported content for now (NNC/aquifer/degenerate cells in or adjacent to boxes) loudly.
2. **Sub-pillars**: per box, lateral interpolation of new pillars on the bilinear surfaces between the four surrounding original pillars (COORD lines are straight; interpolation is linear in the pillar endpoints).
3. **Refined ZCORN**: per parent cell, trilinear interpolation of the 8 zcorn values onto the sub-pillar lattice. Semantics note: this is the *corner-point-native* refinement (refined corners lie on straight sub-pillars — what ECLIPSE LGRs mean). It coincides with the old per-cell trilinear-map refinement **only for vertical pillars**; on inclined pillars the two differ slightly (both are valid; we choose ECLIPSE semantics). Verified properties (grdecl_refinement_test): exact uniform-lattice reproduction; volume partition on vertical-pillar grids (distorted surfaces, faults inside the block); zcorn jumps and ACTNUM inherited. Implemented: `refinement/GrdeclRefinement.{hpp,cpp}` (stages 2+3, fault-ready in the block interior). Input note: `CpGridData::zcorn` is only retained when MINPV created NNCs, so the builder takes COORD/ZCORN explicitly; the integration layer decides whether to retain them at construction (when the deck has LGRs) or pass the EclipseGrid.
4. **Window processing**: run the corner-point preprocessor on the box plus a one-cell halo, producing a `processed_grid` for the refined region. Faults inside or at the box boundary are handled *here* by `findconnections` in milestone 4 — no change to other stages.
5. **Stitch**: identify window-boundary entities with level-0 entities. Conforming case: boundary corners coincide with parent-face lattice points → direct index identification into the D1 pool. Faulted case (milestone 4): parameter-space clipping (review §8.1) produces the leaf face mosaics.
6. **Assemble**: populate level grids and leaf per §2's table, ids per D3.

Property-based acceptance per stage and end-to-end (LESSONS #16): child volumes sum to parent volumes; refined face areas partition parent faces; leaf face_to_cell symmetric and closed; ids stable under rank count; unrefined decks bit-identical to pre-builder runs.

### D5 — Order of capability

1. Serial, single box, unfaulted interior and boundary (milestone 3).
2. Multiple boxes incl. face-touching (shared corners via D1 determinism).
3. Pinched/degenerate parents (collapse handled by the preprocessor in stage 4; stitching unaffected).
4. Fault-adjacent boundaries (stage 5 clipping).
5. Distributed grids: refinement after loadbalance on the distributed view only — D3 makes ids rank-independent, so no global-view duplication (`addLgrs` twice) and no sync step.

## 4. Open questions

- `logical_cartesian_size_` of the leaf when boxes exist: upstream sets it to the level-0 size; output paths assume this. Keep that convention until output is revisited.
- Nested refinement: representable (parent grid name in `BlockRefinement`), deferred until the flat case passes the deck matrix.
- `adapt()`/mark-based general refinement: out of scope for the static builder; the seam can grow a second entry point later.

## 5. What is deliberately *not* built

No LevelHierarchy class for its own sake: §2's table *is* the state, living where it always lived (`CpGridData`), so the simulator stack stays untouched. The builder is the only writer; if a second writer ever appears, extracting the table into a component is mechanical.
