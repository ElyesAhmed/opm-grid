# Design: a fully parallel, locally (dynamically) refinable CpGrid

This updates review §8.4 (forest-of-octrees) with what the static + parallel
rebuild on this branch actually proved. The earlier §8.4 was a sketch; this
is a concrete design whose building blocks are now mostly implemented and
tested. Read alongside [DESIGN-builder.md](DESIGN-builder.md) (state model,
shared corner pool, ids) and [PLAN.md](PLAN.md) (what is built).

## 0. What this session changed about the §8.4 plan

The original §8.4 assumed several things had to be invented. Most now exist
on `strip-lgr` and transfer directly:

| §8.4 assumption | Reality after this session |
|---|---|
| "a leaf-view builder that assembles the DUNE view by forest traversal" | **Built**: `assembleLeafGrid` merges level grids into a leaf — faults (preprocessor), edge/face-sharing boxes (merged corner pool), distributed level zero (off-rank sentinel, inherited partition types). Generalises from blocks to per-cell trees. |
| octant geometry via `refineSingleCell` trilinear map | Confirmed correct; for a single clean cell the per-octant geometry is just the trilinear image of its 8 corners (no sub-pillars needed *within* one cell — simpler than the block case). |
| "stable ids for free (root cell, Morton code)" | Validated as the right scheme — and now strongly motivated: parallel correctness hinges on ids, and rank-independent packed ids remove the id-prediction/winner-selection machinery entirely. The static builder used delegation; the octree should use packed ids. |
| "2:1 balance worth adopting" | Endorsed and now load-bearing: it makes edge-conformalization local and bounds the interface cases (below). |
| "faults stay cheap (refinement-invariant overlap polygons)" | Consistent with the built preprocessor-per-block approach; faults *inside* a refined region are matched by `findconnections`. The open part is faults at a refining cell's boundary (see §3). |
| parallel "exactly what p4est provides" | Partly built without p4est: rank-interior refinement works (np=2/4), partition types inherited, self-comm rank-local grids. But CpGrid's **cell-overlap scatter has real limits** (overlap-1 corner/edge gap; ≥4-rank contracted-region fault) — the octree must manage its **own ghost layer** rather than lean on that scatter (§5). |
| edge-conformal for VEM | **Built** as a generic leaf post-pass (`edgeConformalizeLeaf`); reused unchanged by the octree leaf. 2:1 balance makes it cheap. |

## 1. Data model

- **Macro-grid = level zero**, unchanged: the corner-point grid from the
  preprocessor (faults, pinch-outs, NNCs all already handled). It is the
  *forest connectivity*: each active level-0 cell is a tree root; root-to-root
  adjacency (incl. fault faces) is level zero's `face_to_cell`.
- **Refinement state = a forest of k-ary trees**, one per root cell. A node
  stores its per-direction split factors (octree = 2×2×2; anisotropic =
  rx×ry×rz). Leaves of the forest are the simulation cells. This *replaces*
  the accumulating `CpGridData` levels with compact per-cell data.
- **Leaf view** = a materialized `CpGridData` (the DUNE-facing grid) rebuilt
  from the forest by traversal — exactly what `assembleLeafGrid` does, with
  the block loop replaced by a forest walk. Coarse cells are depth-0 trees.

## 2. Restriction A — what may be refined *within* a cell

(The user's "restriction on some types of refinement in a cell.") Octant
geometry is the trilinear image of the root cell's 8 corners, which assumes a
clean hexahedron. So:

- **Refinable by octree**: clean cells — exactly 8 distinct corners and 6
  logical faces. Their octants are well-defined trilinear sub-boxes.
- **Restricted**: fault-split cells (>6 faces), pinched/degenerate cells
  (<8 distinct corners). Options, in order of effort:
  1. **Disallow refinement** of such cells (mark them non-refinable; a
     refinement request that hits one is rejected or clamped). Cheapest, and
     fine as a first cut — faults are a small fraction of cells.
  2. **Preprocessor refinement** of the single cell (resample its COORD/ZCORN
     box and run `findconnections`), as the static builder already does for a
     fault *inside* a block. More general, more expensive, only where needed.
- This is a real simplification, not a limitation in disguise: dynamic AMR
  typically tracks a moving front in the clean interior, away from faults.

## 3. Restriction B — refinement level *between* cells (2:1 balance)

(The user's "restriction on level of refinement between cells.") Adopt the
standard octree **2:1 balance**: two leaves sharing a face/edge/corner differ
by at most one refinement level. Enforced by a balance pass after marking
(refine the coarser side until satisfied). It buys:

- **Edge-conformalization becomes local and cheap**: at most one hanging node
  per coarse edge per refined neighbour. The built `edgeConformalizeLeaf`
  still works for arbitrary hanging nodes, but 2:1 bounds its work and makes
  it incremental.
- **Bounded interface set**: face/edge transitions are a small enumerable set
  (the coarse face sees a 1- or 2- or k-subdivision mosaic), so the
  leaf-builder's mosaic/face-pairing logic (already written for box
  boundaries) covers them without open-ended cases.
- **Stable transmissibility coupling**: each coarse-fine interface is a clean
  mosaic, matched in pillar parameter space (review §8.1) — refinement-
  invariant, computed once per macro face.

Faults at a 2:1 interface are still the hard case and fall under Restriction
A (don't refine fault-adjacent cells, first cut).

## 4. Identity — construction-stable packed ids (now load-bearing)

A leaf cell's global id = `pack(root cell's level-0 global id, octant path)`,
where the path is the Morton code of the octant within its tree. Faces/
corners: `pack(owning entity id, local lattice code)` with a fixed owner rule
(lowest adjacent root id). Properties:

- **Rank-independent by construction**: an octant owned by rank A and its
  ghost copy on rank B compute the *same* id (it is a function of the
  level-0 id — already globally consistent — and a deterministic local code).
  So the parallel `cell_index_set_` links owner and ghost with **no
  communication, no prediction, no winner-selection** — the machinery that
  dominated the original parallel-LGR code disappears.
- **Persistent under adaptation**: refining/coarsening changes only the
  affected octants' ids; everything else is stable (good for restart, for
  reusing the linear-solver setup).
- 64-bit packing suffices (level-0 id in the high bits, bounded-depth path in
  the low bits); fall back to a deterministic hash if a tree is pathologically
  deep.

This is the D3 scheme from DESIGN-builder, deferred in the static builder
(which used the existing IdSet delegation) but **required** here.

## 5. Parallelism — the octree manages its own ghost layer

The session's key parallel lesson: CpGrid's **cell-overlap scatter is not a
reliable foundation for refinement at rank boundaries** (overlap-1 corner/edge
gap; the ≥4-rank contracted-region `computeFace2Cell` fault). The static
builder sidesteps this with the *rank-interior* restriction (boxes never touch
the overlap). A dynamic octree cannot — the refinement front crosses ranks.

Design: **the forest owns a ghost layer of trees**, p4est-style, rather than
relying on CpGrid to re-scatter refined cells.

- Partition the *roots* (level-0 cells) across ranks — by a space-filling
  curve over (i,j,k), or reuse the existing graph partition. Each rank owns a
  set of root trees plus a **ghost layer**: the neighbour roots (one macro
  layer) replicated read-only.
- A rank refines its owned trees. For each ghost root, it **replicates the
  owner's tree** (exchange the compact per-tree refinement structure — bytes,
  not geometry). Because refinement is deterministic and the corner-point
  input is global, each rank then *reconstructs* the ghost octant geometry
  locally — **no geometry/topology communication**, only the small tree
  structure.
- Build the leaf view per rank from owned + ghost trees. Octant partition type
  = owner of its root tree (interior if owned, overlap/ghost if a ghost root).
  This is the same "inherit partition type from parent/root" rule already
  implemented (`PartitionTypeIndicator` + the leaf `cell_indicator_` setter).
- `cell_index_set_` and interfaces are built from the packed ids (§4): owner
  and ghost octants share ids, so the index set is correct without sync. The
  off-rank sentinel handling in the leaf builder (this session) already covers
  faces to non-replicated cells.
- **Rebalancing after adaptation**: migrate whole root trees between ranks
  (move the compact tree structure; geometry reconstructed locally). SFC
  ordering keeps this cheap and locality-preserving.

Net: the only true communication is (a) exchanging compact tree structures for
the ghost layer and after rebalancing, and (b) the normal solution-field halo
exchange — both small. This is the payoff of deterministic refinement + a
global macro description + construction-stable ids.

## 6. Dynamic adaptation cycle

1. **Mark** leaves (refine/coarsen) from an indicator (saturation front, etc.).
2. **Restrict** marks: drop refinement of non-refinable cells (Restriction A);
   apply 2:1 **balance** (Restriction B); communicate marks on the ghost layer
   so owners and ghosts agree.
3. **Apply** to the forest: grow/prune subtrees (coarsening = delete subtree —
   free, unlike the static builder).
4. **Rebuild** the leaf view by forest traversal (reuse `assembleLeafGrid`
   generalised), then `edgeConformalizeLeaf` if VEM.
5. **Transfer state**: restrict/prolong primary variables across the
   old→new leaf using the parent/child maps (conservative restriction;
   consistent prolongation) — the simulator-side half (review §7), unchanged
   by the grid design.

A correct first implementation can **rebuild the whole leaf** each adaptation
event (cost ∝ leaf size, the forest state is tiny); incremental leaf updates
are an optimisation, not a correctness requirement.

## 7. What transfers directly from this session

- `assembleLeafGrid` (block loop → forest walk), incl. the merged corner pool,
  mosaic face pairing, off-rank sentinel handling, leaf partition types.
- `edgeConformalizeLeaf` — unchanged.
- `GridStateWriter` composition seam — the octree backend is another
  `RefinementBuilder`.
- `assembleBlockLevelGrid`'s geometry/parent-relation logic — simplified to
  per-cell trilinear octants for clean cells.
- The rank-interior partition tooling (`setPartitionCellGroups`, overlap-2
  finding) informs the SFC root partition + ghost layer.
- The preprocessor-per-cell path for fault-adjacent cells (Restriction A
  option 2).

## 8. Risks / open questions

- **CpGrid scatter at scale**: the ≥4-rank contracted-region fault must be
  resolved *or* bypassed by the octree-owned ghost layer (the design choice
  here). Bypassing is cleaner but means the leaf's parallel structures are
  built by the octree code, not CpGrid's scatter — a larger change to how the
  distributed leaf is assembled.
- **Macro-grid overlap vs octree ghost**: reconciling CpGrid's existing
  cell-overlap (used by flow's communication) with the octree ghost layer
  needs care — ideally the octree ghost *is* the overlap, expressed through
  the same index sets.
- **Fault-adjacent refinement**: deferred by Restriction A; lifting it needs
  the per-cell preprocessor path and 2:1 balance across a fault, which is the
  genuinely hard geometry.
- **Anisotropic / k-ary vs strict octree**: CARFIN uses anisotropic factors;
  the trees should be k-ary. 2:1 balance generalises but the bookkeeping grows.

## 9. Incremental path

1. **Octree state + leaf rebuild, serial**: per-cell k-ary trees over level
   zero; mark/balance/apply; rebuild leaf via the generalised assembler;
   coarsening. Reuse edge-conformal. Validate against the static builder
   (a uniform-depth forest must reproduce a CARFIN refinement bitwise).
2. **Packed ids (D3)** replacing delegation — serial first (must stay
   bitwise-equal), then they enable parallel.
3. **Parallel: SFC root partition + ghost tree layer**, owned+ghost leaf
   build, ids-based index sets. 2-rank then 4-rank torture (the cases where
   CpGrid's own scatter failed) — this is where the octree-owned ghost layer
   earns its keep.
4. **Dynamic loop in flow** (`adaptGrid` hook) with conservative state
   transfer; front-tracking deck as the demonstrator.

Steps 1–2 are mostly assembling existing pieces; step 3 is the real new
parallel work; step 4 is the simulator-side half.
