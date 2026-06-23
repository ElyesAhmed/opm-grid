# Design: a fully parallel, locally (dynamically) refinable CpGrid

This updates review §8.4 (forest-of-octrees) with what the static + parallel
rebuild on this branch actually proved. The earlier §8.4 was a sketch; this
is a concrete design whose building blocks are now mostly implemented and
tested. Read alongside [DESIGN-builder.md](DESIGN-builder.md) (state model,
shared corner pool, ids) and [PLAN.md](PLAN.md) (what is built).

**Name:** the class designed here is **`AdaptiveCpGrid`** — heir to `CpGrid` and
to the review's "forest-of-octrees" sketch. "the octree" is used loosely below to
mean `AdaptiveCpGrid`. It is a **new, additive** class: it must not modify or
regress `CpGrid` or the existing corner-point LGR path (those tests stay green);
it reuses the shared kernels (`GrdeclRefinement`, `assembleLeafGrid`,
`edgeConformalizeLeaf`) without changing their behaviour. Three design decisions
were fixed 2026-06-23 — scope (§2), split policy (§13), and AMR backend (§14,
deferred) — see those sections.

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
  stores its per-direction split factors. Per the §13 policy these are a
  (possibly anisotropic, arbitrary) factor on a root's *first* split and
  factor-2-anisotropic (each of i/j/k split by 2 or not) on every level below.
  Leaves of the forest are the simulation cells. This *replaces* the accumulating
  `CpGridData` levels with compact per-cell data (see §10 for the full
  representation).
- **Leaf view** = a materialized `CpGridData` (the DUNE-facing grid) rebuilt
  from the forest by traversal — exactly what `assembleLeafGrid` does, with
  the block loop replaced by a forest walk. Coarse cells are depth-0 trees.

## 2. Restriction A — what may be refined *within* a cell

(The user's "restriction on some types of refinement in a cell.") Octant
geometry is the trilinear image of the root cell's 8 corners, which assumes a
clean hexahedron.

**Scope (fixed 2026-06-23): always refined from a corner-point grid — never
general polyhedral.** Two things that are easy to conflate must be kept apart:

- The **refinable parent** must be a clean hex: exactly **8 distinct corners**,
  so its sub-boxes are well-defined (and resampled on sub-pillars, §10/§13). This
  is the actual restriction.
- The **resulting leaf cells** are still hex *geometry* (8 corners on their
  pillars) but may carry **more than 6 faces** — a logical side abutting finer
  neighbours is split into a conformal sub-face mosaic (§10). They are
  topologically multi-face but never *general* polyhedra (no PEBI/Voronoi). The
  corner-point sub-pillar structure is exactly what keeps the interface cheap: a
  coarse↔fine face is a 2-D clip in pillar parameter space, not 3-D polygon
  clipping. (The old phrasing "clean cells = 8 corners *and* 6 faces" wrongly
  pinned leaf cells at 6 faces — corrected here.)

So:

- **Refinable**: clean-hex parents (8 distinct corners).
- **Restricted**: fault-split parents (non-hex / >6 *logical* faces) and
  pinched/degenerate parents (<8 distinct corners). Options, in order of effort:
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
(refine the coarser side until satisfied). Note this "2:1" presumes the
factor-2 dynamic levels of §13; a root's arbitrary *first* split is the root's
intrinsic subdivision, not a balance boundary (a refined root and an unrefined
neighbour still differ by one level). It buys:

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
where the path encodes the octant within its tree. Per §13 that path is a
**first-split radix code** (the root's arbitrary first factor) followed by
**fixed bits per factor-2 level** (Morton-clean after level 1), so many levels
fit in 64 bits; fall back to a deterministic hash only if a tree is
pathologically deep. Faces/corners: `pack(owning entity id, local lattice code)`
with a fixed owner rule (lowest adjacent root id). Properties:

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

Full-leaf rebuild each adaptation event (cost ∝ leaf size, the forest state is
tiny) is the **correctness oracle** — the first thing to implement and the
reference every later optimisation diffs against. But it is *not* the intended
steady-state path: fast AMR refines a moving front every few timesteps, so
**fast incremental local refinement is a first-class capability**, designed in
from the start and built/validated in serial first. How the leaf vectors are
mutated locally (append + free-list + stable-id identity + stencil-only edits +
lazy compaction) is §12.

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
- **Anisotropic / k-ary vs strict octree**: *resolved* — see §13. Policy is
  factor-2 anisotropic with an arbitrary first split, which keeps 2:1 balance and
  id packing cheap while preserving CARFIN-exact and directional refinement.
- **Forest/AMR backend (build vs reuse)**: *deferred, gated on a spike* — see
  §14. p4est/t8code conflict with §13's anisotropy; in-house forest is the
  working assumption pending the serial spike (§9.1).

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

## 10. Representation: compact persistent state vs derived per-rank leaf

The representation must **not** be full vectors over the global grid, yet must
still represent the multi-face (hex-geometry) leaf cells of §2. This resolves
into two layers; the pleasant surprise is that the multi-face cells need **no new
container**.

**Layer A — persistent / migrated / serialized state (compact, small):**
- **Macro corner-point description**: `COORD/ZCORN/ACTNUM` + dims. This already
  exists as `RetainedCornerPointInput` (retained at construction when the deck
  has LGRs, broadcast to every rank in `CpGrid::addLgrsUpdateLeafView`). It *is*
  level 0 / the forest connectivity.
- **Forest**: one compact tree per root cell (per-node split factors + leaf flag,
  per §13; nothing stored for unrefined roots).

Layer A is the whole source of truth — what you migrate for rebalancing (move
the tree *bytes*) and serialize for restart. No leaf vectors live here.

**Layer B — derived leaf (per-rank, owned+ghost only; never global):**
- **Topology**: CpGrid's `cell_to_face_` is a variable-length `SparseTable`
  (`OrientedEntityTable<0,1>` over `SparseTable`, `OrientedEntityTable.hpp:138`).
  A hex with conformal sub-faces is simply a **longer row** — the "polyhedral"
  cell falls out of the existing CSR topology. The leaf assembler's only job is
  to *generate* the extra faces: for each of a cell's 6 logical sides, consult
  the neighbour root's tree leaves abutting that side; under 2:1 balance that is
  1 face or a bounded k×k mosaic, produced by the parameter-space clip.
- **Geometry**: corners reconstructed from sub-pillars (`root pillars + octant
  code`); the D1 shared corner pool becomes a per-rank cache computed on demand,
  not persisted or communicated (consistent with §5).
- **Ids**: packed D3 (§4).

**Why the global "full vectors" anti-pattern is avoided:** the only full-length
arrays that ever exist are these Layer-B leaf caches, bounded by the local
partition (owned+ghost), never the global grid. Global state stays Layer A
(forest + macro), both compact.

**Consequences:**
- *Rebalance* = migrate tree bytes; the macro is already everywhere → rebuild the
  leaf locally, no leaf-vector communication.
- *Restart* = serialize macro + forest, rebuild deterministically — the real fix
  path for `LGR_GAPS` B3 ("refined-grid restart rejected").
- *Adapt* = grow/prune subtrees (compact), rebuild only the affected local leaf
  (incrementally, §12).

## 11. Review note (2026-06-23): on reducing all refinement to the current LGR kernels

Claim reviewed: "all the machinery for refining any given cell and finding the
intersection between two cells is OK, since all refinement can be achieved by the
current LGR framework."

**Verdict: agree in spirit.** The two kernels exist and are reusable — cell
subdivision is a 1×1×1 box via `assembleBlockLevelGrid`; the conformal interface
is the leaf-assembler mosaic + `edgeConformalizeLeaf`; `Intersection` is a thin
generic reader over `cell_to_face_`/`face_to_cell_` (`Intersection.cpp:33-46`).
So `AdaptiveCpGrid` is mostly assembly of built pieces (§0/§7). Two real caveats
(the earlier "general polyhedral" caveat is **withdrawn** — see §2: the model is
corner-point hex geometry with possibly >6 faces, which is a *simplification*
that makes intersections cheap, not a limitation):

1. **Refinable-parent geometry restriction** (clean hex; fault/pinch parents
   deferred — §2). Reservoir AMR often wants refinement *near* faults, so this is
   more load-bearing than "a small fraction of cells" suggests.
2. **Multi-level interfaces are designed, not proven** (§3): the mosaic is
   verified for a single coarse↔fine jump over a flat box face; general 2:1
   graded, mixed-level interfaces are the riskiest reuse — the nested-LGR
   multi-level leaf assembly was exactly the hard part.

**Overall design assessment:** sound and well-grounded; the load-bearing calls
are right (packed ids §4, octree-owned ghost §5, full-rebuild-as-oracle §6).
Risks to track: (a) mosaic face-pairing for general 2:1 interfaces needs a
genuinely mixed-level test, not just the uniform-depth CARFIN-bitwise check
(§9.1); (b) reconciling CpGrid's overlap (consumed by flow's communication, well
info, and the linear solver's `OwnerOverlapCopyCommunication`) with the octree
ghost layer (§8); (c) root-tree-granularity rebalancing wants a **weighted** SFC
(weight = leaf count per tree) to avoid per-tree imbalance; (d) conservative
black-oil state transfer (§6.5) is real critical-path work.

## 12. Fast (incremental) local refinement: mutating the leaf vectors

The "full-length vectors" are the Layer-B leaf `CpGridData` arrays —
`cell_to_face_`, `face_to_cell_`, `face_to_point_`, `cell_to_point_`,
`geometry_` (codim 0/1/3), `global_cell_`, `face_tag_`, `face_normals_`,
`unique_boundary_ids_`, and the index/id sets — each sized to an entity count and
addressed `[0..N)`.

**The obstacle is the dense numbering, not the size.** Refining one cell into k
children plus new faces/points, if indices must stay dense, means insert-in-
middle + an O(N) renumber that invalidates every index-keyed map. That O(N)-per-
event cost is why §6's baseline is full rebuild. Fast adapt avoids it:

1. **Append, never insert.** Tombstone the refined cell's slot (push to a free
   list); append its k children and the new faces/corners at the end of their
   arrays. Untouched indices never move → cost O(local stencil).
2. **Identity from D3 ids, not from the slot.** Because ids are
   `pack(root id, octant code)` (§4), an entity's identity is independent of its
   array position: `index→id` is a pure function, `id→index` a hash patched
   locally. So append/tombstone/compact need no renumbering. **D3 is therefore
   load-bearing for fast adapt, not only for parallel.**
3. **Edit only the stencil.** Replace each shared face with its conformal
   sub-face mosaic (1 or a bounded k×k under 2:1), set `face_to_cell_` for the
   new faces, edit each neighbour's `cell_to_face_` row in place, and reconstruct
   the new corners/faces from sub-pillars — all within the local neighbourhood.
4. **Coarsening = tombstone the children, revive the parent** from the free list.
   Symmetric and local.
5. **Lazy compaction.** Let tombstones accumulate; run one O(N) compaction +
   index-set rebuild only when fragmentation crosses a threshold (or at
   output/checkpoint). Amortized cost stays near O(touched cells).

**Container implication.** `SparseTable` is build-once
(begin/append/finalize), not in-place-growable. Fast adapt therefore needs either
an extended `SparseTable` (append/tombstone/lazy-compact) or — less invasive as a
first move — a **mutable id-keyed side adjacency** (per-cell small vector of
faces) for the active front, flattened into the leaf arrays on demand.

**Serial-first / locality.** Refining a cell is intrinsically rank-local (the
cell, its face-neighbours, and the new sub-entities). Parallel adds only (a) a
ghost tree-byte exchange for the changed roots and (b) a balance-mark exchange
across rank cuts. So build and validate fast incremental refine/coarsen **in
serial on one partition** first (diffing against full-rebuild, §6), then layer
the ghost exchange on top — matching §9 step 1.

## 13. Split policy (decided 2026-06-23)

**Anisotropic, per-direction, with an arbitrary first split.** A root's *first*
split may use an arbitrary anisotropic factor `(rx,ry,rz)` (so static CARFIN —
including odd factors like 3×3×3 — is reproduced exactly). Every *subsequent*
dynamic level is **factor-2 anisotropic**: each of i/j/k is independently split
by 2 or not (→ 2, 4, or 8 children). A 4× refinement is two binary levels; the
only thing excluded *dynamically* is odd factors beyond the first split. Full
general k-ary refinement remains available in the existing static builder.

Why this hybrid, vs the alternatives considered (strict isotropic 2×2×2; strict
factor-2 anisotropic; full general k-ary):
- Keeps **Morton-clean fixed-width child codes** after level 1 → D3 (§4) packs
  many levels in 64 bits.
- Keeps **textbook 2:1 balance** (§3) and a **bounded ≤2×2 interface mosaic** —
  the part flagged as the riskiest reuse.
- Preserves **CARFIN-exact** static refinement (the arbitrary first split) and
  the **directional** (vertical-only / areal-only) refinement reservoirs need
  (which strict isotropic octree cannot express).
- Geometry is factor-agnostic (sub-pillars handle any factor — `GrdeclRefinement`
  already does arbitrary `rx,ry,rz`), so the cost of generality lands only on
  ids/balance/mosaic, which the hybrid keeps cheap.

## 14. Forest / parallel-AMR backend: build vs reuse (deferred 2026-06-23)

Whether to build the forest/balance/ghost/repartition layer in-house or reuse an
existing AMR library was evaluated and **deferred** (gated on a spike, below).

| Option | Mature parallel AMR | Anisotropic factor-2 (§13) | Irregular CP macro + faults | Custom pillar geometry | Feeds existing `CpGridData` leaf + OPM stack |
|---|---|---|---|---|---|
| **p4est** | ✓ (gold standard) | ✗ isotropic 1→8 only | ~ (1 tree/level-0 cell via connectivity) | ✓ (geometry-agnostic) | ✗ must rebuild leaf + Dune/OPM index sets |
| **t8code** | ✓ + multi-element-type | ✗ schemes still isotropic | ~ | ✓ | ✗ same |
| **dune-alugrid** | ✓ (ParMETIS) | ✗ cube 1→8 | ✗ own grid, no CP semantics | ✗ | ✗ different grid type |
| **dune-spgrid** | ✓ structured, *does* anisotropic | ✓ but structured only | ✗ tensor-product, no faults | partial | ✗ |
| **AMReX / Chombo / SAMRAI** | ✓ block-structured | ~ (ratio 2/4) | ✗ logically Cartesian patches | ✗ | ✗ wrong paradigm |
| **In-house over CpGrid** | must build (borrow p4est algorithms) | ✓ exact | ✓ (level-0 `face_to_cell`) | ✓ (sub-pillars) | ✓ native (reuses `assembleLeafGrid`, `edgeConformalizeLeaf`, `GrdeclRefinement`) |

**Key finding:** p4est/t8code hard-code **isotropic 1→8** refinement, which
conflicts with §13 and with reservoir vertical-/areal-only needs; and even if
accepted they buy only the forest+ghost+partition *half* — the corner-point
leaf/geometry/OPM-index-set half is still in-house. So the backend choice is
**coupled to §13**: anisotropic ⇒ in-house forest; isotropic ⇒ p4est becomes
attractive.

**Decision: DEFERRED — gate on a spike.** Prototype the minimal in-house
**serial** forest (§9 step 1: per-cell trees over level zero; mark/balance/apply;
rebuild leaf via the generalised assembler; validate that a uniform-depth forest
reproduces a CARFIN refinement bitwise) before locking the backend. Reusable
blocks regardless of the final choice: **Zoltan/ParMETIS** (already integrated,
`zoltanGoG`) for the root partition, weighted by per-tree leaf count; a small
**SFC** header for ordering; **Dune comm** primitives for ghost/migration; and
the existing `assembleLeafGrid` / `edgeConformalizeLeaf` / `GrdeclRefinement`
kernels. Port p4est's 2:1-balance and ghost-construction *algorithms* by reading,
not linking. Revisit p4est/t8code only if the spike argues for isotropic
refinement.
