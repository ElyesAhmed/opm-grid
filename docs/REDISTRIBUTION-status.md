# Status: distribution and redistribution, with/without LGR

Two distinct operations are often conflated:

- **Distribution** (`loadBalance`/`scatterGrid`): serial grid → distributed
  grid. One-time, at simulation start.
- **Redistribution / rebalancing**: re-partition a grid that is *already*
  distributed (to fix load imbalance). A separate capability.

This summarises both for the original upstream code and the rebuild on
`strip-lgr`, with and without LGR. Grounded in the `scatterGrid` source of
each (the relevant guards are identical between the two — the strip kept
`scatterGrid` intact).

## Redistribution (rebalancing an already-distributed grid)

**Not supported in any of the four cases.** `scatterGrid` begins with

```
if (!distributed_data_.empty()) {
    std::cerr << "There is already a distributed version of the grid...";
    return {false, {}};
}
```

so a second call on a distributed grid is a no-op (returns false). There is
no `redistribute`/`repartition` method anywhere. This is a **CpGrid-level
gap, identical in original and rebuild, independent of LGR**:

| | redistribute (rebalance) |
|---|---|
| original, no LGR | ✗ (refused) |
| original, with LGR | ✗ (refused; leaf-with-LGR also throws — below) |
| new, no LGR | ✗ (same code) |
| new, with LGR | ✗ (same code) |

## Distribution (one-time, serial → parallel)

`scatterGrid` additionally enforces (same in both codes):
- loadbalancing a refined *level* grid (`level > 0`) → throws;
- loadbalancing the *leaf* of a grid that already has LGRs (`level == -1`,
  `maxLevel() > 0`) → throws "not supported";
- loadbalancing *level zero* of an LGR'd grid → only `zoltanGoG`.

The consequence is that the only supported order is **distribute first, then
refine** — never "refine serially, then distribute the refined grid". Flow
follows this: `CpGridVanguard::loadBalance()` runs before `addLgrs()`.

| | distribute level 0 (no LGR yet) | then refine in parallel | distribute an already-refined grid |
|---|---|---|---|
| original, no LGR | ✓ | n/a | n/a |
| original, with LGR | ✓ | ✓ — general (boxes may cross ranks), via heavy parallel-LGR machinery (id prediction, overlap communication, winner-selection) | ✗ (leaf throws) |
| new, no LGR | ✓ (verified: SPE9 `np=2`) | n/a | n/a |
| new, with LGR | ✓ | ✓ — **rank-interior only** (a box must stay inside one rank, not touch the overlap), verified `np=2`/`np=4`; simpler, no id prediction | ✗ (leaf throws) |

So the two LGR implementations agree on *what is possible* (distribute then
refine; no rebalancing); they differ only in *how* the parallel refinement
works — the original is more general but heavy, the rebuild is restricted but
simple (and is the foundation the octree generalises).

## How this relates to the new locally-refinable (octree) class

Dynamic AMR **requires redistribution** — refining a region creates load
imbalance, so the partition must change between time steps. CpGrid provides
none: not for the original, not for the rebuild, with or without LGR. This is
not a minor gap to patch later; it is a primary reason the octree design
([DESIGN-parallel-octree.md](DESIGN-parallel-octree.md) §5) must **own
partitioning itself** rather than lean on `scatterGrid`:

- Partition the **root trees** (level-0 cells) by a space-filling curve.
- Rebalance by **migrating whole root trees** between ranks — move the
  compact per-tree refinement structure (bytes), then reconstruct octant
  geometry locally (deterministic refinement + global macro input ⇒ no
  geometry communication).
- The leaf view and its parallel index sets are built by the octree from
  owned + ghost trees, using construction-stable ids — *not* by re-scattering
  through CpGrid.

In other words, the octree's ghost-tree + tree-migration machinery is the
redistribution capability that CpGrid lacks, expressed at the root-tree
granularity where it is cheap and communication-free for geometry.

Corollary for the interim static path: the rebuild's parallel LGR (this
session) inherits the "balance once at start, never rebalance" limitation —
fine for static CARFIN runs, fatal for dynamic AMR. So "add redistribution"
and "build the octree" are the same project: there is little value in adding
generic CpGrid rebalancing for the static leaf when the dynamic class needs a
different (root-tree) redistribution anyway.

## One-line summary

Redistribution (rebalancing) is unsupported everywhere — original and
rebuild, with and without LGR — and is a CpGrid-level gap; both LGR
implementations only ever *distribute once then refine*; and providing real
redistribution is intrinsically part of the octree/dynamic-AMR class (via
root-tree migration), not a separate static-grid feature worth adding first.
