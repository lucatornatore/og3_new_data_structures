# Partition benchmark notes - v5

This version replaces the provisional v4 box-leaf builder with a standalone
module adapted from the uploaded `create_boxleaves` prototype.

## Particle generation

Particles are generated as 21-bit integer positions per axis. Keys are then
derived from those positions with the existing Peano-Hilbert routine.

Generation is implemented in `particle_init.[ch]` and supports two modes:

- `OG3_POS_PLAIN`: uniform random integer coordinates.
- `OG3_POS_CLUSTERED`: particles are sampled around random cluster centres.

The uploaded prototype used an OpenMP loop that modified the loop index and
wrote a run of future particles from inside one iteration. That is fine in a
serial loop, but unsafe in an OpenMP `for`: another thread can own those future
indices. The v5 helper keeps the same plain/clustered distinction but uses
independent per-particle RNG streams, so it is safe with or without OpenMP.

The runtime `temperature` is a coldness parameter in `[0,1]`: lower values make
clusters tighter, higher values make them broader.

## Box-leaf creation

Box leaves are built in `boxleaves.[ch]` from globally PH-sorted keys. The
parameters are:

- target particles per leaf;
- tolerance, converted into `[min_particles, max_particles]`;
- maximum side length in 21-bit PH coordinate cells.

The builder follows the uploaded algorithm's core idea: start from a target
particle bunch, choose a PH cube that contains it, refine if overfull, and emit
only non-empty PH cubes. Empty PH space is therefore not represented as leaves;
if a particle later moves into such a gap, the partition benchmark reports it
as `outside old leaves`, which is a rebuild signal.

## Small-N bug fixed

The uploaded implementation can fail when `N` is small compared with the target
leaf occupancy. The most direct problem is initial capacity:

```c
maxNbleaves = (Np / PperBOX) * 2;
```

If `Np < PperBOX`, this can be zero. The code then allocates zero leaves and
writes the first leaf through that pointer. v5 uses a minimum capacity and grows
by `realloc` as needed.

The OpenMP reconstruction path also contains accesses of the form `K[pstart+1]`.
Those need guarding when `pstart` reaches the last particle. The v5 builder has
no unguarded `pstart+1` access.

## Partition semantics

The benchmark now models leaf maintenance between rebuilds:

1. generate positions;
2. compute PH keys;
3. globally sort by PH key;
4. build box leaves;
5. fill the selected layout in sorted order;
6. perturb positions by `+/- perturb_cells` in each coordinate;
7. recompute keys;
8. classify each particle into its post-perturbation old box leaf;
9. counting-sort the layout by leaf;
10. keep gas first inside each leaf;
11. update leaf `[begin,end)` and `[gas_begin,gas_end)` ranges.

The final verification checks both leaf contiguity and the gas subrange inside
each leaf.
