# OG3 C layout benchmark - implementation documentation

This document describes the current pure-C benchmark implementation in this directory. It is a code reference for the v2 draft, not a design note. The design question is the resident data layout for the OpenGadget3 rewrite; this package isolates that question by hard-coding four C layouts and measuring how much work they do when global particle order changes.

The four layouts are:

| Layout | Source file | Short description |
|---|---|---|
| A | `layout_A.c` | One fat AoS array per particle class, plus a logical global slot directory. |
| B' | `layout_Bp.c` | Fine common atoms, fine type-specific atoms, explicit `PLinkage` cross-index. |
| Bc | `layout_Bc.c` | Coarse common atoms, one composite payload struct per type, explicit `PLinkageB` cross-index. |
| C | `layout_C.c` | Fine common atoms, fine type-specific atoms, positional mapping, with resident `type_idx[]` kept in this benchmark. |

The benchmark is intentionally single-process and single-threaded. It does not contain MPI, OpenMP, CUDA, HIP, optional-physics templates, or the C++ `ArrayRegistry`. It is a controlled C harness for measuring the cost of maintaining layout invariants.

## 1. Source tree

The canonical tree is `layouts_bench/og3_bench/`:

```text
Makefile
README.md
VALIDATION.md
DOCUMENTATION.md        # this file
GUIDED_TOUR.md          # bottom-up code walk
bench_common.h
bench_common.c
fields.h
hilbert.h
hilbert.c
layout_api.h
layout_A.c
layout_Bp.c
layout_Bc.c
layout_C.c
bench_reshuffle.c
bench_partition.c
bench_sizes.c
```

The top-level duplicate files under `layouts_bench/` are only convenience copies. For development, read and edit `og3_bench/`.

## 2. Build model

`Makefile` builds 12 executables:

```text
bench_reshuffle_A    bench_reshuffle_Bp    bench_reshuffle_Bc    bench_reshuffle_C
bench_partition_A    bench_partition_Bp    bench_partition_Bc    bench_partition_C
bench_sizes_A        bench_sizes_Bp        bench_sizes_Bc        bench_sizes_C
```

Each executable links exactly one `layout_*.o` because all four layout files export the same C API symbols:

```c
layout_alloc
layout_fill
layout_reshuffle_full
layout_swap_same_type
layout_get_key
layout_set_key
layout_get_pos
layout_set_pos
layout_get_type
layout_get_id
layout_verify_deep
```

This is deliberate. The benchmark driver is identical across layouts; the link line chooses the layout under test.

Default build:

```bash
make clean && make all
```

Debug-marker build:

```bash
make clean && \
CFLAGS='-std=gnu11 -O0 -g -Wall -Wextra -DDEBUG -fno-strict-aliasing' \
make all
```

Pedantic build:

```bash
make clean && \
CFLAGS='-std=gnu11 -O2 -Wall -Wextra -Wpedantic -fno-strict-aliasing' \
make all
```

## 3. Shared scalar and index types

`bench_common.h` defines the benchmark-wide scalar types:

```c
typedef double   pos_t;
typedef double   acc_t;
typedef float    real_t;
typedef uint64_t pkey_t;
typedef uint64_t ppid_t;
typedef int64_t  time_int_t;
typedef uint64_t count_t;
```

`count_t` is used for particle counts, particle-array indices, permutations, type indices, and swap-pair indices. This removes the earlier artificial `uint32_t` indexing limit from the benchmark interface.

Particle types are:

```c
PT_GAS, PT_DM1, PT_DM2, PT_DM3, PT_STAR, PT_BH, NTYPES
```

Dark matter is represented by three DM subtypes because OG3 frequently carries several collisionless species. The macro `PT_IS_DM(t)` recognizes `PT_DM1`, `PT_DM2`, and `PT_DM3`.

## 4. Counts, synthetic keys, types, and IDs

Every benchmark starts by deriving a synthetic population with `counts_from_total(N, &c)` in `bench_common.c`:

```text
gas   = remainder, normally 50.00 percent
dm    = 49.49 percent
star  = 0.50 percent
bh    = 0.01 percent
```

For `N = 1000000`, this gives:

```text
gas=500000  dm=494900  star=5000  bh=100
```

`gen_keys()` generates random 3D integer coordinates and turns them into Peano-Hilbert keys:

```c
keys[p] = ph_key_from_ijk(i, j, k);
```

`ph_key_from_ijk()` calls:

```c
hilbert_key_3d_ijk(i, j, k, 21u);
```

The implementation is in `hilbert.c` / `hilbert.h`. It is a Martin Reinecke/Gadget-style 3D Peano-Hilbert table walk. Coordinates are masked to the requested number of bits, and the benchmark uses 21 bits per dimension, producing a 63-bit effective key.

`assign_types()` first fills a type array with the requested counts, then applies a Fisher-Yates shuffle using the local xoshiro-style PRNG. For small or empty runs, it returns before the unsigned reverse loop, so `N == 0` and `N == 1` do not underflow.

IDs are simply:

```c
ids[p] = (ppid_t)p;
```

The ID is the stable semantic identity used by the verification paths.

## 5. Field atoms in `fields.h`

Layouts B', Bc, and C share the fine atoms declared in `fields.h`. Layout A does not use these declarations directly; it inlines equivalent fields into its fat per-type structs.

Common fine atoms:

| Atom | Role |
|---|---|
| `PCore` | position, mass, PH key, type, flags, time bin, leaf index |
| `PDyn` | velocity, drift displacement, gravity acceleration, PM acceleration, neighbor counters |
| `PTime` | time integration fields and gravity-cost vector |
| `PMeta` | particle ID and auxiliary flags |
| `PLinkage` | `count_t type_idx`, used by cross-index layouts |

Gas atoms:

| Atom | Role |
|---|---|
| `GasCore` | hydro scalars such as density, pressure, entropy, temperature |
| `GasGrad` | predicted velocity, hydro acceleration, div/curl, energy step |
| `GasMetal` | `metals[16]` |
| `GasSF` | star-formation, wind, and BH-swallow bookkeeping |

Star atoms:

| Atom | Role |
|---|---|
| `StarCore` | stellar age, chemical-time metadata, PID, scalar weights |
| `StarMeta` | star metal vector |

BH atoms:

| Atom | Role |
|---|---|
| `BHCore` | BH mass, accretion, swallow ID, local density/entropy |
| `BHEnv` | surrounding cold/hot/gas phase environment |
| `BHRepos` | repositioning, accreted momentum, final spin, minimum-potential data |

The alignment macro is:

```c
#ifndef OG3_ATOM_ALIGN
#define OG3_ATOM_ALIGN 64
#endif
#define OG3_ATOM __attribute__((aligned(OG3_ATOM_ALIGN)))
```

With the default `OG3_ATOM_ALIGN=64`, `bench_sizes_*` reports these actual sizes:

| Atom | `sizeof` |
|---|---:|
| `PCore` | 64 |
| `PDyn` | 128 |
| `PTime` | 64 |
| `PMeta` | 64 |
| `PLinkage` | 64 |
| `GasCore` | 64 |
| `GasGrad` | 64 |
| `GasMetal` | 64 |
| `GasSF` | 64 |
| `StarCore` | 64 |
| `StarMeta` | 64 |
| `BHCore` | 64 |
| `BHEnv` | 64 |
| `BHRepos` | 128 |

Do not rely on raw sizes in comments for performance accounting. Use `bench_sizes_X`, because `sizeof(T)` changes with `OG3_ATOM_ALIGN`.

## 6. Public layout API

`layout_api.h` is the common ABI exported by all four layout files.

### Allocation and destruction

```c
layout_ctx_t *layout_alloc(const pcount_t *c);
void layout_free(layout_ctx_t *ctx);
```

Each layout owns a private `struct layout_ctx`. The benchmark never dereferences it directly.

### Fill

```c
void layout_fill(layout_ctx_t *ctx,
                 const pkey_t  *keys,
                 const uint8_t *types,
                 const ppid_t  *ids);
```

`layout_fill()` walks the synthetic input arrays in input order and constructs the resident layout. All layouts set the benchmark-visible fields `key`, `type`, and `id`; they also write a few representative payload fields such as gas density, star age, and BH mass.

In `DEBUG` builds, type-specific marker fields are filled from the ID:

| Type | Marker field |
|---|---|
| gas | `GasSF.swallow_id`, or `GasFullA.swallow_id`, or `GasAllB.sf.swallow_id` |
| star | `StarCore.pid`, or `StarFullA.pid_star`, or `StarAllB.core.pid` |
| BH | `BHCore.swallow_id`, or `BHFullA.swallow_id_bh`, or `BHAllB.core.swallow_id` |

These markers allow `layout_verify_deep()` to prove that common arrays and type-specific arrays still describe the same particle after a reshuffle or swap.

### Full reshuffle

```c
void layout_reshuffle_full(layout_ctx_t *ctx, const count_t *perm);
```

The permutation convention is gather-style:

```text
after reshuffle: new slot i contains old slot perm[i]
```

This convention is used everywhere: global reshuffle, per-type sub-permutations, and sparse permutations.

### Same-type swap primitive

```c
void layout_swap_same_type(layout_ctx_t *ctx, count_t i, count_t j);
```

The caller guarantees that global slots `i` and `j` have the same type. This benchmark is a primitive for future leaf-partition repair; it is not yet the full box-leaf repartition algorithm.

### Accounting and accessors

```c
size_t layout_reshuffle_bytes(const layout_ctx_t *ctx);
int    layout_reshuffle_streams(const layout_ctx_t *ctx);

size_t layout_swap_bytes(const layout_ctx_t *ctx, uint8_t t);
int    layout_swap_streams(const layout_ctx_t *ctx, uint8_t t);

pkey_t  layout_get_key (const layout_ctx_t *ctx, count_t i);
void    layout_set_key (layout_ctx_t *ctx, count_t i, pkey_t key);

void    layout_get_pos (const layout_ctx_t *ctx, count_t i, pos_t out[3]);
void    layout_set_pos (layout_ctx_t *ctx, count_t i, const pos_t in[3]);

uint8_t layout_get_type(const layout_ctx_t *ctx, count_t i);
ppid_t  layout_get_id  (const layout_ctx_t *ctx, count_t i);
```

The key and position mutators are deliberately part of the layout API because their physical storage differs by layout. In B', Bc, and C they write the common `core[i]` stream. In A they translate global slot `i` through `SlotA` and update the owning fat per-type struct. `bench_reshuffle.c` exercises these accessors immediately after `layout_fill()` with a round-trip set/restore self-test.

`layout_reshuffle_bytes()` reports an estimate of actual dense-copy payload traffic: source read, scratch write, scratch read, destination write. Index-array traffic is estimated separately by `bench_reshuffle.c`.

## 7. Byte-level permutation helper

All layouts call `og3_apply_perm_bytes()` from `bench_common.c`:

```c
void og3_apply_perm_bytes(void *array,
                          const count_t *perm,
                          count_t N,
                          size_t elem_size,
                          void *tmp_bytes,
                          count_t *changed_idx);
```

It chooses between two paths.

### Dense path

If at least 25 percent of slots move, the helper uses a streaming gather plus copy-back:

```text
for i in [0,N): tmp[i] = array[perm[i]]
memcpy(array, tmp, N * elem_size)
```

Estimated payload traffic:

```text
4 * N * elem_size
```

This represents read old array, write scratch, read scratch, write destination.

### Sparse path

If fewer than 25 percent of slots move, the helper copies only changed destinations:

```text
changed_idx[] = all i where perm[i] != i
for m in changed_idx: tmp[m] = array[perm[m]]
for m in changed_idx: array[m] = tmp[m]
```

Estimated payload traffic:

```text
2 * changed * elem_size
```

The first pass reads all sources before any destination write occurs, so cycles and swaps are safe.

## 8. Layout A - fat per-type AoS plus `SlotA[]`

`layout_A.c` declares four per-type fat structs:

```c
GasFullA
DMFullA
StarFullA
BHFullA
```

Each begins with the same common prefix, enforced by `_Static_assert` checks on key offsets such as `offsetof(..., id)`. The helper `fill_common_prefix()` writes the common prefix through a temporary `GasFullA`, then copies the prefix into the actual destination struct.

The resident context is:

```c
struct layout_ctx {
    pcount_t    counts;
    GasFullA   *gas;
    DMFullA    *dm;
    StarFullA  *star;
    BHFullA    *bh;
    SlotA      *slots;
    SlotA      *slot_scratch;
    count_t    *scratch;
    count_t    *scratch_idx;
    void       *scratch_data;
    size_t      scratch_data_bytes;
};
```

`SlotA` is:

```c
typedef struct {
    uint8_t type;
    count_t type_idx;
} SlotA;
```

It is the logical global array. Global slot `i` does not directly contain payload. Instead:

```text
slots[i].type      says which per-type array owns the particle
slots[i].type_idx  says which element inside that per-type array
```

### A fill

`layout_fill()` walks input slot `p` and appends the particle to the correct fat array:

```text
gas input  -> gas[ig],   slots[p] = {PT_GAS,  ig++}
dm input   -> dm[idm],   slots[p] = {PT_DMx,  idm++}
star input -> star[is],  slots[p] = {PT_STAR, is++}
bh input   -> bh[ib],    slots[p] = {PT_BH,   ib++}
```

The per-type arrays are in encounter order after fill. The only global ordering is `SlotA[]`.

### A reshuffle

A global PH permutation cannot be directly applied to one common array, because A has no common array. `layout_reshuffle_full()` therefore does two things in one pass over the global permutation:

1. Builds per-type sub-permutations:

```text
sub_gas[]
sub_dm[]
sub_star[]
sub_bh[]
```

2. Builds the new logical global directory in `slot_scratch[]`.

For every new global slot `i`:

```c
SlotA old = ctx->slots[perm[i]];
```

If `old.type == PT_GAS`, then `old.type_idx` is appended to `sub_gas[]`, and `slot_scratch[i]` is rewritten to point at the new gas-local index produced by that append. The same rule applies to DM, stars, and BHs.

Then A applies four per-type byte permutations:

```text
gas[]  <- sub_gas
dm[]   <- sub_dm
star[] <- sub_star
bh[]   <- sub_bh
```

Finally:

```c
memcpy(ctx->slots, ctx->slot_scratch, N * sizeof(SlotA));
```

### A swap

`layout_swap_same_type()` translates global slots through `SlotA`, asserts the two types match, and swaps one fat struct in one per-type array. `SlotA[]` does not need to change: slot `i` still points at the same type-local index, but the payload at that index has been exchanged with the payload at slot `j`'s type-local index.

## 9. Layout B' - fine atoms plus `PLinkage`

`layout_Bp.c` is the fine-grained cross-index layout.

The resident context contains five common streams of length `N`:

```c
PCore    *core;
PDyn     *dyn;
PTime    *time;
PMeta    *meta;
PLinkage *linkage;
```

and type-specific fine arrays:

```c
GasCore  *gas_core;
GasGrad  *gas_grad;
GasMetal *gas_metal;
GasSF    *gas_sf;

StarCore *star_core;
StarMeta *star_meta;

BHCore   *bh_core;
BHEnv    *bh_env;
BHRepos  *bh_repos;
```

The connector is:

```c
ctx->linkage[p].type_idx
```

For global common slot `p`, the type-specific slot is interpreted according to `ctx->core[p].type`:

```text
PT_GAS  -> gas_*[linkage[p].type_idx]
PT_STAR -> star_*[linkage[p].type_idx]
PT_BH   -> bh_*[linkage[p].type_idx]
PT_DM*  -> no type-specific payload
```

### B' fill

`layout_fill()` writes all common streams at index `p`. For gas/star/BH it also appends payload to the corresponding type arrays and stores the append index in `PLinkage[p].type_idx`.

### B' reshuffle

`layout_reshuffle_full()` applies the same global permutation to exactly the five common streams:

```c
core, dyn, time, meta, linkage
```

The type-specific arrays do not move. Because `linkage` moves together with `core`, the particle's pointer to its payload follows the particle.

### B' swap

Same-type swap exchanges the five common streams, including `PLinkage`. Type-specific payload arrays do not move. The particle identity follows the swapped linkage value.

## 10. Layout Bc - coarse atoms plus `PLinkageB`

`layout_Bc.c` has the same connection scheme as B', but fewer, wider arrays.

The common streams are:

```c
PCoreB    *core;
PDynB     *dyn;
PAuxB     *aux;
PLinkageB *linkage;
```

The type-specific streams are composites:

```c
GasAllB  *gas_all;   // GasCore + GasGrad + GasMetal + GasSF
StarAllB *star_all;  // StarCore + StarMeta
BHAllB   *bh_all;    // BHCore + BHEnv + BHRepos
```

`PLinkageB.type_idx` has the same meaning as `PLinkage.type_idx` in B'.

### Bc fill

`layout_fill()` writes the common streams at global index `p`, appends the type-specific payload to `gas_all`, `star_all`, or `bh_all`, and stores the type-local append index in `PLinkageB[p].type_idx`.

### Bc reshuffle

`layout_reshuffle_full()` applies the global permutation to the four common streams only:

```c
core, dyn, aux, linkage
```

Type-specific composites do not move. This is the main intended bookkeeping advantage of Bc: it has one fewer common stream than B' and does not need to touch the large composite payload arrays during global PH reshuffle.

`layout_alloc()` sizes scratch only for the largest common stream:

```c
max(sizeof(PCoreB), sizeof(PDynB), sizeof(PAuxB), sizeof(PLinkageB)) * N
```

This is intentional; earlier versions over-sized Bc scratch by considering the large type-specific composites even though reshuffle does not permute them.

### Bc swap

Same-type swap exchanges the four common streams, including `PLinkageB`. The composite payload arrays stay still.

## 11. Layout C - fine atoms plus positional mapping and `type_idx[]`

`layout_C.c` uses the same fine atoms as B' but no `PLinkage` stream. Its common streams are:

```c
PCore *core;
PDyn  *dyn;
PTime *time;
PMeta *meta;
```

Its type-specific arrays are the same fine arrays as B':

```c
gas_core, gas_grad, gas_metal, gas_sf
star_core, star_meta
bh_core, bh_env, bh_repos
```

This benchmark intentionally keeps:

```c
count_t *type_idx;
```

`type_idx[p]` is the current type-local ordinal implied by the common stream. It is refreshed by scanning the common `core[p].type` sequence:

```c
ctx->type_idx[p] = counters[type]++;
```

In the future leaf-based form, this positional information should mostly come from leaf descriptors such as `(start,end)` and `(gas_start,gas_end)`. In this temporary C benchmark, the resident `type_idx[]` is kept so that same-type swaps and deep verification are simple and explicit. Its refresh traffic is included in C's accounting.

### C fill

`layout_fill()` writes the four common streams at global index `p`, appends type-specific payload to the type arrays in encounter order, then calls `refresh_type_idx(ctx)`. After fill, `type_idx[p]` tells where the payload for global slot `p` lives inside the relevant type-specific arrays.

### C reshuffle

C must move both common streams and type-specific streams. The order must remain mutually consistent.

First, C derives sub-permutations from the global permutation and the old `type_idx[]`:

```c
for i in new global order:
    src = perm[i]
    t   = core[src].type
    k   = type_idx[src]
    append k to sub_perm_gas/star/bh according to t
```

Then it applies the global permutation to common streams:

```c
core, dyn, time, meta
```

Then it applies the derived sub-permutations to each type-specific stream:

```text
gas_core, gas_grad, gas_metal, gas_sf      <- sub_perm_gas
star_core, star_meta                       <- sub_perm_star
bh_core, bh_env, bh_repos                  <- sub_perm_bh
```

Finally it calls `refresh_type_idx(ctx)` so each new common slot again knows its current type-local ordinal.

### C swap

For a same-type swap of global slots `i` and `j`, C reads:

```c
t  = core[i].type;
ti = type_idx[i];
tj = type_idx[j];
```

It swaps the four common streams at `i` and `j`. Then it swaps the corresponding type-specific payload streams at `ti` and `tj`. It does not swap `type_idx[]`; this is correct because slot `i` continues to point to type-local position `ti`, and the payload at `ti` has been swapped to match the common data now at slot `i`.

## 12. Reshuffle benchmark

Executable form:

```bash
./bench_reshuffle_X [N] [n_iter] [fraction]
```

where `X` is `A`, `Bp`, `Bc`, or `C`.

Default:

```text
N=1000000, n_iter=10, fraction=1.0
```

### Dense mode

If `fraction >= 0.9999`, dense mode runs. The driver:

1. Generates synthetic keys, types, and IDs.
2. Builds a global PH permutation with `build_global_key_perm()`.
3. Verifies that the sorted key array is nondecreasing.
4. Allocates and fills the layout.
5. Applies one warmup reshuffle.
6. Verifies semantic identity against the permutation.
7. Verifies global PH order.
8. Calls `layout_verify_deep()`.
9. Times `n_iter` repeated calls to `layout_reshuffle_full()`.

`build_global_key_perm()` creates `iden[i] = i` and calls the stable radix sorter:

```c
radix_sort_u64_count(keys, iden, sorted_keys, perm_out, N);
```

Because the radix sort is stable, equal PH keys keep their original input order. The original index is therefore the implicit deterministic tie-break.

### Sparse mode

If `fraction < 0.9999`, sparse mode uses:

```c
gen_sparse_permutation(perm, N, fraction, seed)
```

The generated permutation is mostly identity. Only a selected subset of global destinations is shuffled. The layout code is unchanged; the byte-level helper detects the sparse permutation and uses the changed-destination copy path.

Sparse mode does not check global PH sortedness, because the generated sparse permutation is not a PH sort. It checks semantic alignment and, in `DEBUG`, deep payload consistency.

### Output interpretation

Important output fields:

```text
streams permuted        : number of payload streams passed to og3_apply_perm_bytes
payload copy bytes      : dense 4-pass payload traffic estimate
perm-index bytes        : streams * N * sizeof(count_t), coarse upper estimate
sparse payload estimate : 2 * changed * payload, scaled from dense estimate
ms / reshuffle          : measured CPU time per layout_reshuffle_full call
ns / particle           : measured CPU time per particle
throughput estimate     : payload + index estimate divided by time
```

The timed loop repeatedly applies the same permutation to the already reshuffled layout. That is acceptable for timing the movement primitive; the semantic verification is performed immediately after the warmup application.

## 13. Partition benchmark

Executable form:

```bash
./bench_partition_X [N] [fraction] [n_iter]
```

Default:

```text
N=1000000, fraction=0.01, n_iter=20
```

This is currently a same-type swap primitive, not the final leaf repartition benchmark.

The driver:

1. Generates synthetic keys, types, and IDs.
2. Allocates and fills the layout.
3. Builds `n_pairs = fraction * N` random same-type swap pairs using `gen_swap_pairs_same_type()`.
4. Computes per-iteration byte estimates by summing `layout_swap_bytes(ctx, type)` over the generated pairs.
5. Performs one warmup pass over the pairs.
6. Runs `layout_verify_deep()`.
7. Times `n_iter` repeated passes over the same pairs.

The key point is that every layout receives identical global slot pairs. The layout implementation decides what has to move to keep its invariant valid.

## 14. Size benchmark

Executable form:

```bash
./bench_sizes_X [N]
```

It prints:

1. actual `sizeof()` values for shared fine atoms;
2. layout description;
3. particle counts derived from `N`;
4. reshuffle streams;
5. dense payload-copy bytes;
6. swap bytes per type.

With default alignment and `N=1000000`, the current accounting is:

| Layout | Reshuffle streams | Dense payload-copy bytes | Gas swap | DM swap | Star swap | BH swap |
|---|---:|---:|---:|---:|---:|---:|
| A | 4 | 1,425,356,800 | 2,688 | 1,536 | 1,920 | 2,688 |
| B' | 5 | 1,536,000,000 | 2,304 | 2,304 | 2,304 | 2,304 |
| Bc | 4 | 1,280,000,000 | 1,920 | 1,920 | 1,920 | 1,920 |
| C | 13 | 1,802,662,400 | 3,456 | 1,920 | 2,688 | 3,456 |

These are accounting numbers, not performance measurements. Use the timed benchmarks for wall clock.

## 15. DEBUG verification

`layout_verify_deep()` is a zero-cost stub in non-DEBUG builds:

```c
#ifndef DEBUG
return 0;
#endif
```

In `DEBUG` builds, it checks the invariant that the common slot and type-specific payload still belong to the same ID.

For cross-index layouts B' and Bc:

```text
common slot i -> linkage[i].type_idx -> type-specific payload
payload marker must match common ID
```

For C:

```text
common slot i -> type_idx[i] -> type-specific payload
payload marker must match common ID
```

For A:

```text
global slot i -> SlotA{type,type_idx} -> per-type fat payload
payload marker must match payload ID/type
```

This catches precisely the class of errors that a simple sorted-key check would miss: stale linkage, stale positional indices, or type-specific arrays permuted with the wrong sub-permutation.

## 16. Current limitations

The current code intentionally does not yet include:

1. MPI exchange.
2. OpenMP parallel loops.
3. Real box-leaf construction.
4. Leaf descriptor storage.
5. Between-rebuild leaf repartition by leaf ID.
6. A local `distribute.c`-style pack/unpack benchmark.
7. GPU radix sort or GPU tree build.

The next implementation step should be a tree-aware benchmark that consumes the actual leaf builder and records at least:

```text
leaf particle range:        [start, end)
hydro/gas subrange:         [gas_start, gas_end)
```

The future generalization can store more type-specific subranges, but the first useful version only needs the global leaf range plus one gas range for hydro locality.

## 17. Relation to the C++ prototype

This C benchmark intentionally removes the C++ complexity from the resident-layout question. It does not implement compile-time physics policies, constrained accessors, optional fields, arenas, or the C++ `ArrayRegistry`. Those mechanisms remain relevant to the full rewrite, but here the point is to measure the low-level cost of the four resident layouts under identical synthetic operations.

The conceptual mapping is:

| C benchmark | C++ prototype concept |
|---|---|
| `layout_Bp.c` | layout B' fine + cross-index |
| `layout_Bc.c` | layout B coarse + cross-index |
| `layout_C.c` | layout C fine + positional mapping |
| `og3_apply_perm_bytes()` | low-level equivalent of applying one permutation to one registered array |
| `layout_reshuffle_full()` | layout-specific version of reshuffling all arrays that must remain aligned |
| `layout_verify_deep()` | hard invariant check missing from ordinary sorted-key verification |

## 18. Reading checklist

To understand or modify the implementation, read in this order:

```text
1. bench_common.h          types, permutation convention, helper declarations
2. bench_common.c          PRNG, PH-key generation, radix sort, sparse/dense perm helper
3. fields.h                fine atoms shared by B', Bc, C
4. layout_api.h            common layout ABI
5. layout_A.c              fat per-type AoS plus SlotA directory
6. layout_Bp.c             fine cross-index layout
7. layout_Bc.c             coarse cross-index layout
8. layout_C.c              fine positional layout with resident type_idx[]
9. bench_reshuffle.c       global PH sort and reshuffle benchmark
10. bench_partition.c      same-type swap primitive
11. bench_sizes.c          size/accounting utility
```

## v5 partition and box-leaf update

The partition benchmark is no longer a same-type swap primitive. It now starts
from integer particle positions, derives Peano-Hilbert keys, globally sorts the
particles, builds box leaves, perturbs positions, recomputes keys, and partitions
by post-perturbation leaf membership. The code path is in `bench_partition.c`.

`particle_init.[ch]` owns the position generator. It supports plain uniform mode
and clustered mode. Clustered mode is based on the uploaded generator's intent
but avoids its OpenMP race pattern by assigning each particle independently from
thread-local/per-particle RNG streams.

`boxleaves.[ch]` owns the box-leaf builder. The builder is adapted from the
uploaded `create_boxleaves` implementation but fixes the small-N capacity bug and
unguarded `pstart+1` boundary risks. Each `boxleaf_t` stores the PH key interval,
the current common-array `[begin,end)` range, and a gas subrange
`[gas_begin,gas_end)`. During partitioning, gas particles are placed first inside
each leaf so hydro-relevant entries remain locally contiguous.
