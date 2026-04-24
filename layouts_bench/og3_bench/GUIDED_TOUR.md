# Guided tour of the OG3 C layout benchmark

This is a bottom-up walk through the implementation. It follows the code, not the abstract design. The central question is: when the global particle order changes, what does each layout have to move so that all arrays still describe the same particles?

## 0. The executable model

Open `Makefile` first. The important idea is that the benchmark driver is compiled once per layout and linked against exactly one layout object:

```text
bench_reshuffle_A  = bench_reshuffle.o + bench_common.o + hilbert.o + layout_A.o
bench_reshuffle_Bp = bench_reshuffle.o + bench_common.o + hilbert.o + layout_Bp.o
...
```

All four layout objects export the same function names from `layout_api.h`. Therefore `bench_reshuffle.c` does not know whether it is testing A, B', Bc, or C. This prevents accidental benchmark contamination: the same driver code calls the same API, and the layout object determines the storage behavior.

The three benchmark families are:

```text
bench_reshuffle_X   global PH ordering or sparse global permutation
bench_partition_X   same-type swap primitive
bench_sizes_X       actual sizeof/accounting report
```

where `X` is `A`, `Bp`, `Bc`, or `C`.

## 1. `bench_common.h`: the shared vocabulary

`bench_common.h` defines all scalar aliases and the particle-index type:

```c
typedef double   pos_t;
typedef double   acc_t;
typedef float    real_t;
typedef uint64_t pkey_t;
typedef uint64_t ppid_t;
typedef int64_t  time_int_t;
typedef uint64_t count_t;
```

The important one is `count_t`. It is used for:

```text
particle counts
particle indices
permutation entries
type-local indices
swap pair indices
```

The particle type enum is also here:

```c
typedef enum {
    PT_GAS  = 0,
    PT_DM1  = 1,
    PT_DM2  = 2,
    PT_DM3  = 3,
    PT_STAR = 4,
    PT_BH   = 5,
    NTYPES  = 6
} part_type_t;
```

`pcount_t` stores only four counts:

```c
typedef struct {
    count_t n_gas;
    count_t n_dm;
    count_t n_star;
    count_t n_bh;
} pcount_t;
```

DM subtypes are represented in the `types[]` array, but the allocation count is a single total `n_dm`.

## 2. Alignment and allocation

Every field atom uses:

```c
#define OG3_ATOM __attribute__((aligned(OG3_ATOM_ALIGN)))
```

with default:

```c
#define OG3_ATOM_ALIGN 64
```

All dynamic allocations go through:

```c
static inline void *og3_aligned_alloc(size_t bytes)
```

This rounds byte sizes up to a multiple of `OG3_ALLOC_ALIGN` and calls C11 `aligned_alloc`. The layouts allocate each stream independently. There is no arena in this C benchmark.

Why this matters: a struct that has 48 raw bytes may have `sizeof(T) == 64`, and `PDyn` has `sizeof(PDyn) == 128` under the default alignment. The benchmark measures actual `sizeof`, not raw field sums.

## 3. Synthetic population: `counts_from_total`, `gen_keys`, `assign_types`

Every benchmark begins from the same synthetic data pipeline.

### 3.1 Counts

`counts_from_total(N, &c)` computes:

```text
n_dm   = N * 4949 / 10000
n_star = N *   50 / 10000
n_bh   = N *    1 / 10000
n_gas  = N - n_dm - n_star - n_bh
```

So for `N=1000000`, the layout receives:

```text
500000 gas, 494900 DM, 5000 stars, 100 BH
```

### 3.2 Keys

`gen_keys(keys, N, n_cells, seed)` creates random integer coordinates:

```c
i = rng_next() % n_cells;
j = rng_next() % n_cells;
k = rng_next() % n_cells;
```

Then:

```c
keys[p] = ph_key_from_ijk(i, j, k);
```

`ph_key_from_ijk()` is just a wrapper around the real Peano-Hilbert implementation:

```c
return hilbert_key_3d_ijk(i, j, k, 21u);
```

The hard-coded `21u` means 21 bits per dimension, or 63 effective key bits.

### 3.3 Types

`assign_types(types, &c, seed)` first fills the array in blocks:

```text
all gas, then DM1, then DM2, then DM3, then stars, then BH
```

Then it shuffles the whole `types[]` array using Fisher-Yates. The result is a mixed global type sequence, not a type-grouped one.

### 3.4 IDs

The drivers set:

```c
ids[p] = (ppid_t)p;
```

This ID is the semantic identity used to verify that layout reshuffling did not misalign common fields and payload fields.

## 4. Peano-Hilbert key implementation: `hilbert.c`

`hilbert.c` contains:

```c
pkey_t hilbert_key_3d_ijk(uint32_t x, uint32_t y, uint32_t z, unsigned bits)
```

The routine:

1. clamps `bits` to at most 21;
2. masks the input coordinates to `[0, 2^bits)`;
3. starts with `rotation = 0` and `key = 0`;
4. walks from the most significant coordinate bit down to the least significant;
5. builds a 3-bit octant code `pix` from x/y/z at that bit;
6. appends `subpix3[rotation][pix]` to the key;
7. updates `rotation = rottable3[rotation][pix]`.

This is the real PH key path used by the dense reshuffle benchmark. The earlier Morton proxy is gone.

## 5. The permutation convention

The whole codebase uses one convention:

```text
perm[i] = old source slot for new destination slot i
```

Equivalently:

```text
after reshuffle, array[i] contains old array[perm[i]]
```

The byte helper, all layout reshuffles, and all derived sub-permutations follow this convention.

Keep this in mind when reading layout code. The loop usually looks like:

```c
for (count_t i = 0; i < N; ++i) {
    count_t src = perm[i];
    ... use old slot src to fill new slot i ...
}
```

## 6. Stable radix sort: `radix_sort_u64_count`

Dense reshuffle needs a global PH ordering. `bench_reshuffle.c` builds it in `build_global_key_perm()`:

```c
for (count_t i = 0; i < N; ++i) iden[i] = i;
radix_sort_u64_count(keys, iden, sorted_keys, perm_out, N);
```

The key array is `keys[]`. The payload is `iden[]`, so after sorting the payload is the permutation from new sorted position to old input position.

`radix_sort_u64_count()` is a stable LSD radix sort over four 16-bit passes:

```text
pass 0: bits  0..15
pass 1: bits 16..31
pass 2: bits 32..47
pass 3: bits 48..63
```

For each pass:

1. count bucket occupancy;
2. prefix-sum counts into bucket starts;
3. scan source arrays in order and scatter to destination arrays.

Because every pass is stable and the input payload is the original index order, equal keys keep the original order. That gives deterministic tie-breaking without storing a separate `(key, original_idx)` struct.

## 7. Dense and sparse byte permutation: `og3_apply_perm_bytes`

Every layout eventually calls:

```c
og3_apply_perm_bytes(array, perm, N, sizeof(T), scratch, scratch_idx);
```

The helper first constructs the list of changed destinations:

```c
if (perm[i] != i) changed_idx[changed++] = i;
```

If nothing changed, it returns.

### 7.1 Dense path

If `changed * 4 >= N`, at least 25 percent of slots changed. The helper uses full gather plus copy-back:

```c
for i:
    memcpy(tmp + i * elem_size,
           src + perm[i] * elem_size,
           elem_size);
memcpy(src, tmp, N * elem_size);
```

This is a simple streaming path. Accounting uses:

```text
4 * N * elem_size
```

### 7.2 Sparse path

If fewer than 25 percent of slots changed, it moves only changed destinations:

```c
for each changed destination dst:
    tmp[m] = src[perm[dst]];
for each changed destination dst:
    src[dst] = tmp[m];
```

This is safe for cycles because the first pass reads every needed old source value before the second pass writes any destination.

Accounting uses:

```text
2 * changed * elem_size
```

The layouts do not separately implement dense/sparse logic. They only provide arrays and permutations.

## 8. `fields.h`: the fine atoms used by B', Bc, and C

`fields.h` is where the benchmark decides what data exists. The shared fine atoms are:

```text
Common: PCore, PDyn, PTime, PMeta, PLinkage
Gas:    GasCore, GasGrad, GasMetal, GasSF
Star:   StarCore, StarMeta
BH:     BHCore, BHEnv, BHRepos
```

A few examples:

```c
typedef struct OG3_ATOM {
    pos_t      pos[3];
    real_t     mass;
    pkey_t     key;
    uint8_t    type;
    uint8_t    flags;
    int16_t    time_bin;
    int32_t    leaf_idx;
} PCore;
```

```c
typedef struct OG3_ATOM {
    count_t    type_idx;
} PLinkage;
```

```c
typedef struct OG3_ATOM {
    real_t     metals[16];
} GasMetal;
```

Layout A has its own fat structs, but their common prefix and type-specific suffixes mirror this atom decomposition.

## 9. `layout_api.h`: what every layout must implement

The API is intentionally small:

```c
layout_alloc
layout_free
layout_fill
layout_reshuffle_full
layout_swap_same_type
layout_reshuffle_bytes
layout_reshuffle_streams
layout_swap_bytes
layout_swap_streams
layout_get_key
layout_set_key
layout_get_pos
layout_set_pos
layout_get_type
layout_get_id
layout_verify_deep
```

The benchmark drivers only call this API. They never access `struct layout_ctx` fields directly.

This is the C equivalent of testing all candidate layouts under one external interface.

## 10. Layout A tour: fat arrays plus logical slots

Open `layout_A.c`.

### 10.1 What exists in memory

A owns four payload arrays:

```c
GasFullA  *gas;
DMFullA   *dm;
StarFullA *star;
BHFullA   *bh;
```

Each struct contains common fields plus any type-specific fields. For example, `GasFullA` starts with `A_COMMON_FIELDS`, then appends gas density, hydro acceleration, metals, SFR/wind fields, and the DEBUG marker field.

A also owns the logical global order:

```c
SlotA *slots;
SlotA *slot_scratch;
```

`SlotA` is:

```c
typedef struct {
    uint8_t type;
    count_t type_idx;
} SlotA;
```

So global slot `i` means:

```text
slots[i].type tells which of gas/dm/star/bh owns the payload
slots[i].type_idx tells which element inside that per-type array
```

### 10.2 How fill builds it

`layout_fill()` walks synthetic slot `p` once.

For gas:

```c
ctx->slots[p] = (SlotA){ .type = t, .type_idx = ig };
fill_common_prefix(&ctx->gas[ig], k, t, i, p);
ctx->gas[ig].density = ...;
ctx->gas[ig].metals[0] = ...;
++ig;
```

For DM, star, and BH the pattern is identical, with the relevant counter.

The result is:

```text
slots[0..N)       logical global order
gas[0..n_gas)     gas payloads in encounter order
dm[0..n_dm)       DM payloads in encounter order
star[0..n_star)   star payloads in encounter order
bh[0..n_bh)       BH payloads in encounter order
```

### 10.3 How A reshuffles

The driver passes one global permutation. A cannot call `og3_apply_perm_bytes()` on one global payload array, because no such array exists. Instead:

```c
SlotA old = ctx->slots[perm[i]];
```

For every new global slot `i`, A looks up which old type-specific payload should appear there. It appends `old.type_idx` to the appropriate sub-permutation and rewrites `slot_scratch[i]` so that the new global slot points to the new type-local position.

Example for gas:

```c
sub_gas[kg] = old.type_idx;
ctx->slot_scratch[i] = (SlotA){ .type = old.type, .type_idx = kg++ };
```

After the pass, A physically reshuffles each fat array:

```c
og3_apply_perm_bytes(ctx->gas,  sub_gas,  n_gas,  sizeof(GasFullA),  ...);
og3_apply_perm_bytes(ctx->dm,   sub_dm,   n_dm,   sizeof(DMFullA),   ...);
og3_apply_perm_bytes(ctx->star, sub_star, n_star, sizeof(StarFullA), ...);
og3_apply_perm_bytes(ctx->bh,   sub_bh,   n_bh,   sizeof(BHFullA),   ...);
```

Then it installs the new global directory:

```c
memcpy(ctx->slots, ctx->slot_scratch, N * sizeof(SlotA));
```

### 10.4 How A answers accessors

`layout_get_key(ctx, i)` does:

```c
SlotA s = ctx->slots[i];
switch (s.type) {
    case PT_GAS:  return ctx->gas[s.type_idx].key;
    case PT_DM*:  return ctx->dm[s.type_idx].key;
    case PT_STAR: return ctx->star[s.type_idx].key;
    case PT_BH:   return ctx->bh[s.type_idx].key;
}
```

`layout_set_key()` uses the same dispatch but writes the selected fat-struct field.

`layout_get_pos()` and `layout_set_pos()` also dispatch through `SlotA`, then read or write the `pos[3]` prefix in the selected fat struct. This is the important asymmetry relative to B'/Bc/C: A has no `core[i]` common array, so even a common-field update needs the logical-slot directory.

`layout_get_id()` is identical in shape. `layout_get_type()` can return `ctx->slots[i].type` directly.

### 10.5 How A swaps

For same-type swap, A translates both global slots through `SlotA` and swaps exactly one fat struct in one per-type array. It does not alter `slots[]`.

This works because slot `i` still points to the same type-local element, and that element now contains the other particle's payload.

## 11. Layout B' tour: fine arrays plus cross-index

Open `layout_Bp.c`.

### 11.1 What exists in memory

B' owns five common arrays of length `N`:

```c
PCore    *core;
PDyn     *dyn;
PTime    *time;
PMeta    *meta;
PLinkage *linkage;
```

It also owns fine type-specific arrays:

```c
GasCore, GasGrad, GasMetal, GasSF
StarCore, StarMeta
BHCore, BHEnv, BHRepos
```

The connection is explicit:

```text
common slot i -> linkage[i].type_idx -> type-specific slot
```

The type-specific array chosen depends on `core[i].type`.

### 11.2 How fill builds it

For every input slot `p`, B' fills common arrays at the same index `p`:

```c
ctx->core[p].key  = keys[p];
ctx->core[p].type = types[p];
ctx->meta[p].id   = ids[p];
```

If the particle is gas:

```c
ctx->linkage[p].type_idx = ig;
ctx->gas_core[ig].density = ...;
ctx->gas_metal[ig].metals[0] = ...;
++ig;
```

Stars and BHs follow the same append-and-link pattern. DM particles set `type_idx = 0` and have no payload array.

### 11.3 How B' reshuffles

B' reshuffles only the common streams:

```c
og3_apply_perm_bytes(ctx->core,    perm, N, sizeof(PCore),    ...);
og3_apply_perm_bytes(ctx->dyn,     perm, N, sizeof(PDyn),     ...);
og3_apply_perm_bytes(ctx->time,    perm, N, sizeof(PTime),    ...);
og3_apply_perm_bytes(ctx->meta,    perm, N, sizeof(PMeta),    ...);
og3_apply_perm_bytes(ctx->linkage, perm, N, sizeof(PLinkage), ...);
```

The type-specific payload arrays are not moved. This is correct because `PLinkage` moves with the common particle.

After reshuffle:

```text
core[i], dyn[i], time[i], meta[i], linkage[i]
```

all describe the same particle, and `linkage[i].type_idx` still points to the correct old payload slot.

### 11.4 How B' updates common fields through the API

For B', `layout_get_key()` and `layout_set_key()` read/write `ctx->core[i].key`. `layout_get_pos()` and `layout_set_pos()` read/write `ctx->core[i].pos[0..2]`. The type-specific arrays are not involved because key and position are common fields.

### 11.5 How B' swaps

`layout_swap_same_type()` swaps the five common arrays at `i` and `j`. It does not touch the type-specific arrays. Again, the linkage value moves with the common particle, so the payload association remains valid.

## 12. Layout Bc tour: coarse arrays plus cross-index

Open `layout_Bc.c`.

### 12.1 What exists in memory

Bc is the coarse version of B'. It owns four common arrays:

```c
PCoreB    *core;
PDynB     *dyn;
PAuxB     *aux;
PLinkageB *linkage;
```

`PAuxB` combines the identity and time-related fields that B' stores in `PMeta` and `PTime`.

Bc owns one composite payload array per type:

```c
GasAllB  *gas_all;   // GasCore + GasGrad + GasMetal + GasSF
StarAllB *star_all;  // StarCore + StarMeta
BHAllB   *bh_all;    // BHCore + BHEnv + BHRepos
```

### 12.2 How fill builds it

The common fill is similar to B':

```c
ctx->core[p].key  = keys[p];
ctx->core[p].type = t;
ctx->aux[p].id    = ids[p];
```

For gas:

```c
ctx->linkage[p].type_idx = ig;
ctx->gas_all[ig].core.density = ...;
ctx->gas_all[ig].metal.metals[0] = ...;
++ig;
```

The star and BH cases write `star_all[is]` and `bh_all[ib]`.

### 12.3 How Bc reshuffles

Bc only reshuffles:

```c
core, dyn, aux, linkage
```

The large composite arrays stay in place. Since `PLinkageB` moves with the common data, the payload association survives.

This is the cheapest cross-index reshuffle in stream count: four streams instead of B' five streams.

### 12.4 Scratch sizing

`layout_alloc()` computes:

```c
max_common = max(sizeof(PCoreB), sizeof(PDynB), sizeof(PAuxB), sizeof(PLinkageB));
ctx->scratch_bytes = N * max_common;
```

It intentionally ignores `GasAllB`, `StarAllB`, and `BHAllB`, because global reshuffle never calls the byte permutation helper on those arrays.

### 12.5 How Bc updates common fields through the API

For Bc, `layout_get_key()` and `layout_set_key()` read/write `ctx->core[i].key`. `layout_get_pos()` and `layout_set_pos()` read/write `ctx->core[i].pos[0..2]`. The identity accessor is different from B' and C because ID lives in `ctx->aux[i].id` rather than `PMeta`.

### 12.6 How Bc swaps

Same-type swap swaps the four common streams. It does not touch `gas_all`, `star_all`, or `bh_all`. The payload association moves through `PLinkageB`.

## 13. Layout C tour: fine positional layout with resident `type_idx[]`

Open `layout_C.c`.

### 13.1 What exists in memory

C owns four common arrays:

```c
PCore *core;
PDyn  *dyn;
PTime *time;
PMeta *meta;
```

It owns the same fine type-specific arrays as B':

```c
GasCore, GasGrad, GasMetal, GasSF
StarCore, StarMeta
BHCore, BHEnv, BHRepos
```

It also owns bookkeeping arrays:

```c
count_t *type_idx;
count_t *sub_perm_gas;
count_t *sub_perm_star;
count_t *sub_perm_bh;
```

In this C benchmark, `type_idx[p]` is resident. It tells where the payload for global common slot `p` lives in the relevant type-specific arrays.

### 13.2 How `type_idx[]` is built

`refresh_type_idx()` scans the current common order:

```c
count_t counters[NTYPES] = {0};
for (count_t p = 0; p < N; ++p) {
    uint8_t t = ctx->core[p].type;
    ctx->type_idx[p] = counters[t]++;
}
```

So the first gas in common order has `type_idx = 0`, the second gas has `type_idx = 1`, and so on. Same for stars and BHs. DM particles also receive counters, but there is no DM payload array.

### 13.3 How fill builds it

`layout_fill()` writes common arrays at global index `p`, appends gas/star/BH payloads in encounter order, and then calls `refresh_type_idx(ctx)`.

After fill, the type-specific arrays and the common order agree through the implicit ordinal relation stored in `type_idx[]`.

### 13.4 How C derives sub-permutations

C must reshuffle common arrays and type-specific arrays. Therefore, before moving anything, it walks the global permutation:

```c
for (count_t i = 0; i < N; ++i) {
    count_t src = perm[i];
    uint8_t t = ctx->core[src].type;
    count_t k = ctx->type_idx[src];
    append k to sub_perm_gas/star/bh according to t;
}
```

This says: if the new global order sees gas payloads in the order old gas indices `[7, 2, 9, ...]`, then every gas array must be permuted by `[7, 2, 9, ...]`.

### 13.5 How C updates common fields through the API

For C, key and position are common fields, so `layout_set_key()` writes `ctx->core[i].key`, while `layout_get_pos()` / `layout_set_pos()` read and write `ctx->core[i].pos[0..2]`. The resident `type_idx[]` is not touched by these calls because changing key or position does not change the particle type or its type-local ordinal.

### 13.6 How C reshuffles

After deriving sub-permutations, C permutes the common arrays:

```c
core, dyn, time, meta
```

Then it permutes every type-specific stream:

```text
gas_core, gas_grad, gas_metal, gas_sf      with sub_perm_gas
star_core, star_meta                       with sub_perm_star
bh_core, bh_env, bh_repos                  with sub_perm_bh
```

Finally it calls `refresh_type_idx(ctx)` again.

The refresh is not optional in this benchmark. The common type sequence has changed, so the mapping from common slot to type-local ordinal must be recomputed.

### 13.7 How C swaps

For same-type swap, C swaps the common arrays and also swaps the matching payload slots:

```c
uint8_t t = ctx->core[i].type;
count_t ti = ctx->type_idx[i];
count_t tj = ctx->type_idx[j];

swap common arrays at i and j;
swap type arrays at ti and tj;
```

It does not change `type_idx[]`. This is correct because the payloads at `ti` and `tj` were swapped to follow the common particles.

## 14. Deep verification in DEBUG builds

The benchmark has two levels of verification.

The normal reshuffle verifier in `bench_reshuffle.c` checks:

```text
layout_get_id(ctx, i)   == ids[perm[i]]
layout_get_type(ctx, i) == types[perm[i]]
```

Dense mode also checks global key order:

```text
layout_get_key(ctx, i-1) <= layout_get_key(ctx, i)
```

Before the first reshuffle, `bench_reshuffle.c` also exercises the key/position access API by temporarily changing and restoring a few slots:

```text
layout_set_key(ctx, i, new_key)
layout_set_pos(ctx, i, new_pos)
layout_get_key(ctx, i) == new_key
layout_get_pos(ctx, i) == new_pos
```

Those checks do not prove that payload arrays are aligned. That is why each layout implements `layout_verify_deep()` under `#ifdef DEBUG`.

### 14.1 B' and Bc deep verification

For every global slot:

```text
read common type and ID
read linkage/type_idx
read the corresponding payload marker
compare marker against ID-derived expected value
```

A stale linkage fails here.

### 14.2 C deep verification

For every global slot:

```text
read common type and ID
read resident type_idx[p]
read the corresponding payload marker
compare marker against ID-derived expected value
```

A wrong sub-permutation or stale `type_idx[]` fails here.

### 14.3 A deep verification

For every global slot:

```text
read SlotA{type,type_idx}
read the selected fat payload
check type and marker
```

A broken `SlotA[]` rewrite or wrong per-type sub-permutation fails here.

In non-DEBUG builds all of this compiles to a stub returning success.

## 15. Reshuffle benchmark walk-through

Open `bench_reshuffle.c`.

### 15.1 Parse arguments

```c
N        = argv[1] or 1000000
n_iter   = argv[2] or 10
fraction = argv[3] or 1.0
```

`fraction >= 0.9999` selects dense PH-sort mode. Smaller fractions select sparse global permutation mode.

### 15.2 Generate inputs

```c
counts_from_total(N, &c);
gen_keys(keys, N, 1u << 21, 12345ULL);
assign_types(types, &c, 67890ULL);
ids[p] = p;
```

### 15.3 Build permutation

Dense:

```c
build_global_key_perm(keys, N, perm, sorted_keys);
verify_sort(sorted_keys, N);
```

Sparse:

```c
gen_sparse_permutation(perm, N, fraction, 98765ULL);
```

Then:

```c
changed = count_permutation_changed(perm, N);
```

### 15.4 Allocate, fill, warm up

```c
ctx = layout_alloc(&c);
layout_fill(ctx, keys, types, ids);
layout_reshuffle_full(ctx, perm);
```

The warmup is also the correctness target: after it, the layout should contain exactly the particles described by `perm`.

### 15.5 Verify

```c
verify_against_perm(ctx, perm, types, ids, N);
if (dense) verify_global_sorted(ctx, N);
layout_verify_deep(ctx);
```

### 15.6 Time

```c
for (it = 0; it < n_iter; ++it)
    layout_reshuffle_full(ctx, perm);
```

The loop measures the cost of the movement primitive. It does not re-generate keys or re-sort between iterations.

### 15.7 Report

The report combines layout accounting and measured CPU time:

```text
streams permuted
payload copy bytes
perm-index bytes
sparse payload estimate, if sparse
ms / reshuffle
ns / particle
throughput estimate
```

The timing source is `cpu_time_p()`, which wraps `clock_gettime(CLOCK_PROCESS_CPUTIME_ID, ...)`.

## 16. Partition benchmark walk-through

Open `bench_partition.c`.

Despite the name, this is not yet the real box-leaf repartition algorithm. It is a primitive benchmark for same-type swaps.

The flow is:

```text
parse N, fraction, n_iter
generate keys/types/ids
allocate and fill layout
n_pairs = fraction * N
gen_swap_pairs_same_type(pairs, n_pairs, types, N, seed)
sum layout_swap_bytes(ctx, type) over all pairs
warmup: apply all pairs once
DEBUG verify
time n_iter passes over all pairs
print per-type swap cost and timing
```

The important helper is `gen_swap_pairs_same_type()`. It builds per-type pools of global indices, then picks random pairs from pools that have at least two members. Therefore the layout function can assume:

```c
layout_get_type(ctx, i) == layout_get_type(ctx, j)
```

The layouts use this guarantee differently:

```text
A:  swap one fat payload struct through SlotA
B': swap five common streams, including PLinkage
Bc: swap four common streams, including PLinkageB
C:  swap four common streams plus type-specific payload streams
```

## 17. Size benchmark walk-through

Open `bench_sizes.c`.

This program exists because comments and raw field sums are not enough. Alignment changes the actual sizes.

The flow is:

```text
parse N
compute counts
print layout name and description
print shared atom sizeof values
allocate layout
print layout_reshuffle_streams
print layout_reshuffle_bytes
print layout_swap_bytes for gas/dm/star/bh
free layout
```

Use this after changing `fields.h`, `OG3_ATOM_ALIGN`, or any layout struct. It is the quickest way to see whether a field move changed the benchmark's memory traffic.

## 18. One particle through all layouts

Assume input slot `p` is a gas particle with:

```text
key = K
id  = I
type = PT_GAS
```

### In A

```text
slots[p] = { PT_GAS, ig }
gas[ig].key = K
gas[ig].id  = I
gas[ig].density = ...
```

Global access always goes through `slots[p]`.

### In B'

```text
core[p].key = K
core[p].type = PT_GAS
meta[p].id = I
linkage[p].type_idx = ig
gas_core[ig].density = ...
gas_metal[ig].metals[0] = ...
```

Global access reads common arrays directly and payload through `linkage[p]`.

### In Bc

```text
core[p].key = K
core[p].type = PT_GAS
aux[p].id = I
linkage[p].type_idx = ig
gas_all[ig].core.density = ...
gas_all[ig].metal.metals[0] = ...
```

Same cross-index idea as B', but payload fields are packed into one composite struct.

### In C

```text
core[p].key = K
core[p].type = PT_GAS
meta[p].id = I
gas_core[ig].density = ...
gas_metal[ig].metals[0] = ...
```

Then after `refresh_type_idx()`:

```text
type_idx[p] = ordinal of this gas particle in the current common order
```

Global access uses `type_idx[p]` to find payload.

## 19. Why the reshuffle costs differ

The same global permutation has different physical consequences.

### A

A moves four fat arrays, but only within each type. It also copies the global `SlotA[]` directory.

Cost character:

```text
few streams, fat elements, one extra global directory
```

### B'

B' moves all common fine streams and linkage. It does not move payload arrays.

Cost character:

```text
five global streams, fine common atoms, cheap mutation bookkeeping
```

### Bc

Bc moves four common coarse streams and linkage. It does not move payload composites.

Cost character:

```text
fewest cross-index reshuffle streams, coarser common atoms
```

### C

C moves common fine streams and all type-specific fine streams. It also refreshes `type_idx[]`.

Cost character:

```text
many streams, no cross-linkage stream, positional consistency must be rebuilt
```

This is why `bench_sizes_C` reports 13 reshuffle streams under the current field set.

## 20. What is not implemented yet

This package does not yet implement the tree-aware operations that should be the next benchmark layer.

Missing pieces:

```text
real box-leaf builder integration
leaf descriptors with [start,end) and [gas_start,gas_end)
between-rebuild repartition by leaf ID
single-rank distribute.c-style local pack/unpack benchmark
OpenMP parallel loops
MPI exchange
GPU radix/tree code
```

The current reshuffle benchmark is now globally PH-sorted, so it is the right first primitive for the tree rebuild. The current partition benchmark is still only a swap primitive; it should be replaced or supplemented once the leaf builder is connected.

## 21. Suggested reading order for changes

For a correctness change to reshuffling:

```text
bench_common.c: og3_apply_perm_bytes
bench_reshuffle.c: permutation construction and verification
layout_C.c: hardest reshuffle case
layout_A.c: most unusual global-to-per-type translation
layout_Bp.c and layout_Bc.c: cross-index reference cases
```

For a field-layout change:

```text
fields.h
layout_A.c fat structs
layout_Bc.c composite structs
bench_sizes.c output
layout_verify_deep marker paths
```

For future leaf partition work:

```text
layout_swap_same_type in all four layouts
bench_partition.c
then add a new bench_leaf_partition.c rather than overloading the primitive benchmark
```
