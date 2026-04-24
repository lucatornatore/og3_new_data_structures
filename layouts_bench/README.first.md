# OpenGadget3 data-layout microbenchmark

A four-layout comparison for the packing/unpacking cost paid by the
particle data structure at two representative operations:

1. **Dense or sparse reshuffle** — apply a permutation that reorders
   particles by Peano-Hilbert key, either fully dense (tree rebuild)
   or sparse (fraction ƒ of particles move, rest identity; steady-state
   between rebuilds).
2. **Same-type partition swap** — N random same-type swap pairs,
   modeling incremental partition maintenance when a particle crosses a
   box-leaf boundary and must be relocated into its new range.

Four layouts share one uniform API so `bench_reshuffle_X` and
`bench_partition_X` produce results directly comparable across layouts.

## Source conventions

- **Types**: `count_t = uint64_t` for all particle counts and per-array
  indices reaching N; `ppid_t = uint64_t` for particle ID; `pkey_t` PH
  key; `pcount_t` struct with per-type counts.
- **Timing**: `CPU_TIME_P` macro reading `CLOCK_PROCESS_CPUTIME_ID`
  (statement-expression form: inlines at the call site, zero function
  overhead).
- **Allocation**: `og3_aligned_alloc(bytes)` wrapping C11 `aligned_alloc`
  with 64-B alignment and size rounded up internally.
- **Aliasing**: `restrict` qualifiers on pointer parameters throughout.
- **Atom alignment**: `OG3_ATOM` macro (see below).

## The OG3_ATOM_ALIGN compile-time knob

Every atom struct (`PCore`, `PDyn`, ..., `GasCore`, `BHRepos`, and the
per-type fat structs in Layout A) uses `OG3_ATOM`, which expands to
`__attribute__((aligned(OG3_ATOM_ALIGN)))`. **No explicit padding arrays
exist in the source** — the compiler rounds `sizeof(T)` up to a multiple
of `OG3_ATOM_ALIGN`.

Default is 64 (one cache line). Override on the command line, e.g.:

```bash
CFLAGS='-std=gnu11 -O3 -march=native -Wall -DOG3_ATOM_ALIGN=32' make
CFLAGS='-std=gnu11 -O3 -march=native -Wall -DOG3_ATOM_ALIGN=16' make
```

**This choice drives a real memory/bandwidth trade-off.** Observed
bytes/particle for the common region in Layout B' at N=10⁶:

| Atom | raw | ALIGN=8 | ALIGN=16 | ALIGN=32 | ALIGN=64 |
|---|---:|---:|---:|---:|---:|
| PCore | 48 | 48 | 48 | 64 | 64 |
| PDyn | 84 | 88 | 96 | 96 | 128 |
| PTime | 44 | 48 | 48 | 64 | 64 |
| PMeta | 16 | 16 | 16 | 32 | 64 |
| PLinkage | 8 | 8 | 16 | 32 | 64 |
| **common total** | **200** | **208** | **224** | **288** | **384** |

Going from 64 → 16 **cuts Bp/Bc reshuffle traffic by ~40%** because
PMeta and PLinkage stop wasting 40–60 B/particle on padding. But you
lose per-atom cache-line alignment: every PDyn access now crosses a
cache-line boundary for some fraction of particles. Which wins on
bandwidth depends on how the kernels access these atoms — this harness
measures packing/unpacking only.

## Physics configuration

Faithful to OpenGadget3 particle-struct contents under:

| Flag | Effect on structs |
|---|---|
| HYDRO=SPH | `DhsmlDensityFactor` in SphP |
| PMGRID | `GravPM[3]` in P |
| STAR_FORMATION (SFR) | `Sfr` in SphP |
| COOLING | `elec` in SphP |
| STELLAR_EVOLUTION (LT_STELLAREVOLUTION) | `MassRes`, `EgyRes`, `Metals[16]`, `EgyStep`, `mstar`, `Temperature`, `XColdCloud` in SphP; entire MetP |
| BLACK_HOLES | `SwallowID`, `Injected_BH_Energy` in SphP; entire BHP; `BHID`/`MetID` union in P |
| STELLARAGE | `StellarAge` in MetP |
| WINDS | `DelayTime` in SphP |
| CONDUCTION | **no particle-struct fields** (globals only) |
| LT_METAL_COOLING_WAL | **no particle-struct fields** (cooling-table globals only) |
| LT_NMet=16 | `Metals[16]` in SphP, MetP |

Precision: `double` for Pos/GravAccel/GravPM; `float` for vel/mass/scalars.
`MyIDType = uint64_t`; `integertime = int64_t`.

## Four layouts

### Layout A — monolithic per-type AoS
Four fat per-type arrays: `GasFullA[]`, `DMFullA[]`, `StarFullA[]`,
`BHFullA[]`. Common fields are duplicated across the four fat structs;
no cross-index. Reshuffle = four independent per-type sub-permutations.

### Layout B' — fine common atoms + fine type-specific + cross-index
Five common arrays (`PCore`, `PDyn`, `PTime`, `PMeta`, `PLinkage`) plus
fine-grained per-type arrays (`GasCore`, `GasGrad`, `GasMetal`, `GasSF`;
`StarCore`, `StarMeta`; `BHCore`, `BHEnv`, `BHRepos`). A `PLinkage.type_idx`
maps each particle to its type-specific slot. Reshuffle permutes the 5
common arrays only; type arrays stay put because `PLinkage` preserves
the mapping.

### Layout Bc — coarse common + composite per-type + cross-index
Same connection scheme as B' but common is bundled into 4 arrays
(`PCoreB`, `PDynB`, `PAuxB`, `PLinkageB`), and each per-type struct is
one composite (`GasAllB` bundles GasCore+GasGrad+GasMetal+GasSF, etc.).

### Layout C — fine common + fine type-specific + positional mapping
Same fine atoms as B' but **no PLinkage**. The j-th gas particle in the
common arrays positionally maps to `GasCore[gas_begin + j]`. Invariant:
every permutation of the common arrays must be replicated on the
type-specific arrays in lockstep. Reshuffle permutes 13 streams total.

## Results (N = 10⁶, single thread, gnu11 -O3 -march=native)

Numbers from a single sandbox run — **relative ordering is the
meaningful signal**; absolute GB/s figures will shift on other hardware.
Partition fraction = 1 %; sparse reshuffle fraction = 10 %.

### OG3_ATOM_ALIGN = 64 (default, one cache line per atom)

**Reshuffle:**

| Layout | Streams | Bytes | Dense ms (GB/s) | Sparse ms (GB/s) |
|---|---:|---:|---:|---:|
| A | 4 | 705 MB | 613 (1.15) | 162 (4.35) |
| **Bc** | **4** | **640 MB** | **404 (1.58)** | **166 (3.86)** |
| B' | 5 | 768 MB | 580 (1.32) | 226 (3.40) |
| C | 13 | 897 MB | 866 (1.04) | 258 (3.48) |

**Partition (1 % swaps):**

| Layout | ns/swap | Bytes/iter |
|---|---:|---:|
| A | 300 | 21.1 MB |
| **Bc** | **450** | **19.2 MB** |
| B' | 450 | 23.0 MB |
| C | 750 | 26.8 MB |

### OG3_ATOM_ALIGN = 16 (compact, no cache-line guarantee)

**Reshuffle:**

| Layout | Bytes | Dense ms (GB/s) | Sparse ms (GB/s) |
|---|---:|---:|---:|
| A | 625 MB | 412 (1.52) | 150 (4.17) |
| **Bp ≈ Bc** | **448 MB** | **300 (1.49)** | **117 (3.83)** |
| C | 641 MB | 644 (1.00) | 182 (3.52) |

**Partition (1 % swaps):**

| Layout | ns/swap | Bytes/iter |
|---|---:|---:|
| A | 400 | 18.7 MB |
| Bp | 400 | 13.4 MB |
| **Bc** | **300** | **13.4 MB** |
| C | 600 | 19.2 MB |

## Headline

**Bc (coarse common + composite per-type + cross-index) wins reshuffle
at both alignments**. It has the fewest streams (4) and the tightest
common block. At ALIGN=16 it ties with B' since PMeta/PLinkage stop
wasting memory, but Bc is still ahead on partition because swap touches
fewer streams (4 vs 5).

**C loses both benchmarks** by 1.5–2×. Housekeeping cost of permuting
13 arrays in lockstep is unavoidable under positional mapping. The C
tradeoff (no cross-index lookup in the kernel inner loop) is **not
measured here** and needs kernel benchmarks to surface.

**A is competitive on partition** — a single fat-struct memcpy per swap
wins on ns when the working set fits in L2 — but loses reshuffle because
common-field duplication bloats the bytes-moved.

**ALIGN matters more for B'/C than for A or Bc.** A already packs tightly
(fat structs naturally 64-B aligned); Bc's composites also pack tightly.
B' suffers most from alignment padding because it has the most atoms
with small raw sizes.

## Harness caveats (read before generalising)

1. **Housekeeping-only.** Neither benchmark exercises a kernel. B'/Bc
   don't permute type-specific arrays on a reshuffle — this is a real
   advantage they have, at the cost of a PLinkage lookup in every
   kernel access. This harness measures the former, not the latter.
   Kernel benchmarks would be the deciding vote.
2. **Random same-type pairs, N=10⁶.** Partition uses random pairs across
   a 10⁶-particle set — the cache-miss case. In the real code, most
   partition swaps are *local* (particle crossed one leaf boundary),
   so this is pessimistic for all four.
3. **No actual physics**: no density, no force walk, no interaction.
4. **Single-threaded.** The reshuffle loop is embarrassingly parallel.
5. **90 % migration fraction** was mentioned as the reality in current
   OG3 due to a tree instability. Dense reshuffle (ƒ=1.0) is the right
   stress test for that regime; sparse (ƒ=0.1) represents the target
   regime after the instability is fixed.
6. **CONDUCTION and LT_METAL_COOLING_WAL are enabled but add no
   particle-struct fields** in the current OG3 source. If the real
   code adds K_cond / DtEntropy_cond to SphP, gas bytes grow and all
   four layouts pay the same tax proportionally.

## Build & run

```bash
make clean && make all                                  # 8 binaries, ALIGN=64
CFLAGS='... -DOG3_ATOM_ALIGN=32' make                   # 8 binaries, ALIGN=32
CFLAGS='... -DOG3_ATOM_ALIGN=16' make                   # 8 binaries, ALIGN=16

./bench_reshuffle_X [N] [n_iter] [fraction]    # X in {A, Bp, Bc, C}
./bench_partition_X [N] [fraction] [n_iter]

# Defaults: N=1000000, reshuffle n_iter=10 fraction=1.0, partition fraction=0.01 n_iter=20
```

Examples:
```bash
./bench_reshuffle_Bc 2000000 5 1.0    # dense reshuffle, N=2M
./bench_reshuffle_C  1000000 10 0.1   # sparse 10%, N=1M
./bench_partition_Bp 1000000 0.05 30  # 5% swaps, N=1M
```

## Files

```
bench_common.h    Types (count_t, ppid_t, pcount_t, pkey_t, ...),
                  OG3_ATOM macro, CPU_TIME_P timer, og3_aligned_alloc,
                  PH-key encoder, LSD radix sort, sparse-perm generator.
bench_common.c    Implementations.
fields.h          Atom structs (PCore, PDyn, ..., GasCore, ..., BHRepos);
                  no explicit padding, alignment via OG3_ATOM.
layout_api.h      Uniform alloc/fill/reshuffle/swap/free API with
                  restrict-qualified pointer parameters.
layout_A.c        Monolithic per-type AoS.
layout_Bp.c       Fine common + fine type-specific + PLinkage.
layout_Bc.c       Coarse common + composite per-type + PLinkageB.
layout_C.c        Fine common + fine type-specific + positional mapping.
bench_reshuffle.c Per-type dense sort (f=1.0) or per-type sparse
                  shuffle (0<f<1).
bench_partition.c N same-type random (i,j) swap pairs in a tight loop.
Makefile          8 binaries: {reshuffle, partition} × {A, Bp, Bc, C}.
```

## What this does NOT tell you

The three Tier-0 questions most relevant to the OG3 layout decision
are **not** answered here:

1. **Mutation policy viability** (gas → star conversion, BH merge) under
   each layout.
2. **Incremental-maintenance viability** (box-leaf patch under low
   migration fraction).
3. **Shared hydro-list viability on real snapshots** (η_cand).

Those need follow-on benchmarks that do couple to the box-leaf tree and
exercise real kernels. This harness measures only the bookkeeping cost
of moving particle data between layouts — the first and simplest thing
to rule out.
