# OG3 layout microbenchmarks, C draft v2

This directory contains a pure-C benchmark harness for four resident-layout
candidates for the OpenGadget3 rewrite:

| Layout | File | Meaning |
|---|---|---|
| A | `layout_A.c` | one fat per-type AoS array plus a logical global slot directory |
| B' | `layout_Bp.c` | fine common atoms + fine type-specific arrays + `PLinkage` |
| Bc | `layout_Bc.c` | coarse common atoms + composite type-specific arrays + `PLinkageB` |
| C | `layout_C.c` | fine common atoms + fine type-specific arrays + positional mapping, with resident `type_idx[]` kept for this benchmark |

The harness is deliberately C-only. It does not use the C++ `ArrayRegistry`,
physics-policy machinery, or optional-field templates. It hard-codes one
representative physics configuration and measures the cost of physically
maintaining layout invariants.

## What changed in v2

The reshuffle benchmark now uses a **single global Peano-Hilbert ordering**.
Particles are not grouped by type before sorting. Dense mode builds a global
permutation by stable radix-sorting `(PH_key, original_idx)`; stable LSD radix
sort preserves the input order for equal keys, so the original index is the
implicit deterministic tie-break.

Layout A cannot physically expose one global particle array, because its data
owners are the four per-type fat arrays. It now carries a logical `SlotA[]`
directory:

```c
typedef struct {
    uint8_t type;
    count_t type_idx;
} SlotA;
```

A global PH reshuffle reorders this directory and also reorders each per-type
fat array in the order in which that type appears in the global PH sequence.
The tree can therefore address global slot intervals, while kernels that need
payload fields pay one directory lookup.

Other implemented changes:

- `count_t` is used for particle counts, permutations, swap pairs, type indices,
  and every particle-array index. It is currently `uint64_t`.
- The Morton proxy was replaced with a real 3D Peano-Hilbert key using the
  Reinecke/Gadget lookup-table implementation in `hilbert.c` / `hilbert.h`.
- Sparse reshuffle now uses a real changed-destination copy path when fewer than
  25% of slots move. Dense permutations still use the streaming gather + copy
  path.
- DEBUG builds add type-payload markers and `layout_verify_deep()`. Non-DEBUG
  builds compile the deep verification to a stub.
- The public layout API now includes `layout_set_key()`, `layout_get_pos()`,
  and `layout_set_pos()`, so benchmarks and future tree code can update common
  particle fields without knowing whether they live in a common stream or in a
  layout-A fat per-type struct.
- Byte reporting now labels dense payload copy traffic as a 4-pass estimate:
  source read + scratch write + scratch read + destination write. A coarse
  permutation-index traffic estimate is reported separately.
- Bc scratch allocation was fixed: reshuffle scratch is sized only for the
  largest common stream, not for the large type-specific composites.
- `CPU_TIME_P` was replaced by `static inline double cpu_time_p(void)`.
- `assign_types()` handles `N < 2` without unsigned underflow.
- `bench_sizes_{A,Bp,Bc,C}` prints actual `sizeof` values and per-layout byte
  accounting for the current `OG3_ATOM_ALIGN`.

## Build

```bash
make clean && make all
```

The default flags are:

```bash
-std=gnu11 -O3 -march=native -Wall -Wextra -fno-strict-aliasing
```

The code also builds cleanly with:

```bash
CFLAGS='-std=gnu11 -O2 -Wall -Wextra -Wpedantic -fno-strict-aliasing' make all
```

For DEBUG marker verification:

```bash
CFLAGS='-std=gnu11 -O0 -g -Wall -Wextra -DDEBUG -fno-strict-aliasing' make all
```

## Binaries

```bash
./bench_reshuffle_X [N] [n_iter] [fraction]   # X in A Bp Bc C
./bench_partition_X [N] [fraction] [n_iter]
./bench_sizes_X [N]
```

Defaults:

- `bench_reshuffle`: `N=1000000`, `n_iter=10`, `fraction=1.0`
- `bench_partition`: `N=1000000`, `fraction=0.01`, `n_iter=20`
- `bench_sizes`: `N=1000000`

Examples:

```bash
./bench_reshuffle_Bc 1000000 10 1.0   # dense global PH sort reshuffle
./bench_reshuffle_C  1000000 10 0.1   # sparse global reshuffle
./bench_partition_Bp 1000000 0.05 30  # same-type swap primitive
./bench_sizes_C 1000000               # actual sizeof/accounting report
```

## Current benchmark scope

This is still a **single-MPI-task, single-threaded** benchmark. It does not run
MPI and is not OpenMP-parallel yet. The domain-exchange benchmark should be a
second step that reuses the region/intersection logic from `distribute.c` but
measures only local pack/unpack and self-copy costs before adding actual MPI.

The partition benchmark is still intentionally a primitive same-type swap test.
The next tree-specific benchmark should consume the real box-leaf builder and
measure:

1. global PH sort + full layout reshuffle;
2. leaf descriptor construction with `(start,end)` for all particles and at
   least one extra `(gas_start,gas_end)` range for hydro locality;
3. between-rebuild repartitioning by leaf ID without full key sorting.

## Notes on interpretation

Do not compare v1 and v2 reshuffle numbers directly. v1 sorted within type
ranges; v2 sorts all particles globally by PH key. That is the operation needed
for a PH-key tree whose leaves correspond to contiguous global particle slots.

For layout C, this benchmark keeps `type_idx[]` by design, per the current
experiment choice. Its memory footprint and refresh pass are therefore part of
C's measured bookkeeping cost.


## v5 partition update

The partition benchmark now uses position-driven particle generation and a
standalone box-leaf builder adapted from the uploaded `create_boxleaves`
prototype. `particle_init.[ch]` provides plain and clustered 21-bit coordinate
generation; `boxleaves.[ch]` builds PH-cube leaves from globally PH-sorted keys
using target occupancy, tolerance, and maximum side length. `bench_partition.c`
perturbs particle positions, recomputes keys, repartitions particles by their
new old-leaf membership, and keeps gas contiguous inside each leaf via the
updated `[gas_begin, gas_end)` descriptors.
