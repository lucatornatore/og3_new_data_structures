# opg_layout_b_coarse — Layout B (coarse variant)

A **coarse-grained** counterpart to `opg_layout_b` (B'). Same connection
scheme — explicit cross-index via `PLinkageB` — but with large type-specific
structs instead of fine-grained per-module arrays.

## Why this variant exists

The v4 design note treated granularity and connection scheme as bundled:
- B' was "Common + clustered hot/warm/cold macro-blocks, cross-index"
- C was "Common + per-type AoSS, positional mapping"

The v2 prototype implemented both B' and C at the same **fine** granularity
(~17 arrays each), isolating the *connection scheme* (cross-index vs positional).
That left granularity untested.

This layer adds the third candidate:

| Variant | Granularity | Connection | Array count |
|---|---|---|---|
| v2 B' (`opg_layout_b`)            | fine    | cross-index | ~17 |
| v2 C  (`opg_layout_c`)            | fine    | positional  | ~17 × 4 registries |
| v2 B coarse (`opg_layout_b_coarse`) | coarse  | cross-index | ~7 |

Running a benchmark pipeline against all three disentangles the
granularity axis (Q1) from the connection axis (Q2) from the maintenance
discipline axis (Q3). See §6 of design_note_v5.docx.

## Struct layout under B (coarse)

**Common side — 3-way split**

| Struct | Size (mixed prec.) | Role |
|---|---|---|
| `PCoreB`     | 64 B | pos, mass, key, type, flags, time_bin, leaf_idx. Read by everything. |
| `PDynB<Cfg>` | 64 B | vel, grav_acc, old_acc + optional leapfrog predictor. |
| `PAuxB<Cfg>` | 64 B | id, timestep fields, neighbour counts, cost + optional potential. |
| `PLinkageB`  |  8 B | `type_idx` — slot in {Gas,Star,BH}AllB. Local-only (not serialised). |

**Type-specific side — one struct per type**

| Struct | Size (StandardSPH) | Contents |
|---|---|---|
| `GasAllB<Cfg>`  | 256 B | `GasCore` + optional `GasGrad`, `GasMag`, `GasMetal<N>`, `GasSF` via `optional_field`. |
| `StarAllB<Cfg>` | 128 B | `StarCore<N>` + `StarMeta`. |
| `BHAllB<Cfg>`   | 320 B | `BHCore` + `BHEnv` + `BHRepos` + optional `BHSpin`, `BHKinFB`. |

**The 2-way common split is a trivial variant**: merge `PDynB` and `PAuxB`
into a single `PRestB` struct. One registry call disappears. That variant
is not in this repo but the edit is one file.

## Under the bonnet

- **PH-reshuffle** permutes only the four common arrays (Core, Dyn, Aux,
  Linkage). Type-specific arrays are not moved — `PLinkageB` values ride
  with the common slot they describe, so cross-indices remain valid.
  Registry size = 4, vs C's 4 registries × up to 17 arrays.

- **Mutation** (gas → star) is cheaper than under C because only the four
  common arrays participate in the in-leaf sub-sort repair. Orphans
  accumulate in `GasAllB` until the next `compact_type_arrays()` call.

- **Pack/unpack** uses a fixed-size `ParticleRecord<Cfg>` as the MPI
  transfer unit. Bytes flow through explicit `memcpy` to survive alignas
  padding round-trip. No special treatment of DM particles — they carry an
  unused typed-payload region. (The grouped-by-type alternative is less
  bandwidth-wasteful; it is left for a follow-up.)

## Files

```
opg_layout_b_coarse/include/opg_layout_b_coarse/
├── particles.hpp              PCoreB, PDynB<Cfg>, PAuxB<Cfg>, PLinkageB,
│                              GasAllB<Cfg>, StarAllB<Cfg>, BHAllB<Cfg>
├── particle_container.hpp     ParticleContainer<Cfg> with one common registry
├── reshuffle.hpp              reshuffle_common, do_reshuffle_by_key
├── pack_unpack.hpp            ParticleRecord<Cfg>, pack_record, unpack_record,
│                              pack_range, pack_scatter, unpack_range
├── mutation.hpp               mutate_gas_to_star (4-array in-leaf swap)
└── compact.hpp                compact_type_arrays — reclaim orphans
```

Example: `examples/example_layout_b_coarse.cpp` — allocates two containers,
seeds 37 mixed-type particles, runs PH-key reshuffle, round-trips everything
through `pack_record`/`unpack_record`, and verifies bit-for-bit equality.

## Bridge to `distribute.c`

`distribute.c` exchanges fixed-size `data_type` records:

```c
typedef struct { int data[DATA_SIZE]; } data_type;
```

where `DATA_SIZE` is set at compile time with `-DDATA_SIZE=N`. To run the
same harness with layout-B particles as the unit of exchange:

1. Build the `print_record_size` helper:
   ```bash
   cmake --build .  # builds the helper alongside the examples
   ```

2. Query the size for your chosen physics preset:
   ```bash
   ./print_record_size StandardSPH ints
   # prints e.g. 128  — that is DATA_SIZE for distribute.c
   ```

3. Build distribute.c with the matching payload:
   ```bash
   mpicc -O2 -DDATA_SIZE=128 distribute.c -o distribute_b_sph
   ```

4. Run with the exchange function you want to benchmark:
   ```bash
   mpirun -n 8 ./distribute_b_sph -g3 64 -i 2 2 2 -e 4 2 1 -x 0
   #            grid 64^3, initial 2×2×2 tasks, final 4×2×1, hypercube send
   ```

   The payload content is arbitrary in this benchmark; what you are
   measuring is the **cost of moving a fixed-size record** of the shape a
   layout-B particle would have. To compare against B' or C, repeat with
   `DATA_SIZE = sizeof(the largest fine-grained sub-struct)` — e.g.
   `DATA_SIZE = sizeof(PCore) / sizeof(int) = 16` — and multiply the
   per-record send count by the number of arrays in that layout. That is
   the fair "array-count scaling" measurement for the three variants.

5. For the "pack from real particle arrays" test (Tier 2 — not yet in this
   repo), the harness would replace distribute.c's `data_in` array with a
   `ParticleContainer` instance, call `pack_range()` to materialise a
   `ParticleRecord[]` send buffer, and feed that into the existing
   `exchange_hypercube` code path. This requires minor surgery on
   `distribute.c` to route its allocations through the container; see
   §8 of design_note_v5.docx for the workflow.

## Verified

On `g++ 13.3 -std=c++20 -Wall -Wextra -Wpedantic -O2`:

- Clean compile of all three layouts (B', C, B-coarse) and both examples
  plus `print_record_size`.
- Example round-trip: 37 mixed-type particles pack→unpack with bit-for-bit
  equality across all four common arrays and every type-specific record.
- Reshuffle correctness: PH-key sort produces non-decreasing keys, and
  permutation is applied consistently across all four common arrays.
- Mutation correctness: `mutate_gas_to_star` on a 3-gas-2-DM-1-star leaf
  produces the expected type sub-sort (gas,gas,DM,DM,star,star) and correct
  linkage.
- Compaction correctness: orphan removal from `compact_type_arrays` leaves
  live gas densities at the right slots after a mutation.

## Not yet implemented

- Grouped-by-type pack/unpack (Tier 2). The fixed-size `ParticleRecord`
  approach matches `distribute.c` cleanly; grouped-by-type is the
  bandwidth-optimal alternative and is a follow-up.
- BH merger mutation — same pattern as gas→star; straightforward.
- GPU radix-sort swap-in for `PHKeySorter`. Interface is GPU-ready.
- NUMA-aware allocation — `ArenaConfig` has the hook; actual path is Phase 2.
