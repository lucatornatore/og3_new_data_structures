# OpenGadget3 C++ Particle Data Structure — v2 Prototype

The v2 resident-layer prototype for the C++ rewrite of OpenGadget3. Three
**independent** layout implementations sharing a small common foundation,
all in header-only C++20.

## Structure

```
opg_v2/
├── CMakeLists.txt
├── opg_common/include/opg_common/           # shared foundation — NO particle structs
│   ├── types/scalar_types.hpp               # precision aliases, ParticleType enum,
│   │                                          OPG_ASSERT_PARTICLE_COMPONENT
│   ├── types/physics_config.hpp             # PhysicsConfig (NTTP, consteval), presets,
│   │                                          concepts (HasSPH, HasMagnetic, HasMetals, ...),
│   │                                          optional_ptr, optional_field, opt_if, opt_read
│   ├── memory/memory_arena.hpp              # single-alloc arena with ArenaFillPolicy
│   ├── permutation/permutation.hpp          # apply_permutation, ArrayRegistry, PHKeySorter
│   └── tree/box_leaf_base.hpp               # BoxLeafBase, SortHelper
│
├── opg_layout_b/include/opg_layout_b/       # B' — fine + cross-index   (namespace opg::layout_b)
│   ├── particles.hpp                        # PCore, PDyn, PLeap, PTime, PMeta,
│   │                                          GasCore, GasGrad, GasMag, GasMetal<N>, GasSF,
│   │                                          GasMFM<D>, GasCR<Np,Ne>, StarCore<N>, StarMeta,
│   │                                          BHCore, BHEnv, BHRepos, BHSpin, BHKinFB
│   ├── box_leaf.hpp                         # BoxLeaf + PLinkage (cross-index)
│   ├── particle_container.hpp               # 1 common registry, ~17 fine arrays
│   ├── reshuffle.hpp
│   └── mutation.hpp
│
├── opg_layout_c/include/opg_layout_c/       # C — fine + positional mapping (namespace opg::layout_c)
│   ├── particles.hpp                        # same atoms as B' (no PLinkage)
│   ├── box_leaf.hpp                         # BoxLeaf + per-type base indices + staged counts
│   ├── particle_container.hpp               # 4 registries (common, gas, star, bh)
│   ├── reshuffle.hpp                        # PermutationBundle, reshuffle_all
│   └── mutation.hpp                         # stage_gas_to_star with in-leaf swap
│
├── opg_layout_b_coarse/include/opg_layout_b_coarse/  # B coarse + cross-index (namespace opg::layout_b_coarse)
│   ├── particles.hpp                        # atoms + coarse composites (PCoreB, PDynB<Cfg>,
│   │                                          PAuxB<Cfg>, PLinkageB, GasAllB<Cfg>,
│   │                                          StarAllB<Cfg>, BHAllB<Cfg>)
│   ├── particle_container.hpp               # 1 common registry of 4 arrays
│   ├── reshuffle.hpp
│   ├── pack_unpack.hpp                      # ParticleRecord<Cfg> for distribute.c bridge
│   ├── mutation.hpp                         # mutate_gas_to_star (4-array in-leaf swap)
│   └── compact.hpp                          # compact_type_arrays — reclaim orphans
│
└── examples/
    ├── example_layout_b.cpp
    ├── example_layout_c.cpp
    ├── example_layout_b_coarse.cpp
    └── print_record_size.cpp                # DATA_SIZE helper for distribute.c
```

## The three layouts

| Layout | Namespace | Granularity | Connection | Array count (FullMHD) |
|---|---|---|---|---|
| **B'** | `opg::layout_b` | fine (1 struct per access pattern) | explicit `PLinkage` cross-index | ~17 |
| **C**  | `opg::layout_c` | fine (same atoms as B') | positional mapping via leaf descriptor | ~17 across 4 registries |
| **B coarse** | `opg::layout_b_coarse` | coarse (few wider structs) | explicit `PLinkageB` cross-index | 7 |

Benchmarking B vs. B' isolates the granularity axis at fixed connection
scheme. Benchmarking B' vs. C isolates the connection axis at fixed
granularity. See design_note_v5.docx for the decision framework.

## Dependency policy

- `opg_common/` is shared. It contains **no particle structs** — only generic
  plumbing (scalar types, physics config, memory, permutation primitives).
- Each layout OWNS its particle types inside its own namespace. There is no
  type sharing across layouts: `opg::layout_b::PCore`, `opg::layout_c::PCore`,
  `opg::layout_b_coarse::PCoreB` are three distinct types that the compiler
  will never silently confuse.
- This means identical-looking fields are duplicated across the three
  `particles.hpp` files. That is deliberate — the cost of the duplication
  (a few hundred lines) is lower than the cost of the confusion that came
  from trying to share them.

## Build and run

### With CMake

```bash
cmake -S . -B build
cmake --build build -j
./build/example_layout_b
./build/example_layout_c
./build/example_layout_b_coarse
```

### Without CMake

```bash
for LAY in b c b_coarse; do
  g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic \
      -I opg_common/include -I opg_layout_${LAY}/include \
      examples/example_layout_${LAY}.cpp -o example_layout_${LAY}
done
```

### Build-time configuration

- `-DOPG_NMET=<N>` — metal species count. Default 11. Propagates to all
  presets automatically (mirrors the `-DNMET=` flag of the C codebase).
- `-DOPG_CR_PROTON_BINS=<N>`, `-DOPG_CR_ELECTRON_BINS=<N>` — cosmic-ray bins.
- `-DOPG_PRECISION=0|1|2` — single / mixed (default) / double.

Requires a C++20 compiler. Tested with `g++ 13.3`. GCC ≥ 10 and Clang ≥ 12
should also work (class-type NTTPs and `consteval` predicates are needed).

## Verified

On `g++ 13.3 -std=c++20 -Wall -Wextra -Wpedantic -O2`:

- All three layouts compile clean.
- All three examples pass their round-trip / reshuffle / mutation checks.
- `example_layout_b_coarse` round-trips 37 mixed-type particles through
  `pack_record`/`unpack_record` with bit-for-bit equality.
- `-DOPG_NMET=7` override works across all three layouts.
- `B-coarse` output is bitwise-identical across repeated runs
  (determinism preserved for restart reproducibility).

## Design principles

1. **Independence over DRY.** Three prototypes, three namespaces, one
   shared foundation. Opening a layout's directory shows you everything
   that layout is: its data model, its container, its reshuffle, its
   mutation. No hops to `opg_common/` to understand what a particle looks
   like.

2. **Compile-time physics, zero runtime overhead.** `PhysicsConfig` is a
   C++20 NTTP; predicates are `consteval`; disabled features are absent
   from the class interface (not null pointers). Misuse is a compile
   error, not a crash.

3. **Authoritative reshuffle.** Each layout's container registers every
   component array at construction. Reshuffling traverses the registry;
   there is no way to add an array without including it in reshuffle.

4. **Determinism by construction.** The PH-key sort tiebreaks on
   `original_idx`, so equal keys produce a reproducible order across
   runs. Restart files are bitwise-reproducible.

5. **Padding matters for MPI.** `pack_record`/`unpack_record` use
   explicit `memcpy`, not struct assignment — so `alignas` trailing
   padding flows end-to-end.

## TODOs (Phase 2)

- GPU radix-sort swap-in for `PHKeySorter`. Interface is GPU-ready
  (`SortHelper` is 16 B with key-first layout).
- `integrate_staged_mutations` (C) — needs tree-rebuild infrastructure.
- Transient AoSoA tile layer for GPU kernels.
- MPI wiring (shared-window intra-node; node-master inter-node).
- NUMA-aware allocation — hooks exist in `ArenaConfig`.
- BH merger mutation (same pattern as gas→star).
- Grouped-by-type pack/unpack as a bandwidth-optimal alternative to the
  fixed-size `ParticleRecord<Cfg>` (Tier 2 — current version matches
  `distribute.c`'s data model directly).
