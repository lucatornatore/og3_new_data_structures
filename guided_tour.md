# Guided Tour — OpenGadget3 C++ v2 Prototype, Bottom Up

*A walk through the v2 prototype from the most primitive pieces upward,
explaining for each what happens at compile time vs at run time.*

*Three layouts are covered. Sections 1–7 describe the shared foundation
(`opg_common/`); sections 8–10 cover layout B' (fine + cross-index);
sections 11–13 cover layout C (fine + positional); sections 14–17 cover
layout B (coarse + cross-index).*

*Each section is self-contained: you can read only the sections for the
layout you care about after reading the foundation. Cross-references are
explicit.*

---

## Part I — Shared foundation (`opg_common/`)

Every layout depends on these. None of them says anything about particles.

### 1. `scalar_types.hpp` — the vocabulary

**What it is.** Precision aliases, the `ParticleType` enum, a 3-vector
POD, the `ParticleComponent` concept, and the
`OPG_ASSERT_PARTICLE_COMPONENT` macro.

**What happens at compile time.** The aliases resolve to concrete types
based on the `OPG_PRECISION` build flag:

```cpp
#if OPG_PRECISION == 0                    // single
using real_t    = float;
using pos_t     = float;
using acc_t     = float;
#elif OPG_PRECISION == 2                  // double
using real_t    = double;
using pos_t     = double;
using acc_t     = double;
#else                                      // mixed (default)
using real_t    = float;
using pos_t     = double;                 // positions always double
using acc_t     = double;                 // accelerations always double
#endif
```

`ParticleType` is a plain enum class with `Gas, DM1, DM2, DM3, Star, BH`
mapped to values 0..5. `NTYPES = 6`.

The `ParticleComponent<T>` concept:

```cpp
template<typename T>
concept ParticleComponent =
    std::is_trivially_copyable_v<T> &&
    std::is_standard_layout_v<T>;
```

This is what every particle struct must satisfy. Trivially-copyable
because `memcpy` must work (MPI_BYTE, arena reshuffle). Standard-layout
because struct field offsets must be portable across translation units.

The convenience macro:

```cpp
#define OPG_ASSERT_PARTICLE_COMPONENT(...)                                  \
    static_assert(std::is_trivially_copyable_v<__VA_ARGS__>,                \
                  #__VA_ARGS__ " must be trivially copyable ...");          \
    static_assert(std::is_standard_layout_v<__VA_ARGS__>,                   \
                  #__VA_ARGS__ " must be standard-layout ...")
```

Applied after every particle struct. Failure gives a pointed diagnostic
("GasCore must be trivially copyable") rather than the uninformative
"ParticleComponent not satisfied".

**What happens at runtime.** Nothing — these are all type definitions.

---

### 2. `physics_config.hpp` — compile-time physics policy

**What it is.** `PhysicsConfig`, the preset configurations, the
physics concepts, and the `optional_ptr` / `optional_field` primitives
with their ergonomic helpers.

**What happens at compile time.**

#### 2.1 `PhysicsConfig` is an NTTP of class type

```cpp
struct PhysicsConfig {
    bool hydro_sph         = true;
    bool hydro_mfm         = false;
    bool magnetic          = false;
    bool stellar_evolution = true;
    int  n_metal_species   = OPG_NMET_DEFAULT;
    bool black_holes       = true;
    /* ... ~20 fields ... */

    consteval bool has_metals() const noexcept {
        return stellar_evolution || (star_formation && n_metal_species > 0);
    }
    consteval bool has_leapfrog_predictor() const noexcept {
        return hydro_sph || hydro_mfm;
    }
};
```

This is a C++20 "structural" type — aggregate, public members, no
user-declared copy/move/dtor. It can therefore be used as a non-type
template parameter. `template<PhysicsConfig Cfg> class X {};` is legal,
and inside `X`, `Cfg` is a compile-time value.

#### 2.2 The presets

Four `inline constexpr PhysicsConfig` values live in namespace
`opg::common::physics_configs`: `GravityOnly`, `StandardSPH`, `FullMHD`,
`MFMHydro`. Every field is set explicitly. Values that must be consistent
across presets (like `n_metal_species`) read from the build-time macro
`OPG_NMET_DEFAULT`.

#### 2.3 The physics concepts

```cpp
template<PhysicsConfig Cfg> concept HasSPH         = Cfg.hydro_sph;
template<PhysicsConfig Cfg> concept HasMagnetic    = Cfg.magnetic;
template<PhysicsConfig Cfg> concept HasMetals      = Cfg.has_metals();
template<PhysicsConfig Cfg> concept HasBHSpin      = Cfg.black_holes && Cfg.bh_spin;
template<PhysicsConfig Cfg> concept HasLeapfrog    = Cfg.has_leapfrog_predictor();
/* ... etc ... */
```

These are used two ways: as guards for `if constexpr` and as `requires`
clauses on member functions.

#### 2.4 `optional_ptr` and `optional_field`

```cpp
template<bool Enabled, typename T>
struct optional_ptr { T* p = nullptr; };

template<typename T>
struct optional_ptr<false, T> {};            // empty

template<bool Enabled, typename T>
struct optional_field { T data; };

template<typename T>
struct optional_field<false, T> {};           // empty
```

Combined with `[[no_unique_address]]`, the disabled specialisations
contribute zero bytes to their containing class.

- `optional_ptr<E, T>` holds a `T*`. Used by containers that own arrays
  of `T`.
- `optional_field<E, T>` holds a `T` value. Used inside particle structs
  to carry a conditionally-present sub-block.

#### 2.5 The ergonomic helpers

`opt_read(field, fallback)` and `opt_if(field, lambda)` hide the
`if constexpr` gate and the `.data` hop in the common cases:

```cpp
real_t B = opt_read(record.mag, bfield3_t{});

opt_if(record.mag, [&](auto& m){
    m.B     = new_B;
    m.div_B = 0.0f;
});
```

`opt_if` uses `if constexpr (Enabled) fn(f.data);`, so the lambda is
never instantiated for the disabled specialisation (which has no `.data`
member). `opt_read` returns the fallback when disabled.

**What happens at runtime.**

- `optional_ptr<true, T>` is a `T*`; `.p` is read and used like any other
  pointer. Runtime cost: one load.
- `optional_ptr<false, T>` has no data members; the compiler emits no
  load and no store.
- `opt_if` on the enabled specialisation inlines to the lambda body;
  on the disabled specialisation, inlines to a no-op.

---

### 3. `memory_arena.hpp` — single-allocation bump arena

**What it is.** A very simple arena allocator: one `posix_memalign` at
construction, `allocate<T>(n)` returns a pointer into the block and
advances a watermark. No individual `free`; the arena is deallocated
wholesale at destruction.

**Why.** Three reasons:
1. No per-allocation overhead, no fragmentation.
2. Memory layout is predictable: everything the container owns lives
   in one contiguous block.
3. Easy to swap for MPI-3 shared-window allocation later: just change
   the backing `posix_memalign` to `MPI_Win_allocate_shared`.

**What happens at compile time.** Template type deduction on
`allocate<T>(n)` produces one specialisation per `T`. The alignment
computation picks `max(alignof(T), config_.alignment)` at compile time.

**What happens at runtime.**

```cpp
MemoryArena arena(total_bytes, ArenaConfig{});        // one big allocation
PCore*   core  = arena.allocate<PCore>(N);            // bump pointer by N*sizeof(PCore)
GasCore* gas   = arena.allocate<GasCore>(N_gas);      // bump pointer further
```

Each `allocate<T>` rounds up to the alignment and returns a typed pointer.
A watermark (`save_point()`/`restore(p)`) lets us unwind for reuse.

---

### 4. `box_leaf_base.hpp` — the leaf descriptor base

**What it is.** `BoxLeafBase` (the data every layout's leaf needs) and
`SortHelper` (the 16-byte sort element for the PH-key sort).

```cpp
struct BoxLeafBase {
    idx_t common_begin;           // first common slot in this leaf
    idx_t common_end;             // one-past-last common slot
    idx_t type_offset[NTYPES];    // offset within [begin, end) of each type's run
    idx_t type_count [NTYPES];    // length of each type's run
    uint8 flags;
    int16 _reserved;
};

struct alignas(16) SortHelper {
    pkey_t  key;                  // 8 B — Peano-Hilbert key
    uint8   type;                 // 1 B
    uint8   _pad[3];
    idx_t   original_idx;         // 4 B — position before sort
};
static_assert(sizeof(SortHelper) == 16);
```

`SortHelper` is 16 bytes with the key first: radix-sort optimal.

**What happens at compile time.** Struct layout is fixed. `NTYPES`
propagates from the enum.

**What happens at runtime.**

- Each leaf stores its range in the common arrays (`[common_begin,
  common_end)`) and the per-type sub-ranges inside that range
  (`type_offset[]`, `type_count[]`).
- A kernel that wants to iterate over all gas in a leaf walks
  `[common_begin + type_offset[Gas], common_begin + type_offset[Gas] +
  type_count[Gas])`.
- Each layout extends `BoxLeafBase` with its own specifics (B' adds a
  `PLinkage` array pointer; C adds per-type array base indices).

---

### 5. `permutation.hpp` — the R8 mechanism

**What it is.** `apply_permutation`, `ArrayRegistry`, and `PHKeySorter`.
This is the file that solves the "reshuffle all arrays in lockstep"
problem.

#### 5.1 `apply_permutation`

```cpp
template<ParticleComponent T>
void apply_permutation(T* array, const idx_t* perm, count_t n, T* temp) noexcept {
    for (count_t i = 0; i < n; ++i) temp[i] = array[perm[i]];
    std::memcpy(array, temp, n * sizeof(T));
}
```

Constrained on `ParticleComponent`: we're going to `memcpy`, so the type
must be trivially copyable. Out-of-place gather into `temp`, then one
`memcpy` back. `perm[i] = j` means: after the call, `array[i]` holds
what was at `array[j]` before.

#### 5.2 `ArrayRegistry` — type erasure via stateless lambda

```cpp
template<std::size_t MaxArrays = 32>
class ArrayRegistry {
    struct Entry {
        void*  base;
        size_t elem_size;
        count_t elems;
        const char* name;
        void (*permute)(void*, const idx_t*, count_t, void*);
    };
    Entry entries_[MaxArrays]{};
    std::size_t count_ = 0;

public:
    template<ParticleComponent T>
    void register_array(T* base, count_t n, const char* name = "") noexcept {
        entries_[count_++] = Entry{
            base, sizeof(T), n, name,
            [](void* b, const idx_t* p, count_t n, void* tmp) {
                apply_permutation(static_cast<T*>(b), p, n, static_cast<T*>(tmp));
            }
        };
    }

    void reshuffle_all(const idx_t* perm, count_t n, void* scratch) const noexcept {
        for (std::size_t i = 0; i < count_; ++i) {
            entries_[i].permute(entries_[i].base, perm, n, scratch);
        }
    }
};
```

The key line is the lambda inside `register_array`. It captures nothing,
but it knows the concrete type `T` at registration time. The compiler
turns it into a plain function pointer that knows how to cast `void*` to
`T*` and call `apply_permutation<T>`.

**Compile time:** one specialisation of `register_array<T>` per particle
struct type. Each produces a function pointer baked for that `T`.

**Runtime:** `reshuffle_all` is a tight loop over `entries_`, calling each
function pointer. GCC can inline through the indirect call at `-O2` when
the call site is hot.

#### 5.3 `PHKeySorter`

Three static helpers: `build_helpers` (walk particles, build `SortHelper`
buffer), `sort_by_key` (the comparator), `extract_permutation` (produce
the `idx_t[]` that the registry consumes). Plus
`subsort_by_type_within_leaf` which is used after the PH-key sort to
group particles by type inside each leaf.

The `sort_by_key` comparator tiebreaks on `original_idx` so equal-key
ties resolve deterministically — this is how we guarantee restart files
are bitwise-reproducible across runs.

---

## Part II — Layout B' (fine + cross-index)

Namespace `opg::layout_b`. Each access-pattern cluster is its own small
struct with its own array. Cross-index via `PLinkage`.

### 6. `opg_layout_b/particles.hpp` — all the fine atoms

**What it is.** One file containing every particle struct B' needs.
Compile-time only (pure type definitions).

Common-side structs (all trivially copyable, all standard layout):

| Struct  | alignas | Content |
|---|---|---|
| `PCore`  | 64 | pos, mass, key, type, flags, time_bin, leaf_idx |
| `PDyn`   | 32 | vel, grav_acc, old_acc |
| `PLeap`  | 32 | dp, vel_pred (leapfrog predictor) |
| `PTime`  | 32 | ti_begin, ti_current, dt_step, grav_cost |
| `PMeta`  | 32 | id, num_ngb, true_ngb |
| `PPotential` | — | potential, pm_potential |

Gas-side structs (template parameters hidden from callers via container
type aliases):

| Struct  | alignas | Content |
|---|---|---|
| `GasCore`  | 32 | hsml, density, pressure, entropy, hydro_acc, sound_speed |
| `GasGrad`  | 32 | div/curl v, ∇v matrix, α, f_balsara |
| `GasMag`   | 32 | B, B_pred, div_B, phi, clean_vel |
| `GasMetal<N>` | 32 | metals[N], temperature, mass/egy res, z_smooth |
| `GasSF`    | 32 | sfr, delay_time, egy_sn, multiphase |
| `GasMFM<D>` | 32 | num_dens, alpha_slope[D+2], internal_energy |
| `GasCR<Np,Ne>` | 32 | cr proton/electron norm, slope, cut, pressure |

Plus `StarCore<N>`, `StarMeta`, `BHCore`, `BHEnv`, `BHRepos`, `BHSpin`,
`BHKinFB`.

All go through `OPG_ASSERT_PARTICLE_COMPONENT` — a layout that forgot one
would fail to compile.

---

### 7. `opg_layout_b/box_leaf.hpp` — `BoxLeaf` + `PLinkage`

`BoxLeaf` extends `BoxLeafBase` with nothing extra — under B', the leaf
needs only the common-range + per-type-offset bookkeeping that the base
already provides. The type-specific arrays are reached via `PLinkage`,
not via the leaf.

`PLinkage` is the cross-index: one `idx_t` per common slot. Meaning
depends on `PCore[j].type`:

```cpp
struct PLinkage { opg::common::idx_t type_idx; };

// Interpretation:
//   PCore[j].type == Gas   → PLinkage[j].type_idx indexes GasCore[]
//   PCore[j].type == Star  → PLinkage[j].type_idx indexes StarCore[]
//   PCore[j].type == BH    → PLinkage[j].type_idx indexes BHCore[]
//   PCore[j].type == DM*   → PLinkage[j].type_idx unused
```

---

### 8. `opg_layout_b/particle_container.hpp` — `ParticleContainer<Cfg>`

**What happens at compile time.**

The container is templated on `PhysicsConfig Cfg`. Inside the class body,
`Cfg` is a compile-time constant. Dispatch uses three mechanisms
combined:

```cpp
template<PhysicsConfig Cfg>
class ParticleContainer {
public:
    using GasMetalT = GasMetal<Cfg.n_metal_species>;        // pinned once
    using GasMFMT   = GasMFM<3>;
    using GasCRT    = GasCR<Cfg.cr_proton_bins, Cfg.cr_electron_bins>;
    using StarCoreT = StarCore<Cfg.n_metal_species>;

    // Always-present accessors
    PCore*   core()     noexcept { return core_;   }
    PDyn*    dyn()      noexcept { return dyn_;    }
    PTime*   time()     noexcept { return time_;   }
    PMeta*   meta()     noexcept { return meta_;   }
    PLinkage* linkage() noexcept { return linkage_;}

    // Conditional accessors — absent when physics disabled
    PLeap*     leap()      noexcept requires HasLeapfrog<Cfg>;
    GasMag*    gas_mag()   noexcept requires HasMagnetic<Cfg>;
    GasMetalT* gas_metal() noexcept requires HasMetals<Cfg>;
    /* ... etc ... */

private:
    PCore*    core_    = nullptr;
    /* ... always-allocated pointers ... */
    PLinkage* linkage_ = nullptr;

    [[no_unique_address]] optional_ptr<HasLeapfrog<Cfg>,  PLeap>      leap_storage_;
    [[no_unique_address]] optional_ptr<HasMagnetic<Cfg>,  GasMag>     gas_mag_storage_;
    [[no_unique_address]] optional_ptr<HasMetals<Cfg>,    GasMetalT>  gas_metal_storage_;
    /* ... etc ... */

    ArrayRegistry<> reg_common_;     // ONE registry for common arrays only
};
```

The compile-time cascade:

1. `ParticleContainer<StandardSPH>` instantiates with `Cfg.magnetic =
   false`.
2. `HasMagnetic<StandardSPH>` evaluates to `false`.
3. The `requires HasMagnetic<Cfg>` on `gas_mag()` fails → `gas_mag` is
   not a member of the class.
4. `optional_ptr<false, GasMag>` specialises to empty → `gas_mag_storage_`
   is 0 bytes.
5. The constructor's `if constexpr (HasMagnetic<Cfg>)` branch is not
   taken → `arena.allocate<GasMag>(...)` is never called.

A write to `container.gas_mag()` is a compile error with diagnostic
`'class opg::layout_b::ParticleContainer<...>' has no member named
'gas_mag'`. Not a null-pointer dereference.

**What happens at runtime.**

```cpp
ArenaConfig acfg{};
MemoryArena arena(64 * 1024 * 1024, acfg);
Capacity cap{ .n_max = 10000, .n_max_gas = 5000, ... };
ParticleContainer<StandardSPH> c(arena, cap);
```

Construction:

1. Arena bumps the watermark for each enabled array
   (`core_`, `dyn_`, `time_`, `meta_`, `linkage_`, `gas_core_`,
   `gas_grad_`, `gas_metal_`, ...).
2. `register_all_arrays()` calls `reg_common_.register_array<T>(ptr, n,
   name)` for each of the common arrays (6 for StandardSPH, 7 for
   FullMHD). Type-specific arrays are NOT registered — under B' they
   don't participate in PH-key reshuffle.
3. Runtime counts (`n_part_`, `n_gas_`, ...) start at 0. Caller fills
   the arrays from an IC file or from MPI migration.

Reshuffle:

```cpp
// Caller has: keys[], types[], scratch space.
PHKeySorter::build_helpers(keys, types, n, helpers);
PHKeySorter::sort_by_key(helpers, n);
PHKeySorter::extract_permutation(helpers, n, perm, nullptr);
c.common_registry().reshuffle_all(perm, n, scratch);
```

After the reshuffle, all 6 (or 7) common arrays are in PH-key order,
consistently permuted. `PLinkage[j]` values rode along with `PCore[j]`,
so the cross-indices still point to the right (un-permuted) type-
specific slots.

The type-specific arrays (`GasCore[]`, `StarCore[]`, etc.) are NOT
touched by the reshuffle. If memory locality for them matters, a
separate (layout-specific) `compact_type_arrays` pass can reorder them;
it runs much less often than the common reshuffle.

---

### 9. `opg_layout_b/reshuffle.hpp` and `mutation.hpp`

**`reshuffle.hpp`** defines `reshuffle_common(container, perm, scratch)`
as a thin wrapper over `container.common_registry().reshuffle_all(...)`.

**`mutation.hpp`** sketches `mutate_gas_to_star` (stub at v6 — the
pattern is the same as B-coarse but with more common arrays to swap;
implementable in ~60 lines when needed).

---

### 10. `examples/example_layout_b.cpp` — end-to-end runtime

The example:

1. Allocates an arena and a `ParticleContainer<StandardSPH>`.
2. Reports the registry size (should be 6 or 7 for SPH).
3. Seeds a small particle set with declining keys: `[100, 90, 80, 70, 60]`.
4. Calls `reshuffle_common` with permutation `[4, 3, 2, 1, 0]`.
5. Verifies keys afterwards are `[60, 70, 80, 90, 100]`.
6. Demonstrates that `container.gas_mag()` on a `GravityOnly` container
   is a compile error (commented out).

Run output:

```
=== B' prototype demonstration ===
Layout: B' (Common + clustered type-specific, explicit cross-index)

StandardSPH container:
  registry size   : 6  (common arrays; type arrays don't move under B')
GravityOnly container:
  registry size   : 5
  (gas_mag/leap/star_core are absent at compile time for GravityOnly)

Before reshuffle: keys = 100 90 80 70 60
Permutation     : 4 3 2 1 0
After reshuffle : keys = 60 70 80 90 100
```

---

## Part III — Layout C (fine + positional mapping)

Namespace `opg::layout_c`. Same atoms as B', no `PLinkage`. Connection
via box-leaf descriptor base indices. Four registries.

### 11. `opg_layout_c/particles.hpp` — the atoms (sans PLinkage)

Same file as B' minus the `PLinkage` struct. The type-specific arrays
are reached via leaf-descriptor indexing, not via cross-index. This is
the key qualitative difference from B'.

### 12. `opg_layout_c/box_leaf.hpp` — `BoxLeaf` with per-type bases

```cpp
struct BoxLeaf : public BoxLeafBase {
    idx_t gas_begin,  gas_end;     // range in the gas arrays
    idx_t star_begin, star_end;    // range in the star arrays
    idx_t bh_begin,   bh_end;      // range in the BH arrays
    idx_t staged_star_count;       // for the mutation protocol
    idx_t staged_bh_count;
    /* ... */
};
```

The positional invariant: the j-th gas particle in the leaf lives at
`common_begin + type_offset[Gas] + j` in the common arrays AND at
`gas_begin + j` in the gas arrays. Same j on both sides. No stored
cross-index; the descriptor *is* the mapping.

This invariant is powerful and fragile: any operation that permutes
arrays must permute all of them atomically in a way that keeps the
matching. Break it once and every subsequent kernel call silently
reads wrong data.

### 13. `opg_layout_c/particle_container.hpp` — FOUR registries

**Why four?** Because the four array-length classes have four different
lengths:

- Common arrays have `n_part` elements.
- Gas arrays have `n_gas` elements.
- Star arrays have `n_star` elements.
- BH arrays have `n_bh` elements.

A reshuffle that permutes `PCore[0..n_part)` has no meaning applied to
`GasCore[0..n_gas)`. The two arrays have different lengths. Each
length-class gets its own `ArrayRegistry`, and each registry is reshuffled
with its own permutation.

```cpp
template<PhysicsConfig Cfg>
class ParticleContainer {
    ArrayRegistry<> reg_common_;
    ArrayRegistry<> reg_gas_;
    ArrayRegistry<> reg_star_;
    ArrayRegistry<> reg_bh_;

    ArrayRegistry<>& common_registry() noexcept { return reg_common_; }
    ArrayRegistry<>& gas_registry()    noexcept { return reg_gas_;    }
    ArrayRegistry<>& star_registry()   noexcept { return reg_star_;   }
    ArrayRegistry<>& bh_registry()     noexcept { return reg_bh_;     }
};
```

### 13a. `opg_layout_c/reshuffle.hpp` — `PermutationBundle`

```cpp
struct PermutationBundle {
    idx_t* perm_common;   // length n_part
    idx_t* perm_gas;      // length n_gas
    idx_t* perm_star;     // length n_star
    idx_t* perm_bh;       // length n_bh
};

void reshuffle_all(ParticleContainer<Cfg>&, const PermutationBundle&, scratch...);
```

The four permutations are derived from ONE sorted `SortHelper[]` buffer
by `derive_permutations`. They must agree on the "correct relative
order of the particle types" so the positional invariant survives the
reshuffle.

### 13b. `opg_layout_c/mutation.hpp` — the staged protocol

Gas → star under C is the hardest case. The protocol:

1. Claim a star slot `s`; copy common data to `StarCore[N_star++]`
   staging.
2. Mark `PCore[j].type` with the staged sentinel bit (high bit).
3. **In-leaf swap**: move the flagged common slot to the END of the
   leaf's gas sub-range across ALL registered common AND gas arrays.
4. `leaf.type_count[Gas] -= 1`; `leaf.staged_star_count += 1`.
5. Pending-state particles are excluded from SPH loops by the descriptor
   range, not a per-particle branch.
6. At next rebuild, `integrate_staged_mutations` promotes staged entries
   into the proper star sub-range and renumbers.

Without step 3, the next kernel that iterates the leaf's gas sub-range
would walk off the positional invariant: `GasCore[gas_begin+j]` would
still contain the original gas data, but `PCore[common_begin+j]` would
now have `type == Star`. The in-leaf swap is the surgery that keeps the
invariant intact.

### Runtime cost

Mutation: ~10–17 array swaps per mutation under FullMHD. Expensive.
This is a real Tier-0 question: does the elegance of positional mapping
survive realistic mutation rates?

---

## Part IV — Layout B (coarse + cross-index)

Namespace `opg::layout_b_coarse`. Coarse granularity (7 arrays), same
connection scheme as B' (cross-index). Adds pack/unpack for MPI
exchange and compaction for orphan reclamation.

### 14. `opg_layout_b_coarse/particles.hpp` — fewer, wider structs

**Composition principle.** The coarse structs CONTAIN the fine atoms as
named members rather than redeclaring their fields. This means adding a
field to (for example) `GasCore` automatically propagates to `GasAllB`.

Common side:

```cpp
struct alignas(64) PCoreB {            // 64 B, not templated
    pos3_t  pos;     mass_t  mass;
    pkey_t  key;     uint8   type;
    uint8   flags;   int16   time_bin;
    idx_t   leaf_idx;
};

template<PhysicsConfig Cfg>
struct alignas(64) PDynB {             // 64 B with or without leapfrog
    vel3_t  vel;
    acc3_t  grav_acc;
    real_t  old_acc;
    real_t  _pad0;
    struct LeapState { vel3_t vel_pred; };
    [[no_unique_address]] optional_field<HasLeapfrog<Cfg>, LeapState> leap;
};

template<PhysicsConfig Cfg>
struct alignas(64) PAuxB {             // 64 B
    pid_t      id;
    time_int_t ti_begin, ti_current;
    int32      dt_step;
    real_t     grav_cost;
    real_t     num_ngb;
    idx_t      true_ngb;
    [[no_unique_address]] optional_field<HasPotentialOutput<Cfg>, real_t> potential;
};

struct alignas(8) PLinkageB {
    uint32_t type_idx;
    uint32_t _pad;
};
```

Type-specific side — composition from the fine atoms:

```cpp
template<PhysicsConfig Cfg>
struct alignas(64) GasAllB {
    GasCore core;
    [[no_unique_address]] optional_field<HasSPH<Cfg>,           GasGrad>                       grad;
    [[no_unique_address]] optional_field<HasMagnetic<Cfg>,      GasMag>                        mag;
    [[no_unique_address]] optional_field<HasMetals<Cfg>,        GasMetal<Cfg.n_metal_species>> metals;
    [[no_unique_address]] optional_field<HasStarFormation<Cfg>, GasSF>                         sfr;
};

template<PhysicsConfig Cfg>
struct alignas(64) StarAllB {
    StarCore<Cfg.n_metal_species> core;
    StarMeta                      meta;
};

template<PhysicsConfig Cfg>
struct alignas(64) BHAllB {
    BHCore  core;
    BHEnv   env;
    BHRepos repos;
    [[no_unique_address]] optional_field<HasBHSpin<Cfg>,      BHSpin>  spin;
    [[no_unique_address]] optional_field<HasBHKineticFB<Cfg>, BHKinFB> kinfb;
};
```

**Sizes under StandardSPH** (mixed precision):

| Struct | Size |
|---|---|
| `PCoreB`  | 64 B  (one cache line) |
| `PDynB`   | 64 B  (one cache line, with or without leapfrog) |
| `PAuxB`   | 64 B  (one cache line) |
| `PLinkageB` | 8 B |
| `GasAllB` | 256 B |
| `StarAllB` | 128 B |
| `BHAllB`  | 320 B |

---

### 15. `opg_layout_b_coarse/particle_container.hpp`

The container holds seven arrays for FullMHD:

- Four common arrays (`core`, `dyn`, `aux`, `linkage`).
- One `gas_all_` array (always allocated, may be size 0 under GravityOnly).
- Optional `star_storage_`, `bh_storage_` via `optional_ptr`.

**One common-side registry of 4 arrays.** That is the structural advantage
over B': reshuffle touches 4 arrays, not ~17.

Compile-time: same pattern as B' — `requires`-constrained accessors,
`optional_ptr` for conditionally-allocated storage, `if constexpr` in
the constructor.

Runtime:

```cpp
ParticleContainer<StandardSPH> c(arena, cap);
// reshuffle_common(c, perm, scratch)  →  4-array permutation
// type-specific arrays not moved; linkage values ride with common entries
```

---

### 16. `pack_unpack.hpp` — `ParticleRecord<Cfg>` and the distribute.c bridge

**Fixed-size POD for MPI exchange.**

```cpp
template<PhysicsConfig Cfg>
struct alignas(64) ParticleRecord {
    PCoreB       core;
    PDynB<Cfg>   dyn;
    PAuxB<Cfg>   aux;
    alignas(64) std::byte typed_payload[detail::typed_payload_bytes<Cfg>()];
};
```

`typed_payload_bytes<Cfg>()` is a `consteval` that returns the max of
`sizeof(GasAllB)`, `sizeof(StarAllB)`, `sizeof(BHAllB)`. Under
StandardSPH: 320 B (BHAllB is largest). Total record: 512 B.

**pack_record** copies PCoreB, PDynB, PAuxB from the container's arrays
into the record via `std::memcpy` (not struct assignment — see
§8 of implementation_notes.md for why). The typed payload is zeroed
first, then overwritten with the matching `GasAllB`/`StarAllB`/`BHAllB`
content based on `record.core.type`. DM particles leave the typed
payload as zero.

**unpack_record** is the inverse: copy common fields to the destination
container's next slot; claim a slot in the appropriate type-specific
array; rewrite linkage.

**Distribute.c bridge.** The existing C harness takes a `-DDATA_SIZE=N`
flag meaning "N ints per exchange cell". `print_record_size StandardSPH
ints` prints 128 (= 512 B / 4 B), which the C build then uses directly.

---

### 17. `mutation.hpp` and `compact.hpp`

**Mutation.** `mutate_gas_to_star`:

1. Claim star slot s; init `StarAllB[s]`.
2. Set `PCoreB[j].type = Star`; set `PLinkageB[j].type_idx = s`.
3. The old `GasAllB` slot is orphaned.
4. **In-leaf swap** across the 4 common arrays to keep the leaf's
   type sub-ordering correct (gas, DM, star, BH).

Four-array swap per mutation vs C's ~17-array swap. Cheaper by a
factor of ~4.

**Compaction.** `compact_type_arrays(container)`:

1. Walk all common entries. For each type, mark which slots in the
   type-specific array are still "live" (pointed to by some `PLinkageB`).
2. Left-pack the live slots into a compact array.
3. Rewrite `PLinkageB` values so the compacted indices remain correct.

Runs much less often than PH-reshuffle — every N rebuilds, or when the
orphan fraction exceeds a threshold. O(n_part) scan + O(n_gas + n_star +
n_bh) rewrite.

---

## Summary: a single picture of compile-time vs runtime

| What | When | How |
|---|---|---|
| `PhysicsConfig` presets | compile | `inline constexpr PhysicsConfig StandardSPH{...};` |
| Disabled struct → 0 bytes | compile | `optional_ptr<false, T> = {};` + `[[no_unique_address]]` |
| Disabled accessor → absent | compile | `requires HasX<Cfg>` |
| `gas_metal()` returns `GasMetal<11>*` | compile | `using GasMetalT = GasMetal<Cfg.n_metal_species>` |
| Type erasure for reshuffle | compile-to-function-ptr | stateless lambda captured at registration |
| Actual reshuffle | runtime | for-loop over `entries_`, one indirect call per array |
| MPI pack/unpack | runtime | `std::memcpy` of fixed-size record (Layout B) |
| Gas → star mutation | runtime | in-leaf swap across common arrays; orphan typed slot |
| PH-key sort determinism | runtime | tiebreak on `original_idx` in the comparator |

---

## How to explore the code in order

If you want to read the code from the bottom up in the order this tour
presents:

```
1. opg_common/include/opg_common/types/scalar_types.hpp
2. opg_common/include/opg_common/types/physics_config.hpp
3. opg_common/include/opg_common/memory/memory_arena.hpp
4. opg_common/include/opg_common/tree/box_leaf_base.hpp
5. opg_common/include/opg_common/permutation/permutation.hpp

Then, for B' (fine + cross-index):
6. opg_layout_b/include/opg_layout_b/particles.hpp
7. opg_layout_b/include/opg_layout_b/box_leaf.hpp
8. opg_layout_b/include/opg_layout_b/particle_container.hpp
9. opg_layout_b/include/opg_layout_b/reshuffle.hpp
10. examples/example_layout_b.cpp

Then, for C (fine + positional):
11. opg_layout_c/include/opg_layout_c/particles.hpp  (same atoms, no PLinkage)
12. opg_layout_c/include/opg_layout_c/box_leaf.hpp   (per-type base indices)
13. opg_layout_c/include/opg_layout_c/particle_container.hpp  (four registries)
14. opg_layout_c/include/opg_layout_c/reshuffle.hpp
15. opg_layout_c/include/opg_layout_c/mutation.hpp
16. examples/example_layout_c.cpp

Then, for B (coarse + cross-index):
17. opg_layout_b_coarse/include/opg_layout_b_coarse/particles.hpp
18. opg_layout_b_coarse/include/opg_layout_b_coarse/particle_container.hpp
19. opg_layout_b_coarse/include/opg_layout_b_coarse/reshuffle.hpp
20. opg_layout_b_coarse/include/opg_layout_b_coarse/pack_unpack.hpp
21. opg_layout_b_coarse/include/opg_layout_b_coarse/mutation.hpp
22. opg_layout_b_coarse/include/opg_layout_b_coarse/compact.hpp
23. examples/example_layout_b_coarse.cpp
24. examples/print_record_size.cpp
```

---

*End of guided tour.*
