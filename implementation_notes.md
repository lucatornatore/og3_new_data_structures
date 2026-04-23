# Implementation Notes — OpenGadget3 C++ v2 Prototype

*Companion to design_note_v6.docx. This document covers the C++ implementation
details: techniques used, choices made, and the engineering rationale that
shaped the resident-layer prototype.*

---

## 0. Scope and format

This is an engineering document, not a design document. The design note
(v6.0) answers *what* the architecture is and *why*. This document answers
*how* the C++ realises it, in enough detail that a collaborator picking up
the code can read it without re-deriving the choices.

Structure:

1. Foundation — what lives in `opg_common/` and why
2. The three layouts — what they share, what they don't
3. The C++20 feature set used, and why each one
4. Type erasure via stateless lambdas (the `ArrayRegistry` trick)
5. Zero-cost optional features — the five-layer mechanism
6. Ergonomic access to optional fields
7. Sort determinism and comparator pitfalls
8. MPI pack/unpack and the padding trap
9. Mutation: the three layouts' cost profiles
10. Build-time configuration
11. Notes on compiler behaviour
12. What's deliberately not here

---

## 1. Foundation: `opg_common/`

`opg_common/` is the shared plumbing. It is a *header-only* library with
five files and no particle structs. Its job is to provide the vocabulary
that every layout speaks and the primitives that every layout uses.

```
opg_common/include/opg_common/
├── types/scalar_types.hpp    # precision aliases, ParticleType enum,
│                               OPG_ASSERT_PARTICLE_COMPONENT, ParticleComponent concept
├── types/physics_config.hpp  # PhysicsConfig NTTP, presets, physics concepts,
│                               optional_ptr/field, opt_if/opt_read helpers
├── memory/memory_arena.hpp   # single-allocation bump arena
├── permutation/permutation.hpp  # apply_permutation, ArrayRegistry, PHKeySorter
└── tree/box_leaf_base.hpp    # BoxLeafBase, SortHelper, leaf_flags
```

Everything in `opg_common/` lives in `namespace opg::common`. Nothing
particle-specific lives here. Nothing layout-specific lives here. A layout
depends on `opg_common/` but not the reverse.

### Why header-only?

Three reasons:
- The templates (`PhysicsConfig Cfg`-parameterised classes) cannot be
  pre-compiled without knowing which presets will be used; separating
  headers and .cpp files would force a Phase-2 choice here.
- The code is small enough that compile time is not a bottleneck.
- It keeps the dependency graph trivial: `target_link_libraries(x INTERFACE
  opg_common)` is all that's needed.

---

## 2. The three layouts

Each layout lives in `opg_layout_*/` and owns its own particle types in its
own namespace:

| Directory | Namespace | Granularity | Connection |
|---|---|---|---|
| `opg_layout_b/`        | `opg::layout_b`        | fine (~17 atoms) | explicit `PLinkage` cross-index |
| `opg_layout_c/`        | `opg::layout_c`        | fine (same atoms, no `PLinkage`) | positional mapping via leaf descriptor |
| `opg_layout_b_coarse/` | `opg::layout_b_coarse` | coarse (7 arrays) | explicit `PLinkageB` cross-index |

Three distinct C++ types result from this: `opg::layout_b::PCore`,
`opg::layout_c::PCore`, `opg::layout_b_coarse::PCoreB`. The compiler will
**never** silently convert between them. This is by design. It's the
answer to the question "how do I prevent a helper from accidentally working
across layouts in ways that make benchmark contamination possible?"

### Cost: duplicated structs

`PCore`, `PDyn`, `GasCore`, etc. appear in both `opg_layout_b/particles.hpp`
and `opg_layout_c/particles.hpp` with identical field sets. That is the
cost of independence. If you find a bug in `BHCore`, you fix it in two or
three places.

The alternative — sharing the atoms via a common header — recreates the
confusion that led to the restructure (see v5 design note for the history).
For a benchmark prototype whose job is to produce three fully comparable
candidates, the duplication cost is smaller than the confusion cost.

---

## 3. C++20 feature set used

Five C++20 features do the heavy lifting. Each is used because nothing in
C++17 could replace it without runtime overhead.

### 3.1 Non-type template parameter (NTTP) of class type

```cpp
struct PhysicsConfig {
    bool hydro_sph = true;
    bool magnetic  = false;
    int  n_metal_species = OPG_NMET_DEFAULT;
    /* ... ~20 fields ... */
    consteval bool has_metals() const noexcept { ... }
};

template<PhysicsConfig Cfg>
class ParticleContainer { /* ... */ };

ParticleContainer<physics_configs::StandardSPH> c1;
ParticleContainer<physics_configs::FullMHD>     c2;   // distinct type
```

Before C++20, non-type template parameters could be integers, enums, or
pointers — not `struct`s. C++20 lifted this for "structural" types:
aggregates with public members, no user-declared copy/move/dtor, every
member itself structural. `PhysicsConfig` is carefully built to meet this.

The consequence is that `Cfg` is a compile-time value inside the container
body. `if constexpr (Cfg.magnetic)` works. `requires Cfg.magnetic` works.
Two containers with different physics configs have different types, so
mismatched kernels refuse to link.

### 3.2 `consteval` functions

```cpp
consteval bool has_metals() const noexcept {
    return stellar_evolution || (star_formation && n_metal_species > 0);
}
```

`consteval` is strictly stronger than `constexpr`. A `constexpr` function
*may* execute at compile time; a `consteval` function *must*. Calling a
`consteval` function on a non-constant-evaluated instance is a compile
error, not silent runtime evaluation.

This is used for predicates on `PhysicsConfig`. We never want
`cfg.has_metals()` to become a runtime branch at a hot path. `consteval`
guarantees the compiler will complain if someone writes such code.

### 3.3 Concepts

```cpp
template<PhysicsConfig Cfg> concept HasMagnetic = Cfg.magnetic;
template<PhysicsConfig Cfg> concept HasSPH      = Cfg.hydro_sph;
template<PhysicsConfig Cfg> concept HasMetals   = Cfg.has_metals();
```

These are the "feature is enabled" predicates. They're used both as
boolean guards for `if constexpr` and as `requires` clauses on member
functions.

### 3.4 `requires`-constrained member functions

```cpp
template<PhysicsConfig Cfg>
class ParticleContainer {
public:
    GasCore* gas_core() noexcept;                            // always present
    GasMag*  gas_mag()  noexcept requires HasMagnetic<Cfg>;  // absent when disabled
    GasMetalT* gas_metal() noexcept requires HasMetals<Cfg>;
    BHSpin*  bh_spin()   noexcept requires HasBHSpin<Cfg>;
};
```

When the constraint is not satisfied, the member function is **genuinely
absent** from the class interface. Calling `container.gas_mag()` on a
`StandardSPH` container (magnetic=false) is not a runtime null dereference;
it's a compile error with the diagnostic `no member named 'gas_mag'`.

This is the mechanism that turns "disabled physics is impossible to use
incorrectly" from a convention into a type-system guarantee.

### 3.5 `[[no_unique_address]]` + empty optional specialisations

```cpp
template<bool Enabled, typename T>
struct optional_ptr { T* p = nullptr; };

template<typename T>
struct optional_ptr<false, T> {};                            // empty

class ParticleContainer {
    [[no_unique_address]]
    optional_ptr<HasMagnetic<Cfg>, GasMag> gas_mag_storage_;  // 0 bytes when disabled
};
```

`[[no_unique_address]]` (C++20) allows an empty subobject to share its
address with an adjacent member. Combined with the empty specialisation,
disabled optional storage contributes zero bytes to the containing class.

The same trick is used for value-holding members via `optional_field`.

---

## 4. Type erasure via stateless lambdas (`ArrayRegistry`)

This is the single bit of runtime polymorphism in the whole codebase.
Everywhere else, dispatch is compile-time.

### The problem

The PH-key reshuffle must apply the same permutation to every component
array. The set of arrays depends on `PhysicsConfig` (compile-time) but
varies per container instantiation (so it isn't a fixed variadic pack).
We need a way to say: "the container knows which arrays exist; the
reshuffle should visit each one and call `apply_permutation<T>` with the
correct `T`".

### The solution

```cpp
struct Entry {
    void*       base;
    size_t      elem_size;
    count_t     elems;
    const char* name;
    void (*permute)(void*, const idx_t*, count_t, void*);
};

template<ParticleComponent T>
void register_array(T* base, count_t n, const char* name) noexcept {
    entries_[count_++] = Entry{
        base, sizeof(T), n, name,
        [](void* b, const idx_t* p, count_t n, void* tmp) {
            apply_permutation(static_cast<T*>(b), p, n, static_cast<T*>(tmp));
        }
    };
}
```

The lambda `[](void* b, ...){ apply_permutation(static_cast<T*>(b), ...); }`
captures nothing. It's therefore convertible to a plain function pointer
— no vtable, no heap allocation, no `std::function` overhead.

`T` is captured at **registration time** (when we know it) and baked into
the function pointer. At call time we have `void*` and a function pointer
that knows how to interpret it. Type erasure without polymorphism cost.

### Cost

One indirect call per registered array per reshuffle. With 5–20 arrays,
the reshuffle loop is a tight `for` with a cache-hot array of function
pointers. GCC inlines across the indirect call at `-O2`.

### Why not variadic templates?

Because the pack would need to vary per `PhysicsConfig`:
`std::tuple<PCore*, PDyn*, optional<GasMag*>, optional<GasMetal<N>*>, ...>`.
We'd spend all our compile budget on std::tuple machinery for a problem
that doesn't need it.

---

## 5. Zero-cost optional features — five layers

To make "disabled physics = invisible" work end to end, five layers
cooperate:

| Layer | Mechanism | Effect when disabled |
|---|---|---|
| 1. Config | `PhysicsConfig::magnetic = false` | bool in NTTP |
| 2. Concept | `HasMagnetic<Cfg> = Cfg.magnetic` | evaluates to `false` |
| 3. Accessor | `gas_mag() requires HasMagnetic<Cfg>` | member function absent |
| 4. Storage | `optional_ptr<HasMagnetic<Cfg>, GasMag>` + `[[no_unique_address]]` | 0 bytes |
| 5. Allocation | `if constexpr (HasMagnetic<Cfg>) arena.allocate<GasMag>(...)` | never called |

If any one layer is missing, the feature leaks cost somewhere. All five
are necessary for a truly zero-overhead disabled feature.

---

## 6. Ergonomic access to optional fields

For an `optional_field<Enabled, T>` nested inside a particle struct, the
access pattern is inherently conditional:

```cpp
if constexpr (HasMagnetic<Cfg>) {
    record.mag.data.B = new_B;     // .data.B hop required
}
```

The `if constexpr` gate is mandatory when the feature may be disabled —
without it, `record.mag.data.B` would fail to compile because the
disabled specialisation has no `.data` member. The `.data.` hop is
language plumbing, not intent.

Two helpers in `physics_config.hpp` hide the gate in the common cases:

### `opt_read(field, fallback)` — read with default

```cpp
real_t B = opt_read(record.mag, bfield3_t{});
// == record.mag.data when enabled
// == bfield3_t{} (zero-initialised) when disabled
```

### `opt_if(field, lambda)` — write / mutate conditionally

```cpp
opt_if(record.mag, [&](auto& m){
    m.B     = new_B;
    m.div_B = 0.0f;
});
```

The lambda is only instantiated (and only called) when `Enabled` is true.
Under the disabled specialisation this expands to a no-op; the lambda is
never compiled against the missing `.data` member.

### Why not a macro?

I deliberately did **not** wrap these in a macro. Macros encourage use;
the gate is still genuine documentation of "this field may not exist in
this build". The lambda version keeps the gate visible as a word
(`opt_if`) rather than hidden, while eliminating the verbose `.data`
hop and the type-existence check.

### Limitation

If your code needs to hold a *reference* to the enabled type
(`auto& m = record.mag.data;` followed by many mutations), the outer
`if constexpr` gate cannot be eliminated without macros or reflection.
That's genuinely how the language works; we accept it.

---

## 7. Sort determinism and the comparator trap

The PH-key sort is used every rebuild. Two concerns combined into one fix.

### Concern 1: strict weak ordering, not three-way compare

`std::sort`'s `Compare` parameter must be a **strict weak ordering** —
`comp(a, b)` returns `true` iff `a` must come before `b`. It is NOT a
three-way comparator like `qsort`'s `int compar(a, b)` returning -1/0/+1.

The mistake to avoid:

```cpp
std::sort(helpers, helpers + n, [](auto& a, auto& b) {
    return (a.key > b.key) - (a.key < b.key);   // WRONG
});
```

This returns -1, 0, or +1. Both -1 and +1 cast to `true`, so
`comp(a, b) && comp(b, a)` is true for distinct `a, b`. That violates
irreflexivity and asymmetry, and triggers **undefined behaviour** inside
`std::sort` — in practice an infinite loop or silently corrupted sort.

The correct boolean form:

```cpp
std::sort(helpers, helpers + n, [](auto& a, auto& b) {
    return a.key < b.key;
});
```

### Concern 2: determinism across runs

`a.key < b.key` is a valid strict weak ordering, but it is **not a total
order** when keys collide. `std::sort` can reorder collided keys
arbitrarily — and "arbitrarily" means the result depends on the library
implementation. Two runs of the same binary over the same input can
produce different permutations, which breaks restart reproducibility.

### Fix: tiebreak on `original_idx`

```cpp
std::sort(helpers, helpers + n, [](const SortHelper& a, const SortHelper& b) {
    if (a.key != b.key) return a.key < b.key;
    return a.original_idx < b.original_idx;
});
```

At 21 bits/dim, PH-key collisions are rare, so the tiebreak branch almost
never fires. Cost: a single predictable branch. Restart files are now
bitwise-reproducible across runs.

### Bonus: why not a `cmov`?

At `-O3 -march=native`, GCC emits a branch for this comparator, not a
`cmov`. That's correct: the tiebreak branch is highly predictable
(collisions are rare), and a predicted-not-taken branch is faster than a
`cmov` because the latter introduces a data dependency and forces an
unnecessary load on the common path. Assembly inspection confirms
`if/return`, ternary, and `std::tie` source forms all produce identical
16-instruction output at `-O2/-O3/-O3 -march=native`. The compiler
normalises the source form and chooses the branch based on predictability
heuristics.

---

## 8. MPI pack/unpack and the padding trap

`opg_layout_b_coarse/pack_unpack.hpp` uses a fixed-size POD
`ParticleRecord<Cfg>` as the MPI exchange unit. During initial testing,
the round-trip failed on `PAuxB` — `src.id == dst.id` yet `memcmp`
reported a mismatch.

### Diagnosis

C++ trivially-copyable aggregate assignment (`dst_struct = src_struct`) is
**not** obliged by the standard to copy trailing `alignas` padding. GCC
was performing member-wise assignment that left `dst`'s arena-originated
padding untouched. For in-process use this is invisible; for `MPI_BYTE`
transfer, the wire sees every byte the sender emits, including padding.
Two hosts with the same `alignas` layout but different arena initialisation
disagree on padding bytes.

### Fix

Use `std::memcpy` everywhere in `pack_record` / `unpack_record`:

```cpp
void pack_record(const ParticleContainer<Cfg>& container, idx_t j,
                 ParticleRecord<Cfg>& record) noexcept {
    std::memcpy(&record.core, &container.core()[j], sizeof(PCoreB));
    std::memcpy(&record.dyn,  &container.dyn()[j],  sizeof(record.dyn));
    std::memcpy(&record.aux,  &container.aux()[j],  sizeof(record.aux));
    /* ... typed payload ... */
}
```

### Generalisation

Any place POD structs move around with bitwise-reproducibility requirements
(MPI, checksums, restart files, distributed hash tables) must use
`memcpy`, not `=`. Struct assignment is fine for in-process copies; for
anything that crosses a trust boundary or a wire, be explicit.

---

## 9. Mutation: three layouts, three cost profiles

Gas → star mutation is the canonical hard case for any split layout. Each
layout has a different cost structure.

### Layout C (fine, positional)

Under C, the positional alignment is a global invariant: the j-th gas
particle in a leaf lives at `CommonCore[common_begin+j]` and
`GasCore[gas_begin+j]`. Mutating one gas particle can corrupt that
relationship for every surviving gas particle in the same leaf, unless we
perform an in-leaf swap that preserves it.

Protocol:
1. Claim star slot s; copy common data to `StarCore[N_star++]` staging.
2. Mark `CommonCore[j].type` with the staged sentinel bit.
3. **In-leaf swap**: move the flagged slot to the END of the leaf's gas
   sub-range across ALL ~10–17 common AND gas arrays.
4. Decrement `leaf.gas_count`; increment `leaf.staged_star_count`.
5. Full structural consistency restored at next rebuild (integrate-staged).

Cost: ~10–17 array swaps per mutation + one typed copy.

### Layout B' (fine, cross-index)

Under B', `PLinkage[j].type_idx` connects common to type-specific. Mutating
gas → star requires updating `PCore[j].type` and rewriting `PLinkage[j]`
to point at the new star slot. The old gas slot is orphaned until
compaction.

Cost: ~6 common-array swaps per mutation (to repair the leaf's type
sub-sort) + one typed copy + one linkage update.

### Layout B (coarse, cross-index)

Under B (coarse), the common arrays are 4 (Core, Dyn, Aux, Linkage). The
type-specific arrays (GasAllB, StarAllB, BHAllB) are not permuted by the
in-leaf swap — they are reached via `PLinkageB`.

Cost: **4 common-array swaps** per mutation + one typed copy + linkage
update + orphan left behind for compaction.

### Reclamation

B' and B (coarse) accumulate orphaned type-specific slots after
mutations. `compact_type_arrays()` is the periodic rewrite that reclaims
them: O(n_part) scan + O(n_gas + n_star + n_bh) rewrite.

---

## 10. Build-time configuration

### `OPG_NMET` — metal species count

Mirrors the C code's `-DNMET=...` convention:

```cpp
#ifndef OPG_NMET
#define OPG_NMET 11
#endif

inline constexpr int OPG_NMET_DEFAULT = OPG_NMET;
```

All `PhysicsConfig` presets read from `OPG_NMET_DEFAULT`. No preset
carries an independent 11 (the presets I wrote first did, which was a
mistake caught in code review).

### `OPG_CR_PROTON_BINS`, `OPG_CR_ELECTRON_BINS`

Same pattern for cosmic-ray bin counts.

### `OPG_PRECISION`

0 = single, 1 = mixed (default), 2 = double. Controls which float type
`real_t`, `pos_t`, `vel_t`, `acc_t`, `mass_t` map to in
`scalar_types.hpp`.

### Honest caveat

`OPG_NMET` is binary-wide. If you ever want two containers with different
metal counts in one binary, a preset must override explicitly:
`PhysicsConfig{.n_metal_species = 7, ...}`. The C code has the same
limitation, so this is not a regression.

---

## 11. Compiler behaviour notes

These are things that came up while the prototype was being written. None
are language-required; all are `g++ 13.3` at `-O2`/`-O3` on `x86_64`.

### 11.1 Source form doesn't matter for simple comparators

Tested on the PH-key comparator: `if/return`, ternary, and `std::tie`
forms produce identical 16-instruction assembly at `-O2`, `-O3`, and
`-O3 -march=native -mtune=native`. The compiler normalises the source
to an internal SSA form before code generation.

### 11.2 `cmov` is not always better than a branch

GCC chose a branch for the PH-key tiebreak, not a `cmov`. Correct: the
branch is highly predictable (collisions are rare), and `cmov` creates
a data dependency that can't be broken by speculation. Predictable
branches beat `cmov` in practice; unpredictable branches lose to `cmov`.

### 11.3 `alignas` padding under aggregate assignment

As covered in §8: don't assume `dst = src` copies padding for trivially-
copyable types. Use `memcpy` when bitwise reproducibility matters.

### 11.4 Stateless lambda → function pointer

Guaranteed by the standard (since C++11 for regular captureless lambdas,
explicitly tightened in C++17). Used aggressively in `ArrayRegistry`.

### 11.5 Warnings we care about

`-Wall -Wextra -Wpedantic` is the standard. `-Wpedantic` specifically
catches non-ASCII characters inside comments, which bit us once when em-
dashes were used in comment blocks. ASCII hyphens only, in comments that
will be read on any locale.

---

## 12. What's deliberately not here

The v2 prototype is resident-layer only. These are out of scope:

- **Transient AoSoA tile layer.** Phase 2. The prototype builds the input
  to this layer; the layer itself is where kernel-level GPU optimisation
  lives, and it's a different design problem (warp width, register
  pressure, coalescing) from the resident layer.
- **Tree walk, SPH/MFM kernels, gravity.** Phase 1/2. The prototype ensures
  the data model can feed them efficiently; the algorithms themselves are
  a separate workstream.
- **GPU radix sort.** The `PHKeySorter` interface is GPU-ready (`SortHelper`
  is 16 B with key-first for radix). Current implementation is `std::sort`.
- **MPI wiring.** Shared-window intra-node + node-master inter-node. The
  `ArenaConfig` has `use_mpi_shared` hooks; the actual code is Phase 2.
- **NUMA-aware allocation.** Same — hooks present, implementation deferred.
- **`integrate_staged_mutations` (C) body.** Needs tree-rebuild
  infrastructure.
- **Full SFR physics, BH merger mutation, stellar feedback recipes.**
  Content, not architecture. Will be filled in when the kernels come
  online.

---

*End of implementation notes.*
