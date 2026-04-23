/**
 * @file permutation.hpp
 * @brief Permutation primitives, type-erased ArrayRegistry, PH-key sort helper.
 *
 * ============================================================================
 * WHAT THIS FILE IS FOR
 * ============================================================================
 *
 * When the simulation advances, particles drift in space, so the memory
 * ordering that was "sorted by PH key" goes stale. Periodically we need to
 * REORDER all the particle arrays back to PH-key order. That is one problem.
 *
 * The data layout is split across many arrays: under layout C we have one
 * array per small struct (PCore[], PDyn[], PLeap[], ..., GasCore[], GasGrad[],
 * ..., and so on -- roughly a dozen under StandardSPH, up to ~17 under
 * FullMHD). The reorder must touch every one of them, applying the SAME
 * permutation, or the arrays fall out of sync and we corrupt per-particle
 * data. That is the hard part of the problem.
 *
 * This file provides three pieces that together solve both:
 *
 *   1. apply_permutation<T>(array, perm, n, scratch)
 *      the primitive. Given a permutation `perm` of length n, rewrite
 *      array[i] := old_array[perm[i]]. Out-of-place gather into `scratch`,
 *      then a single memcpy back into place. Constrained on
 *      ParticleComponent: the `memcpy` is only safe for trivially-copyable,
 *      standard-layout types.
 *
 *   2. ArrayRegistry<MaxArrays>
 *      a type-erased registry of component arrays. The container's
 *      register_all_arrays() calls registry.register_array<T>(...) once
 *      for every component array; each call records a plain function
 *      pointer (captured from a stateless lambda) that knows how to run
 *      apply_permutation on that specific T. At reshuffle time,
 *      registry.reshuffle_all(perm, n, scratch) fans out over every
 *      registered array, each with its own type-aware callback.
 *
 *      This is the only runtime polymorphism in the otherwise compile-time-
 *      dispatched architecture. The registry exists because the set of
 *      enabled arrays depends on PhysicsConfig, which is compile-time but
 *      varies per container instantiation. The cost is one indirect call
 *      per array per reshuffle, inlined by the compiler at the use site.
 *
 *   3. PHKeySorter
 *      the PH-key sort reference implementation. Walks the particles,
 *      builds a SortHelper buffer of (key, type, original_idx) triples,
 *      sorts by key, and produces the permutation that `ArrayRegistry`
 *      consumes. Currently a CPU std::sort; the interface is GPU-ready
 *      (SortHelper is 16 B with the key in the first 8 B, optimal for
 *      radix sort).
 *
 * ============================================================================
 * TYPICAL WORKFLOW FROM A CONTAINER'S POINT OF VIEW
 * ============================================================================
 *
 *   // At container construction:
 *   register_all_arrays();   // one line per component array, in one function
 *
 *   // At reshuffle time (every N timesteps, configurable):
 *   PHKeySorter::build_helpers(keys, types, n, helpers);
 *   PHKeySorter::sort_by_key(helpers, n);
 *   PHKeySorter::extract_permutation(helpers, n, perm, nullptr);
 *   registry.reshuffle_all(perm, n, scratch);
 *   // now every registered array is in PH-key order, all in lockstep.
 *
 * ============================================================================
 * WHY THE REGISTRY RATHER THAN A VARIADIC TEMPLATE PACK
 * ============================================================================
 *
 * The alternative "expand a fixed parameter pack of arrays at compile time"
 * approach would work if the set of arrays were statically known, but it
 * depends on PhysicsConfig (compile-time but per-instantiation), so the
 * pack would have to be assembled differently for each container type.
 * The registry sidesteps that by collecting type information at
 * registration time (which IS per-instantiation) and erasing it into a
 * plain function pointer. Zero vtable, zero heap allocation, one indirect
 * call per array per reshuffle, inlinable.
 */

#ifndef OPG_COMMON_PERMUTATION_PERMUTATION_HPP
#define OPG_COMMON_PERMUTATION_PERMUTATION_HPP

#include "../types/scalar_types.hpp"
#include "../tree/box_leaf_base.hpp"
#include <algorithm>
#include <cstddef>
#include <cstring>

namespace opg::common {

// =============================================================================
// apply_permutation — out-of-place gather into temp, memcpy back
// =============================================================================
//
// perm[i] = j means: after the call, array[i] holds what was in array[j].

template<ParticleComponent T>
void apply_permutation(T* array, const idx_t* perm, count_t n, T* temp) noexcept {
    for (count_t i = 0; i < n; ++i) {
        temp[i] = array[perm[i]];
    }
    std::memcpy(array, temp, static_cast<size_t>(n) * sizeof(T));
}

// =============================================================================
// ArrayRegistry<MaxArrays>
// =============================================================================

template<std::size_t MaxArrays = 32>
class ArrayRegistry {
public:
    ArrayRegistry() = default;

    template<ParticleComponent T>
    void register_array(T* base, count_t elem_count, const char* name = "") noexcept {
        if (count_ >= MaxArrays) return;     // silent no-op for now
        entries_[count_].base      = static_cast<void*>(base);
        entries_[count_].elem_size = sizeof(T);
        entries_[count_].elems     = elem_count;
        entries_[count_].name      = name;
        entries_[count_].permute   =
            [](void* b, const idx_t* p, count_t n, void* tmp) {
                apply_permutation(static_cast<T*>(b), p, n, static_cast<T*>(tmp));
            };
        ++count_;
    }

    /**
     * Apply the same permutation to every registered array.
     *
     * The permutation has length `n` and acts on the FIRST `n` entries of
     * each registered array; entries past position `n` are untouched. This
     * is the normal case for active-count < capacity workloads.
     *
     * Arrays whose registered capacity is < n are a caller bug (permutation
     * would run out of bounds); they are skipped silently here. In debug
     * builds you may wish to add an assert.
     *
     * `scratch` must be sized ≥ max_elem_size() * n bytes.
     */
    void reshuffle_all(const idx_t* perm, count_t n, void* scratch) const noexcept {
        for (std::size_t i = 0; i < count_; ++i) {
            const auto& e = entries_[i];
            if (e.elems < n) continue;   // caller bug: n > capacity
            e.permute(e.base, perm, n, scratch);
        }
    }

    std::size_t size()  const noexcept { return count_; }
    bool        empty() const noexcept { return count_ == 0; }
    void        clear() noexcept       { count_ = 0; }

    /** Largest per-element size across registered arrays, in bytes.
     *  Scratch for reshuffle_all(perm, n, ...) should be ≥ max_elem_size() * n. */
    size_t max_elem_size() const noexcept {
        size_t m = 0;
        for (std::size_t i = 0; i < count_; ++i) {
            if (entries_[i].elem_size > m) m = entries_[i].elem_size;
        }
        return m;
    }

    /** Bytes of the largest registered FULL array. Useful if you intend to
     *  reshuffle at full capacity. For partial reshuffles, use
     *  max_elem_size() * n instead. */
    size_t max_elem_bytes() const noexcept {
        size_t m = 0;
        for (std::size_t i = 0; i < count_; ++i) {
            size_t b = entries_[i].elem_size * entries_[i].elems;
            if (b > m) m = b;
        }
        return m;
    }

private:
    struct Entry {
        void*       base;
        size_t      elem_size;
        count_t     elems;
        const char* name;
        void (*permute)(void*, const idx_t*, count_t, void*);
    };

    Entry        entries_[MaxArrays]{};
    std::size_t  count_ = 0;
};

// =============================================================================
// PHKeySorter — serial reference sorter (GPU-ready interface)
// =============================================================================

class PHKeySorter {
public:
    // Populate helpers[i] = (keys[i], types[i], i)
    static void build_helpers(const pkey_t* keys, const uint8* types,
                              count_t n, SortHelper* helpers) noexcept {
        for (count_t i = 0; i < n; ++i) {
            helpers[i].key          = keys[i];
            helpers[i].type         = types[i];
            helpers[i]._pad[0]      = 0;
            helpers[i]._pad[1]      = 0;
            helpers[i]._pad[2]      = 0;
            helpers[i].original_idx = static_cast<idx_t>(i);
        }
    }

    // =========================================================================
    // sort_by_key — primary PH-key sort, deterministic
    // =========================================================================
    //
    // Note on the comparator: std::sort (and std::stable_sort) expect a
    // STRICT WEAK ORDERING predicate — a boolean `comp(a, b)` that returns
    // true iff `a` must come before `b`. It is NOT the three-way `strcmp`
    // convention (-1 / 0 / +1). So
    //
    //      [](a, b){ return a.key < b.key; }                    // correct
    //      [](a, b){ return (a.key > b.key) - (a.key < b.key);} // WRONG
    //
    // The three-way form returns -1, 0, +1 cast to bool — so -1 and +1
    // both evaluate to `true`, which violates irreflexivity and asymmetry
    // and leads to undefined behaviour inside std::sort (in practice: an
    // infinite loop or a corrupted sort).
    //
    // Determinism. The primary sort is stable by construction only on the
    // (key, original_idx) pair. A plain `a.key < b.key` comparator is NOT
    // stable on equal keys: std::sort may reorder them arbitrarily, which
    // means the resulting permutation depends on the library version. For
    // the PH-key reshuffle that matters: two runs of the same binary with
    // the same input should produce the same order, to make restart files
    // bit-reproducible. We therefore tiebreak on original_idx; the extra
    // cost is negligible because collisions are rare at 21 bits/dim.
    static void sort_by_key(SortHelper* helpers, count_t n) noexcept {
        std::sort(helpers, helpers + n,
            [](const SortHelper& a, const SortHelper& b) {
                if (a.key != b.key) return a.key < b.key;
                return a.original_idx < b.original_idx;
            });
    }

    // Extract the common-array forward and inverse permutations.
    static void extract_permutation(const SortHelper* helpers, count_t n,
                                    idx_t* perm, idx_t* inverse) noexcept {
        for (count_t i = 0; i < n; ++i) {
            perm[i] = helpers[i].original_idx;
            if (inverse) inverse[helpers[i].original_idx] = static_cast<idx_t>(i);
        }
    }

    // Stable sub-sort by type within a single leaf, then fill the leaf's
    // type_offset[] and type_count[] arrays.
    static void subsort_by_type_within_leaf(SortHelper* helpers,
                                            BoxLeafBase& leaf) noexcept {
        const idx_t begin = leaf.common_begin;
        const idx_t end   = leaf.common_end;

        std::stable_sort(helpers + begin, helpers + end,
            [](const SortHelper& a, const SortHelper& b) {
                if (a.type != b.type) return a.type < b.type;
                return a.key < b.key;
            });

        for (int t = 0; t < NTYPES; ++t) {
            leaf.type_offset[t] = 0;
            leaf.type_count [t] = 0;
        }
        idx_t current_offset = 0;
        int   current_type   = -1;
        for (idx_t i = begin; i < end; ++i) {
            const int t = helpers[i].type;
            if (t != current_type) {
                leaf.type_offset[t] = current_offset;
                current_type        = t;
            }
            leaf.type_count[t]++;
            ++current_offset;
        }
    }
};

} // namespace opg::common

#endif // OPG_COMMON_PERMUTATION_PERMUTATION_HPP
