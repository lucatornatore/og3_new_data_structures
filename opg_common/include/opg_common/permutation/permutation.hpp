/**
 * @file permutation.hpp
 * @brief Permutation primitives, type-erased ArrayRegistry, PH-key sort helper.
 *
 * The authoritative low-level rules in this file are:
 *   - component moves are byte-wise (`memcpy`), not aggregate assignment;
 *   - scratch is untyped bytes, so alignment of the caller's scratch buffer
 *     does not leak into correctness;
 *   - element swaps are also byte-wise, so padding moves deterministically.
 */

#ifndef OPG_COMMON_PERMUTATION_PERMUTATION_HPP
#define OPG_COMMON_PERMUTATION_PERMUTATION_HPP

#include "../types/scalar_types.hpp"
#include "../tree/box_leaf_base.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>

namespace opg::common {

// =============================================================================
// apply_permutation -- out-of-place gather into byte scratch, memcpy back
// =============================================================================
//
// perm[i] = j means: after the call, array[i] holds what was in array[j].
//
// The scratch buffer is untyped bytes on purpose:
//   1. callers may pass arena scratch or std::byte storage without having to
//      prove alignment for T;
//   2. the gather is byte-wise, so padding bytes are preserved exactly.

template<ParticleComponent T>
void apply_permutation(T* array,
                       const idx_t* perm,
                       count_t n,
                       void* scratch) noexcept
{
    auto* temp = static_cast<std::byte*>(scratch);
    for (count_t i = 0; i < n; ++i) {
        std::memcpy(temp + static_cast<size_t>(i) * sizeof(T),
                    &array[perm[i]],
                    sizeof(T));
    }
    std::memcpy(array, temp, static_cast<size_t>(n) * sizeof(T));
}

template<ParticleComponent T>
void apply_permutation(T* array,
                       const idx_t* perm,
                       count_t n,
                       T* temp) noexcept
{
    apply_permutation(array, perm, n, static_cast<void*>(temp));
}

// =============================================================================
// swap_elements -- byte-wise element swap preserving padding
// =============================================================================

template<ParticleComponent T>
void swap_elements(T* array, idx_t a, idx_t b) noexcept
{
    if (a == b) return;
    alignas(T) std::byte tmp[sizeof(T)];
    std::memcpy(tmp,       &array[a], sizeof(T));
    std::memcpy(&array[a], &array[b], sizeof(T));
    std::memcpy(&array[b], tmp,       sizeof(T));
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
        assert(count_ < MaxArrays);
        if (count_ >= MaxArrays) std::abort();

        entries_[count_].base      = static_cast<void*>(base);
        entries_[count_].elem_size = sizeof(T);
        entries_[count_].elems     = elem_count;
        entries_[count_].name      = name;
        entries_[count_].permute   =
            [](void* b, const idx_t* p, count_t n, void* tmp) noexcept {
                apply_permutation(static_cast<T*>(b), p, n, tmp);
            };
        entries_[count_].swap_two  =
            [](void* b, idx_t a, idx_t c) noexcept {
                swap_elements(static_cast<T*>(b), a, c);
            };
        ++count_;
    }

    void reshuffle_all(const idx_t* perm, count_t n, void* scratch) const noexcept {
        if (n == 0) return;
        assert(scratch != nullptr);
        if (scratch == nullptr) std::abort();

        for (std::size_t i = 0; i < count_; ++i) {
            const auto& e = entries_[i];
            assert(e.elems >= n);
            if (e.elems < n) std::abort();
            e.permute(e.base, perm, n, scratch);
        }
    }

    void swap_all(idx_t a, idx_t b) const noexcept {
        if (a == b) return;
        for (std::size_t i = 0; i < count_; ++i) {
            entries_[i].swap_two(entries_[i].base, a, b);
        }
    }

    std::size_t size()  const noexcept { return count_; }
    bool        empty() const noexcept { return count_ == 0; }
    void        clear() noexcept       { count_ = 0; }

    size_t max_elem_size() const noexcept {
        size_t m = 0;
        for (std::size_t i = 0; i < count_; ++i) {
            if (entries_[i].elem_size > m) m = entries_[i].elem_size;
        }
        return m;
    }

    size_t max_elem_bytes() const noexcept {
        size_t m = 0;
        for (std::size_t i = 0; i < count_; ++i) {
            const size_t b = entries_[i].elem_size * entries_[i].elems;
            if (b > m) m = b;
        }
        return m;
    }

private:
    struct Entry {
        void*       base      = nullptr;
        size_t      elem_size = 0;
        count_t     elems     = 0;
        const char* name      = "";
        void (*permute)(void*, const idx_t*, count_t, void*) noexcept = nullptr;
        void (*swap_two)(void*, idx_t, idx_t) noexcept = nullptr;
    };

    Entry       entries_[MaxArrays]{};
    std::size_t count_ = 0;
};

// =============================================================================
// PHKeySorter -- serial reference sorter (GPU-ready interface)
// =============================================================================

class PHKeySorter {
public:
    static void build_helpers(const pkey_t* keys,
                              const uint8* types,
                              count_t n,
                              SortHelper* helpers) noexcept {
        for (count_t i = 0; i < n; ++i) {
            helpers[i].key          = keys[i];
            helpers[i].type         = types[i];
            helpers[i]._pad[0]      = 0;
            helpers[i]._pad[1]      = 0;
            helpers[i]._pad[2]      = 0;
            helpers[i].original_idx = static_cast<idx_t>(i);
        }
    }

    static void sort_by_key(SortHelper* helpers, count_t n) noexcept {
        std::sort(helpers, helpers + n,
            [](const SortHelper& a, const SortHelper& b) {
                if (a.key != b.key) return a.key < b.key;
                return a.original_idx < b.original_idx;
            });
    }

    static void extract_permutation(const SortHelper* helpers,
                                    count_t n,
                                    idx_t* perm,
                                    idx_t* inverse) noexcept {
        for (count_t i = 0; i < n; ++i) {
            perm[i] = helpers[i].original_idx;
            if (inverse) inverse[helpers[i].original_idx] = static_cast<idx_t>(i);
        }
    }

    static void subsort_by_type_within_leaf(SortHelper* helpers,
                                            BoxLeafBase& leaf) noexcept {
        const idx_t begin = leaf.common_begin;
        const idx_t end   = leaf.common_end;

        std::stable_sort(helpers + begin, helpers + end,
            [](const SortHelper& a, const SortHelper& b) {
                if (a.type != b.type) return a.type < b.type;
                if (a.key  != b.key)  return a.key  < b.key;
                return a.original_idx < b.original_idx;
            });

        for (int t = 0; t < NTYPES; ++t) {
            leaf.type_offset[t] = 0;
            leaf.type_count[t]  = 0;
        }

        idx_t current_offset = 0;
        int   current_type   = -1;
        for (idx_t i = begin; i < end; ++i) {
            const int t = helpers[i].type;
            if (t != current_type) {
                leaf.type_offset[t] = current_offset;
                current_type        = t;
            }
            leaf.type_count[t] += 1;
            ++current_offset;
        }
    }
};

} // namespace opg::common

#endif // OPG_COMMON_PERMUTATION_PERMUTATION_HPP
