#ifndef OPG_COMMON_PERMUTATION_PERMUTATION_HPP
#define OPG_COMMON_PERMUTATION_PERMUTATION_HPP

#include "../types/scalar_types.hpp"
#include "../tree/box_leaf_base.hpp"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iostream>

namespace opg::common {

// FIX: Byte-wise copy to prevent alignment UB
template<ParticleComponent T>
void apply_permutation(T* array, const idx_t* perm, count_t n, void* scratch) noexcept {
    std::byte* tmp_bytes = static_cast<std::byte*>(scratch);
    for (count_t i = 0; i < n; ++i) {
        std::memcpy(tmp_bytes + i * sizeof(T), &array[perm[i]], sizeof(T));
    }
    std::memcpy(array, tmp_bytes, static_cast<size_t>(n) * sizeof(T));
}

template<std::size_t MaxArrays = 32>
class ArrayRegistry {
public:
    ArrayRegistry() = default;

    template<ParticleComponent T>
    void register_array(T* base, count_t elem_count, const char* name = "") noexcept {
        if (count_ >= MaxArrays) {
            std::cerr << "CRITICAL: ArrayRegistry overflow for " << name << std::endl;
            std::abort();
        }
        entries_[count_].base      = static_cast<void*>(base);
        entries_[count_].elem_size = sizeof(T);
        entries_[count_].elems     = elem_count;
        entries_[count_].name      = name;
        
        entries_[count_].permute   = [](void* b, const idx_t* p, count_t n, void* tmp) {
            apply_permutation(static_cast<T*>(b), p, n, tmp);
        };

        // FIX 7: swap_fn for automated mutation handling
        entries_[count_].swap_fn   = [](void* b, count_t i, count_t j) {
            T* arr = static_cast<T*>(b);
            std::swap(arr[i], arr[j]);
        };
        
        ++count_;
    }

    void reshuffle_all(const idx_t* perm, count_t n, void* scratch) const noexcept {
        for (std::size_t i = 0; i < count_; ++i) {
            const auto& e = entries_[i];
            if (e.elems < n) {
                std::cerr << "CRITICAL: Reshuffle bounds violation for " << e.name << std::endl;
                std::abort();
            }
            e.permute(e.base, perm, n, scratch);
        }
    }

    // FIX 7: Generic swap_all dispatch
    void swap_all(count_t i, count_t j) const noexcept {
        for (std::size_t k = 0; k < count_; ++k) {
            entries_[k].swap_fn(entries_[k].base, i, j);
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
            size_t b = entries_[i].elem_size * entries_[i].elems;
            if (b > m) m = b;
        }
        return m;
    }

private:
    struct Entry {
        void* base;
        size_t      elem_size;
        count_t     elems;
        const char* name;
        void (*permute)(void*, const idx_t*, count_t, void*);
        void (*swap_fn)(void*, count_t, count_t);
    };

    Entry        entries_[MaxArrays]{};
    std::size_t  count_ = 0;
};

class PHKeySorter {
public:
    static void build_helpers(const pkey_t* keys, const uint8* types, count_t n, SortHelper* helpers) noexcept {
        for (count_t i = 0; i < n; ++i) {
            helpers[i].key          = keys[i];
            helpers[i].type         = types[i];
            helpers[i].original_idx = static_cast<idx_t>(i);
        }
    }

    static void sort_by_key(SortHelper* helpers, count_t n) noexcept {
        std::sort(helpers, helpers + n, [](const SortHelper& a, const SortHelper& b) {
            if (a.key != b.key) return a.key < b.key;
            return a.original_idx < b.original_idx;
        });
    }

    static void extract_permutation(const SortHelper* helpers, count_t n, idx_t* perm, idx_t* inverse) noexcept {
        for (count_t i = 0; i < n; ++i) {
            perm[i] = helpers[i].original_idx;
            if (inverse) inverse[helpers[i].original_idx] = static_cast<idx_t>(i);
        }
    }
};

} // namespace opg::common

#endif // OPG_COMMON_PERMUTATION_PERMUTATION_HPP
