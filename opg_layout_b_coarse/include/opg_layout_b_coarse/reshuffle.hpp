/**
 * @file reshuffle.hpp  (layout B, coarse)
 * @brief PH-key reshuffle driver for layout B.
 *
 * Under B the PH-key reshuffle permutes ONLY the common arrays (Core, Dyn,
 * Aux, Linkage). Type-specific arrays (GasAllB, StarAllB, BHAllB) are NOT
 * permuted here — they are reached through PLinkageB[j].type_idx, and the
 * values inside PLinkageB ride with the common entries under the
 * permutation, so the cross-indices remain valid.
 *
 * This is B's structural advantage over C: a single permutation of four
 * arrays. The registry remains authoritative for group membership and swap_all,
 * but the hot reshuffle loop itself can be a compile-time direct walk.
 */

#ifndef OPG_LAYOUT_B_COARSE_RESHUFFLE_HPP
#define OPG_LAYOUT_B_COARSE_RESHUFFLE_HPP

#include <algorithm>
#include <cstring>

#include "opg_common/permutation/permutation.hpp"
#include "opg_common/tree/box_leaf_base.hpp"
#include "particle_container.hpp"

namespace opg::layout_b_coarse {

using namespace opg::common;

// =============================================================================
// reshuffle_common — permute all common arrays by `perm`
// =============================================================================
// `perm` must be a length-`n` permutation of [0, n): the element at output
// position i takes the value at input position perm[i].
//
// `scratch` must be at least `max_elem_size() * n` bytes (= 64 B × n for B's
// largest common struct). The caller may supply an arena-allocated scratch
// pointer; the registry uses it as a temporary during each per-array
// apply_permutation call.

template<PhysicsConfig Cfg>
void reshuffle_common(ParticleContainer<Cfg>& container,
                      const idx_t* perm,
                      count_t n,
                      void* scratch) noexcept
{
    if (n == 0) return;
    container.for_each_common_array([&](auto* base, count_t elems, const char*) {
        (void)elems;
        apply_permutation(base, perm, n, scratch);
    });
}

// =============================================================================
// do_reshuffle_by_key — full driver from SortHelpers to applied permutation
// =============================================================================
// The v2 reference implementation of the PH-key sort uses CPU std::sort on a
// helper buffer of (key, type, original_idx) triples. A GPU radix-sort swap-in
// does not require any change to this interface — helpers are the same 16-byte
// struct, and the permutation is extracted the same way.

template<PhysicsConfig Cfg>
count_t do_reshuffle_by_key(ParticleContainer<Cfg>& container,
                            SortHelper* helpers,   // length == n_part
                            idx_t*      perm,      // output, length == n_part
                            void*       scratch)
{
    const count_t n = container.n_part();
    if (n == 0) return 0;

    // 1. Pull keys and types from Core into the helper buffer.
    const PCoreB* core = container.core();
    for (count_t i = 0; i < n; ++i) {
        helpers[i].key          = core[i].key;
        helpers[i].type         = core[i].type;
        helpers[i].original_idx = static_cast<idx_t>(i);
    }

    // 2. Sort helpers with the shared deterministic comparator.
    PHKeySorter::sort_by_key(helpers, n);

    // 3. Extract the permutation.
    PHKeySorter::extract_permutation(helpers, n, perm, nullptr);

    // 4. Apply to all common arrays via the single registry.
    reshuffle_common(container, perm, n, scratch);

    return n;
}

} // namespace opg::layout_b_coarse

#endif // OPG_LAYOUT_B_COARSE_RESHUFFLE_HPP
