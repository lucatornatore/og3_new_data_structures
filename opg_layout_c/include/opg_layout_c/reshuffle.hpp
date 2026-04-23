/**
 * @file reshuffle.hpp  (layout C)
 * @brief Authoritative reshuffle driver for C — four permutations, four groups.
 *
 * The central and most error-prone operation in Option C (design note v4 § 8).
 * It must permute, atomically and in lockstep:
 *   - Every common array  by perm_common.
 *   - Every gas    array  by perm_gas.
 *   - Every star   array  by perm_star.
 *   - Every bh     array  by perm_bh.
 *
 * The four permutations are derived from a single sorted SortHelper buffer
 * so that they are internally consistent. If any component array is not in
 * its length-matched registry, the positional mapping breaks and every
 * subsequent access becomes silent data corruption — hence the N2 discipline
 * of declaring registration exactly once, inside `register_all_arrays()`.
 */

#ifndef OPG_LAYOUT_C_RESHUFFLE_HPP
#define OPG_LAYOUT_C_RESHUFFLE_HPP

#include <opg_common/permutation/permutation.hpp>
#include "particle_container.hpp"

namespace opg::layout_c {

using namespace opg::common;

// =============================================================================
// PermutationBundle — the four permutations derived from one sort
// =============================================================================
//
// Caller owns the memory; typically these are arena-allocated scratch.
// perm_common is of length n_common; the per-type permutations are "local"
// permutations of the respective type arrays (length n_gas, n_star, n_bh).

struct PermutationBundle {
    idx_t*  perm_common = nullptr;
    idx_t*  perm_gas    = nullptr;
    idx_t*  perm_star   = nullptr;
    idx_t*  perm_bh     = nullptr;

    count_t n_common = 0;
    count_t n_gas    = 0;
    count_t n_star   = 0;
    count_t n_bh     = 0;
};

// =============================================================================
// derive_permutations
// =============================================================================
//
// Given a SortHelper array *already sorted by PH key*, and an auxiliary array
// `type_idx_before_sort[i]` that maps the i-th particle in the PRE-sort
// common order to its position within its type-array (i.e., the gas-array
// slot if type=Gas, star-array slot if type=Star, bh-array slot if type=BH),
// fill the four permutations.
//
// The auxiliary `type_idx_before_sort` must be built alongside the
// SortHelper when the sort is initiated — it is what lets us, after the key
// sort, recover the original type-array positions. In practice this array is
// a parallel companion to the SortHelper produced at helper-build time.
//
// For DM types, there is no type array, so no entries in perm_gas/star/bh.

template<PhysicsConfig Cfg>
void derive_permutations(const SortHelper* sorted_helpers,
                         const idx_t*      type_idx_before_sort,
                         count_t           n_total,
                         PermutationBundle& out) noexcept
{
    count_t ig = 0, is = 0, ib = 0;

    for (count_t i = 0; i < n_total; ++i) {
        const SortHelper& h = sorted_helpers[i];

        // Common permutation: position i in the sorted order was originally
        // at h.original_idx.
        out.perm_common[i] = h.original_idx;

        const ParticleType pt = static_cast<ParticleType>(h.type);
        const idx_t orig_type_slot = type_idx_before_sort
                                      ? type_idx_before_sort[h.original_idx]
                                      : idx_t{-1};

        if (pt == ParticleType::Gas) {
            if (out.perm_gas) out.perm_gas[ig] = orig_type_slot;
            ++ig;
        } else if constexpr (HasStellarEvolution<Cfg>) {
            if (pt == ParticleType::Star) {
                if (out.perm_star) out.perm_star[is] = orig_type_slot;
                ++is;
            }
        }
        if constexpr (HasBlackHoles<Cfg>) {
            if (pt == ParticleType::BH) {
                if (out.perm_bh) out.perm_bh[ib] = orig_type_slot;
                ++ib;
            }
        }
    }

    out.n_common = n_total;
    out.n_gas    = ig;
    out.n_star   = is;
    out.n_bh     = ib;
}

// =============================================================================
// reshuffle_all — apply the bundle to the container's four registries
// =============================================================================
//
// `scratch_*` are per-registry scratch buffers sized ≥ max_elem_bytes()
// of their respective registry. In practice the container can expose
// `registry_gas().max_elem_bytes()` so the caller allocates exactly enough.

template<PhysicsConfig Cfg>
void reshuffle_all(ParticleContainer<Cfg>&  container,
                   const PermutationBundle& bundle,
                   void* scratch_common,
                   void* scratch_gas,
                   void* scratch_star,
                   void* scratch_bh) noexcept
{
    if (bundle.n_common > 0) {
        container.for_each_common_array([&](auto* base, count_t elems, const char*) {
            (void)elems;
            apply_permutation(base, bundle.perm_common, bundle.n_common, scratch_common);
        });
    }

    if constexpr (HasHydro<Cfg>) {
        if (bundle.n_gas > 0) {
            container.for_each_gas_array([&](auto* base, count_t elems, const char*) {
                (void)elems;
                apply_permutation(base, bundle.perm_gas, bundle.n_gas, scratch_gas);
            });
        }
    }

    if constexpr (HasStellarEvolution<Cfg>) {
        if (bundle.n_star > 0) {
            container.for_each_star_array([&](auto* base, count_t elems, const char*) {
                (void)elems;
                apply_permutation(base, bundle.perm_star, bundle.n_star, scratch_star);
            });
        }
    }

    if constexpr (HasBlackHoles<Cfg>) {
        if (bundle.n_bh > 0) {
            container.for_each_bh_array([&](auto* base, count_t elems, const char*) {
                (void)elems;
                apply_permutation(base, bundle.perm_bh, bundle.n_bh, scratch_bh);
            });
        }
    }
}

} // namespace opg::layout_c

#endif // OPG_LAYOUT_C_RESHUFFLE_HPP
