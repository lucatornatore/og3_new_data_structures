/**
 * @file mutation.hpp  (layout C)
 * @brief Staged gas→star mutation for C, preserving positional alignment.
 *
 * Protocol (design note v4 § 9.1):
 *
 *   1. Copy Common data of Common[j] to a staging slot past the current
 *      star range: Star[n_star + leaf.staged_star_count]. (The slot lives
 *      outside any leaf's star range until the next rebuild.)
 *   2. Mark Common[j] with PendingMutation and flip its type's high bit
 *      (staged_type sentinel). This keeps the position valid for the
 *      arrays it still occupies while flagging the entry as inactive.
 *   3. In-leaf swap: move the flagged entry to the END of the gas
 *      sub-range within its box-leaf. This preserves positional
 *      alignment for the remaining gas particles — each surviving gas
 *      particle's position in the gas-arrays still matches its position
 *      in the leaf's common range.
 *   4. Shrink the leaf: leaf.type_count[Gas] -= 1; leaf.staged_star_count += 1.
 *   5. SPH loops consume only leaf.gas_count() gas particles, so the
 *      staged entry is skipped WITHOUT a per-particle branch and WITHOUT
 *      GPU thread divergence.
 *   6. At the next full rebuild, integrate_staged_mutations() folds the
 *      staged entries into their proper star ranges. (Stubbed here.)
 *
 * The in-leaf swap is the non-obvious part. Naively marking Common[j] as
 * inactive would break positional alignment for *every subsequent gas
 * access in that leaf* because the surviving gas particles' indices in the
 * gas-arrays would no longer match their indices in the leaf's common
 * range. Swapping to the end restores the alignment for the survivors.
 */

#ifndef OPG_LAYOUT_C_MUTATION_HPP
#define OPG_LAYOUT_C_MUTATION_HPP

#include "particle_container.hpp"
#include "box_leaf.hpp"
#include <cstring>

namespace opg::layout_c {

using namespace opg::common;

namespace staged_type {
    // A staged particle has the high bit of its type field set.
    inline constexpr uint8 StagedBit = 0x80;

    constexpr uint8 mark(uint8 t)   noexcept { return static_cast<uint8>(t | StagedBit); }
    constexpr bool  is_staged(uint8 t) noexcept { return (t & StagedBit) != 0; }
    constexpr uint8 unmark(uint8 t) noexcept { return static_cast<uint8>(t & ~StagedBit); }
}

// =============================================================================
// swap helper — swap common-array entries at two indices
// =============================================================================
//
// Only common-array swap: the flagged particle's gas slot (in the gas arrays)
// is NOT swapped, because the flagged particle is about to become a star and
// its gas slot will be released at the integrate step. During the in-flight
// window, the flagged gas slot is left in place but ignored (the leaf's
// shrunk range means it is out-of-bounds for SPH loops).
//
// This swap must be applied to every REGISTERED common array.

template<PhysicsConfig Cfg>
void swap_common_entries(ParticleContainer<Cfg>& container,
                         idx_t a, idx_t b) noexcept
{
    if (a == b) return;
    auto* core = container.core();
    std::swap(core[a], core[b]);
    auto* dyn  = container.dyn();
    std::swap(dyn[a],  dyn[b]);
    auto* time = container.time();
    std::swap(time[a], time[b]);
    auto* meta = container.meta();
    std::swap(meta[a], meta[b]);

    if constexpr (HasLeapfrog<Cfg>) {
        auto* leap = container.leap();
        std::swap(leap[a], leap[b]);
    }
    if constexpr (HasPotentialOutput<Cfg>) {
        auto* pot = container.potential();
        std::swap(pot[a], pot[b]);
    }
}

// =============================================================================
// stage_gas_to_star
// =============================================================================
//
// Stage a gas particle for conversion to a star. Returns the index of the
// new star's staging slot (in star arrays, past the current star range).
//
// `j_within_leaf` is the 0-based index within the leaf's gas sub-range.

template<PhysicsConfig Cfg>
idx_t stage_gas_to_star(ParticleContainer<Cfg>& container,
                        BoxLeaf& leaf,
                        idx_t    j_within_leaf) noexcept
    requires HasStellarEvolution<Cfg>
{
    // Absolute indices.
    const idx_t gas_type_idx_in_common = leaf.common_begin
                                       + leaf.type_offset[static_cast<int>(ParticleType::Gas)]
                                       + j_within_leaf;
    const idx_t gas_array_idx = leaf.gas_array_begin + j_within_leaf;

    // (1) Star staging slot: append past current star range.
    const count_t n_star_before = container.count_star();
    const idx_t   star_staging_slot = static_cast<idx_t>(n_star_before);

    // Copy common data → StarCore[staging_slot] is NOT the move; under C, the
    // staged particle's common fields remain in Common[j] during the
    // in-flight window. At integrate time, the proper StarCore entry is
    // populated from the surviving common fields + the staged payload.
    //
    // We DO record the intent in StarMeta / StarCore of the staging slot so
    // the integrate pass has somewhere to find it. For simplicity here we
    // just zero-initialise the slot and leave the detailed population as a
    // Phase-2 task.
    auto* star_core = container.star_core();
    auto* star_meta = container.star_meta();
    std::memset(&star_core[star_staging_slot], 0, sizeof(*star_core));
    std::memset(&star_meta[star_staging_slot], 0, sizeof(*star_meta));
    container.set_count_star(n_star_before + 1);

    // (2) Mark Common[j] as staged: high bit on type, set flag.
    auto* core = container.core();
    core[gas_type_idx_in_common].type  = staged_type::mark(core[gas_type_idx_in_common].type);
    core[gas_type_idx_in_common].flags |= particle_flags::PendingMutation;

    // (3) In-leaf swap: move staged entry to END of gas sub-range in common.
    const idx_t last_gas_in_common = leaf.common_begin
                                   + leaf.type_offset[static_cast<int>(ParticleType::Gas)]
                                   + leaf.type_count [static_cast<int>(ParticleType::Gas)] - 1;
    if (gas_type_idx_in_common != last_gas_in_common) {
        swap_common_entries(container, gas_type_idx_in_common, last_gas_in_common);
    }

    // The staged entry is now the LAST gas position in the leaf. Shrinking
    // the leaf's gas count by 1 moves it out of the SPH-visible range.
    // (No swap in the gas-arrays: the gas-array slot at gas_array_idx now
    // holds data for a DIFFERENT gas particle — the one just swapped in —
    // so we must swap those too to keep positional alignment for the
    // survivors. This is the correction: the in-leaf swap applies to both
    // common AND gas arrays.)
    const idx_t last_gas_in_gasarr = leaf.gas_array_begin
                                   + leaf.type_count[static_cast<int>(ParticleType::Gas)] - 1;
    if (gas_array_idx != last_gas_in_gasarr) {
        // Swap all registered gas arrays at (gas_array_idx, last_gas_in_gasarr).
        auto* gc = container.gas_core();
        std::swap(gc[gas_array_idx], gc[last_gas_in_gasarr]);

        if constexpr (HasSPH<Cfg>) {
            auto* gg = container.gas_grad();
            std::swap(gg[gas_array_idx], gg[last_gas_in_gasarr]);
        }
        if constexpr (HasMagnetic<Cfg>) {
            auto* gm = container.gas_mag();
            std::swap(gm[gas_array_idx], gm[last_gas_in_gasarr]);
        }
        if constexpr (HasMetals<Cfg>) {
            auto* gmet = container.gas_metal();
            std::swap(gmet[gas_array_idx], gmet[last_gas_in_gasarr]);
        }
        if constexpr (HasStarFormation<Cfg>) {
            auto* gsf = container.gas_sf();
            std::swap(gsf[gas_array_idx], gsf[last_gas_in_gasarr]);
        }
        if constexpr (HasMFM<Cfg>) {
            auto* gmf = container.gas_mfm();
            std::swap(gmf[gas_array_idx], gmf[last_gas_in_gasarr]);
        }
    }

    // (4) Shrink the leaf's gas range; record the staged star count.
    leaf.type_count[static_cast<int>(ParticleType::Gas)] -= 1;
    leaf.staged_star_count += 1;
    leaf.flags |= leaf_flags::HasStaged;

    return star_staging_slot;
}

// =============================================================================
// integrate_staged_mutations — called at next full rebuild (stub)
// =============================================================================
//
// At rebuild time:
//   - Walk the staged star range, populate proper Star* fields from the
//     surviving Common data of the still-staged particles.
//   - Clear the staged sentinel bits on Common[].
//   - During the rebuild's PH-key sort, the (now real) star particles get
//     their proper positions in the new box-leaves' star ranges, and the
//     orphaned gas slots are discarded because their leaves no longer
//     point to them.
//
// Implementation deferred to the rebuild infrastructure.

template<PhysicsConfig Cfg>
void integrate_staged_mutations(ParticleContainer<Cfg>& /*container*/,
                                BoxLeaf*                /*leaves*/,
                                count_t                 /*n_leaves*/) noexcept
{
    // TODO: Phase 2.
}

} // namespace opg::layout_c

#endif // OPG_LAYOUT_C_MUTATION_HPP
