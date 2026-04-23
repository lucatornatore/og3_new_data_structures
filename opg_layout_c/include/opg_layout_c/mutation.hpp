/**
 * @file mutation.hpp  (layout C)
 * @brief Staged gas->star mutation for C, preserving positional alignment.
 *
 * The two invariants that matter here are:
 *   1. the staged common entry must leave the live gas sub-range; and
 *   2. every registered gas-side array must follow the same in-leaf swap,
 *      automatically, so adding a new gas block cannot silently break C.
 */

#ifndef OPG_LAYOUT_C_MUTATION_HPP
#define OPG_LAYOUT_C_MUTATION_HPP

#include "particle_container.hpp"
#include "box_leaf.hpp"
#include <cstring>

namespace opg::layout_c {

using namespace opg::common;

namespace staged_type {
    inline constexpr uint8 StagedBit = 0x80;

    constexpr uint8 mark(uint8 t) noexcept   { return static_cast<uint8>(t | StagedBit); }
    constexpr bool  is_staged(uint8 t) noexcept { return (t & StagedBit) != 0; }
    constexpr uint8 unmark(uint8 t) noexcept { return static_cast<uint8>(t & ~StagedBit); }
}

template<PhysicsConfig Cfg>
idx_t stage_gas_to_star(ParticleContainer<Cfg>& container,
                        BoxLeaf& leaf,
                        idx_t j_within_leaf) noexcept
    requires (HasStarFormation<Cfg> && HasHydro<Cfg>)
{
    if (j_within_leaf >= leaf.type_count[static_cast<int>(ParticleType::Gas)]) {
        return -1;
    }

    const idx_t gas_type_idx_in_common = leaf.type_begin(ParticleType::Gas) + j_within_leaf;
    const idx_t gas_array_idx          = leaf.gas_array_begin + j_within_leaf;

    const count_t n_star_before = container.count_star();
    if (n_star_before >= container.capacity().n_max_star) return -1;
    const idx_t star_staging_slot = static_cast<idx_t>(n_star_before);

    auto* star_core = container.star_core();
    auto* star_meta = container.star_meta();
    std::memset(&star_core[star_staging_slot], 0, sizeof(*star_core));
    std::memset(&star_meta[star_staging_slot], 0, sizeof(*star_meta));
    container.set_count_star(n_star_before + 1);

    auto* core = container.core();
    core[gas_type_idx_in_common].type  = staged_type::mark(core[gas_type_idx_in_common].type);
    core[gas_type_idx_in_common].flags |= particle_flags::PendingMutation;

    const idx_t last_gas_in_common = leaf.type_begin(ParticleType::Gas)
                                   + leaf.type_count[static_cast<int>(ParticleType::Gas)] - 1;
    if (gas_type_idx_in_common != last_gas_in_common) {
        container.registry_common().swap_all(gas_type_idx_in_common, last_gas_in_common);
    }

    const idx_t last_gas_in_gasarr = leaf.gas_array_begin
                                   + leaf.type_count[static_cast<int>(ParticleType::Gas)] - 1;
    if (gas_array_idx != last_gas_in_gasarr) {
        container.registry_gas().swap_all(gas_array_idx, last_gas_in_gasarr);
    }

    leaf.type_count[static_cast<int>(ParticleType::Gas)] -= 1;
    leaf.staged_star_count += 1;
    leaf.flags |= leaf_flags::HasStaged;
    if (leaf.type_count[static_cast<int>(ParticleType::Gas)] == 0) {
        leaf.flags &= ~leaf_flags::HasGas;
    }
    leaf.flags |= leaf_flags::HasStar;

    return star_staging_slot;
}

template<PhysicsConfig Cfg>
void integrate_staged_mutations(ParticleContainer<Cfg>& /*container*/,
                                BoxLeaf* /*leaves*/,
                                count_t /*n_leaves*/) noexcept
{
    // TODO: Phase 2.
}

} // namespace opg::layout_c

#endif // OPG_LAYOUT_C_MUTATION_HPP
