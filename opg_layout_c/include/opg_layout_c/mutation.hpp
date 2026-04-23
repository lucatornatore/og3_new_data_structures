#ifndef OPG_LAYOUT_C_MUTATION_HPP
#define OPG_LAYOUT_C_MUTATION_HPP

#include "particle_container.hpp"
#include "box_leaf.hpp"
#include <cstring>

namespace opg::layout_c {

using namespace opg::common;

namespace staged_type {
    inline constexpr uint8 StagedBit = 0x80;
    constexpr uint8 mark(uint8 t)   noexcept { return static_cast<uint8>(t | StagedBit); }
    constexpr bool  is_staged(uint8 t) noexcept { return (t & StagedBit) != 0; }
    constexpr uint8 unmark(uint8 t) noexcept { return static_cast<uint8>(t & ~StagedBit); }
}

template<PhysicsConfig Cfg>
idx_t stage_gas_to_star(ParticleContainer<Cfg>& container, BoxLeaf& leaf, idx_t j_within_leaf) noexcept
    requires HasStellarEvolution<Cfg> && HasStarFormation<Cfg>
{
    const idx_t gas_type_idx_in_common = leaf.common_begin + leaf.type_offset[static_cast<int>(ParticleType::Gas)] + j_within_leaf;
    const idx_t gas_array_idx = leaf.gas_array_begin + j_within_leaf;

    const count_t n_star_before = container.count_star();
    const idx_t star_staging_slot = static_cast<idx_t>(n_star_before);

    auto* star_core = container.star_core();
    auto* star_meta = container.star_meta();
    std::memset(&star_core[star_staging_slot], 0, sizeof(*star_core));
    std::memset(&star_meta[star_staging_slot], 0, sizeof(*star_meta));
    container.set_count_star(n_star_before + 1);

    auto* core = container.core();
    core[gas_type_idx_in_common].type  = staged_type::mark(core[gas_type_idx_in_common].type);
    
    // FIX: Add pending mutation flag. Physics loops must mask this out!
    core[gas_type_idx_in_common].flags |= particle_flags::PendingMutation;

    const idx_t last_gas_in_common = leaf.common_begin + leaf.type_offset[static_cast<int>(ParticleType::Gas)] + leaf.type_count [static_cast<int>(ParticleType::Gas)] - 1;
    
    // FIX 7: Use Automated swap_all instead of manual logic
    if (gas_type_idx_in_common != last_gas_in_common) {
        container.registry_common().swap_all(gas_type_idx_in_common, last_gas_in_common);
    }

    const idx_t last_gas_in_gasarr = leaf.gas_array_begin + leaf.type_count[static_cast<int>(ParticleType::Gas)] - 1;
    
    if (gas_array_idx != last_gas_in_gasarr) {
        container.registry_gas().swap_all(gas_array_idx, last_gas_in_gasarr);
    }

    leaf.type_count[static_cast<int>(ParticleType::Gas)] -= 1;
    leaf.staged_star_count += 1;
    leaf.flags |= leaf_flags::HasStaged;

    return star_staging_slot;
}

template<PhysicsConfig Cfg>
void integrate_staged_mutations(ParticleContainer<Cfg>&, BoxLeaf*, count_t) noexcept {}

} // namespace opg::layout_c

#endif // OPG_LAYOUT_C_MUTATION_HPP
