/**
 * @file mutation.hpp  (layout B')
 * @brief Gas -> star mutation for B'.
 *
 * B' uses an explicit cross-index (`PLinkage`) like coarse B, but retains the
 * fine common split. Therefore the immediate repair step is: rewrite the
 * common entry and linkage, then move that COMMON entry to the star run inside
 * the leaf. Type-specific gas slots are left orphaned until compaction.
 */

#ifndef OPG_LAYOUT_B_MUTATION_HPP
#define OPG_LAYOUT_B_MUTATION_HPP

#include <cstring>

#include "particle_container.hpp"
#include "box_leaf.hpp"

namespace opg::layout_b {

using namespace opg::common;

template<PhysicsConfig Cfg>
static inline void move_common_entry_right(ParticleContainer<Cfg>& c,
                                           idx_t from,
                                           idx_t to) noexcept
{
    for (idx_t i = from; i < to; ++i) {
        c.registry().swap_all(i, i + 1);
    }
}

template<PhysicsConfig Cfg>
idx_t mutate_gas_to_star(ParticleContainer<Cfg>& c,
                         BoxLeafBase& leaf,
                         idx_t j_within_leaf) noexcept
    requires (HasStarFormation<Cfg> && HasHydro<Cfg>)
{
    const int Gas_i  = static_cast<int>(ParticleType::Gas);
    const int Star_i = static_cast<int>(ParticleType::Star);

    if (j_within_leaf >= leaf.type_count[Gas_i]) return -1;
    if (c.count_star() >= c.capacity().n_max_star) return -1;

    const idx_t j_common = leaf.type_begin(ParticleType::Gas) + j_within_leaf;
    if (c.core()[j_common].type != static_cast<uint8>(ParticleType::Gas)) return -1;

    const idx_t old_gas_slot = static_cast<idx_t>(c.linkage()[j_common].type_idx);
    const idx_t new_star_slot = static_cast<idx_t>(c.count_star());

    auto& star_core = c.star_core()[new_star_slot];
    auto& star_meta = c.star_meta()[new_star_slot];
    std::memset(&star_core, 0, sizeof(star_core));
    std::memset(&star_meta, 0, sizeof(star_meta));
    c.set_count_star(c.count_star() + 1);

    c.core()[j_common].type = static_cast<uint8>(ParticleType::Star);
    c.linkage()[j_common].type_idx = static_cast<uint32_t>(new_star_slot);
    (void)old_gas_slot;

    leaf.type_count[Gas_i]  -= 1;
    leaf.type_count[Star_i] += 1;
    recompute_type_offsets(leaf);

    const idx_t target = leaf.type_begin(ParticleType::Star);
    move_common_entry_right(c, j_common, target);

    if (leaf.type_count[Gas_i] == 0) {
        leaf.flags &= ~leaf_flags::HasGas;
    }
    leaf.flags |= leaf_flags::HasStar;

    c.set_count_gas(c.count_gas() - 1);
    return new_star_slot;
}

} // namespace opg::layout_b

#endif // OPG_LAYOUT_B_MUTATION_HPP
