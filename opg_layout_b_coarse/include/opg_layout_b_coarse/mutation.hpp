/**
 * @file mutation.hpp  (layout B, coarse)
 * @brief Gas -> star mutation for layout B.
 *
 * In B the only ordering that must be repaired immediately is the leaf-local
 * ordering of the COMMON arrays. Type-specific slots are reached through
 * linkage and may accumulate holes until compaction.
 */

#ifndef OPG_LAYOUT_B_COARSE_MUTATION_HPP
#define OPG_LAYOUT_B_COARSE_MUTATION_HPP

#include <cstring>

#include "opg_common/tree/box_leaf_base.hpp"
#include "particle_container.hpp"

namespace opg::layout_b_coarse {

using namespace opg::common;

template<PhysicsConfig Cfg>
static inline void move_common_entry_right(ParticleContainer<Cfg>& c,
                                           idx_t from,
                                           idx_t to) noexcept
{
    for (idx_t i = from; i < to; ++i) {
        c.common_registry().swap_all(i, i + 1);
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
    if (c.n_star() >= c.capacity().n_max_star) return -1;

    const idx_t j_common = leaf.type_begin(ParticleType::Gas) + j_within_leaf;
    if (c.core()[j_common].type != static_cast<uint8>(ParticleType::Gas)) return -1;

    const idx_t old_gas_slot = static_cast<idx_t>(c.linkage()[j_common].type_idx);
    const idx_t new_star_slot = static_cast<idx_t>(c.n_star());

    auto& star = c.star_all()[new_star_slot];
    std::memset(&star, 0, sizeof(star));
    c.set_count_star(c.n_star() + 1);

    c.core()[j_common].type = static_cast<uint8>(ParticleType::Star);
    c.linkage()[j_common].type_idx = static_cast<uint32_t>(new_star_slot);
    (void)old_gas_slot; // orphan reclaimed later by compaction

    leaf.type_count[Gas_i]  -= 1;
    leaf.type_count[Star_i] += 1;
    recompute_type_offsets(leaf);

    const idx_t target = leaf.type_begin(ParticleType::Star);
    move_common_entry_right(c, j_common, target);

    if (leaf.type_count[Gas_i] == 0) {
        leaf.flags &= ~leaf_flags::HasGas;
    }
    leaf.flags |= leaf_flags::HasStar;

    return new_star_slot;
}

} // namespace opg::layout_b_coarse

#endif // OPG_LAYOUT_B_COARSE_MUTATION_HPP
