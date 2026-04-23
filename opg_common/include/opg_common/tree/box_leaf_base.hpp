/**
 * @file box_leaf_base.hpp
 * @brief Base box-leaf descriptor, shared by B' and C prototypes.
 *
 * B' adds nothing (cross-indices live in PLinkage[]).
 * C adds gas_array_begin / star_array_begin / bh_array_begin +
 * staged_star_count / staged_bh_count, for positional mapping and staging.
 */

#ifndef OPG_COMMON_TREE_BOX_LEAF_BASE_HPP
#define OPG_COMMON_TREE_BOX_LEAF_BASE_HPP

#include "../types/scalar_types.hpp"
#include "../types/physics_config.hpp"
#include <type_traits>

namespace opg::common {

namespace leaf_flags {
    inline constexpr uint32 Active      = 1u << 0;
    inline constexpr uint32 NeedsUpdate = 1u << 1;
    inline constexpr uint32 Boundary    = 1u << 2;
    inline constexpr uint32 HasGas      = 1u << 3;
    inline constexpr uint32 HasStar     = 1u << 4;
    inline constexpr uint32 HasBH       = 1u << 5;
    inline constexpr uint32 HasDM       = 1u << 6;
    inline constexpr uint32 HasStaged   = 1u << 7;   // pending mutations
}

struct BoxLeafBase {
    // Key range
    pkey_t  key_begin;
    pkey_t  key_end;

    // Geometry
    pos3_t  center;
    real_t  half_width;

    // Barycenter (modified Barnes-Hut)
    pos3_t  com;
    mass_t  total_mass;

    // Common-array range
    idx_t   common_begin;
    idx_t   common_end;

    // Per-type sub-ranges inside [common_begin, common_end)
    idx_t   type_offset[NTYPES];
    idx_t   type_count [NTYPES];

    // Tree structure
    idx_t   parent;
    idx_t   level;

    // Auxiliary
    real_t  h_max;
    real_t  h_min;
    real_t  v_max;
    time_int_t ti_lastkicked;
    uint32  flags;

    constexpr idx_t count() const noexcept { return common_end - common_begin; }

    constexpr idx_t type_begin(ParticleType pt) const noexcept {
        return common_begin + type_offset[static_cast<int>(pt)];
    }
    constexpr idx_t type_end(ParticleType pt) const noexcept {
        return type_begin(pt) + type_count[static_cast<int>(pt)];
    }

    constexpr idx_t gas_count()  const noexcept { return type_count[static_cast<int>(ParticleType::Gas)];  }
    constexpr idx_t star_count() const noexcept { return type_count[static_cast<int>(ParticleType::Star)]; }
    constexpr idx_t bh_count()   const noexcept { return type_count[static_cast<int>(ParticleType::BH)];   }
};

static_assert(std::is_trivially_copyable_v<BoxLeafBase>);
static_assert(std::is_standard_layout_v<BoxLeafBase>);

inline void recompute_type_offsets(BoxLeafBase& leaf) noexcept {
    idx_t offset = 0;
    for (int t = 0; t < NTYPES; ++t) {
        leaf.type_offset[t] = offset;
        offset += leaf.type_count[t];
    }
    leaf.common_end = leaf.common_begin + offset;
}

// =============================================================================
// SortHelper — the sortable (key, type, origin) triplet
// =============================================================================

struct SortHelper {
    pkey_t  key;
    uint8   type;
    uint8   _pad[3];
    idx_t   original_idx;
};

static_assert(sizeof(SortHelper) == 16);
static_assert(std::is_trivially_copyable_v<SortHelper>);

} // namespace opg::common

#endif // OPG_COMMON_TREE_BOX_LEAF_BASE_HPP
