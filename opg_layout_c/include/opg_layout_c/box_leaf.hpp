/**
 * @file box_leaf.hpp  (layout C)
 * @brief C box-leaf: adds per-type array-base fields for implicit positional
 *        mapping, plus staged_{star,bh}_count for pending mutations.
 */

#ifndef OPG_LAYOUT_C_BOX_LEAF_HPP
#define OPG_LAYOUT_C_BOX_LEAF_HPP

#include <opg_common/tree/box_leaf_base.hpp>

namespace opg::layout_c {

struct BoxLeaf : public opg::common::BoxLeafBase {
    // Per-type array-base indices (C-specific positional mapping).
    opg::common::idx_t gas_array_begin  = 0;
    opg::common::idx_t star_array_begin = 0;
    opg::common::idx_t bh_array_begin   = 0;

    // Mutation staging: pending stars/BHs created from gas in this leaf, not
    // yet integrated into the proper star/bh ranges (integrated at next rebuild).
    opg::common::idx_t staged_star_count = 0;
    opg::common::idx_t staged_bh_count   = 0;

    // Convenience positional accessors.
    constexpr opg::common::idx_t gas_slot_for(opg::common::idx_t j_within_leaf) const noexcept {
        return gas_array_begin + j_within_leaf;
    }
    constexpr opg::common::idx_t star_slot_for(opg::common::idx_t j_within_leaf) const noexcept {
        return star_array_begin + j_within_leaf;
    }
    constexpr opg::common::idx_t bh_slot_for(opg::common::idx_t j_within_leaf) const noexcept {
        return bh_array_begin + j_within_leaf;
    }
};

static_assert(std::is_trivially_copyable_v<BoxLeaf>);
// NOTE: not standard-layout because it inherits from BoxLeafBase AND adds its
// own non-static data members. Trivial copyability is what we need for memcpy
// in MPI pack/unpack and scratch-buffer moves; standard-layout is not required.

} // namespace opg::layout_c

#endif // OPG_LAYOUT_C_BOX_LEAF_HPP
