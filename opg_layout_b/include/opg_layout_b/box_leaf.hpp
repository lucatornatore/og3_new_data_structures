/**
 * @file box_leaf.hpp  (layout B')
 * @brief B' box-leaf descriptor (no extra per-type bases) + PLinkage.
 *
 * Under B', each common index j is connected to its type-specific slot via
 * an explicit cross-index held in PLinkage[j]. The leaf descriptor therefore
 * needs no extra per-type array-base fields.
 */

#ifndef OPG_LAYOUT_B_BOX_LEAF_HPP
#define OPG_LAYOUT_B_BOX_LEAF_HPP

#include <opg_common/tree/box_leaf_base.hpp>

namespace opg::layout_b {

struct BoxLeaf : public opg::common::BoxLeafBase {};

static_assert(std::is_trivially_copyable_v<BoxLeaf>);
static_assert(std::is_standard_layout_v<BoxLeaf>);

// =============================================================================
// PLinkage — explicit cross-index, B'-specific
// =============================================================================
//
// PLinkage[j].type_idx gives the slot in the type-specific array for
// Common[j]. The meaning depends on Common[j].type:
//   type==Gas  → index into Gas*
//   type==Star → index into Star*
//   type==BH   → index into BH*

struct PLinkage {
    opg::common::idx_t type_idx;
};

static_assert(std::is_trivially_copyable_v<PLinkage>);
static_assert(std::is_standard_layout_v<PLinkage>);

} // namespace opg::layout_b

#endif // OPG_LAYOUT_B_BOX_LEAF_HPP
