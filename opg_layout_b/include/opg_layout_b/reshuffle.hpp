/**
 * @file reshuffle.hpp  (layout B')
 * @brief Reshuffle driver for B'.
 *
 * Under B':
 *   - Common arrays (and all registered parallel common arrays) are permuted
 *     by PH-key through a single ArrayRegistry.
 *   - Type-specific arrays (Gas[], Star[], BH[]) do NOT move here; they are
 *     compacted separately via compact_type_arrays() when needed. PLinkage
 *     values are rewritten after that compaction to point at the new slots.
 *
 * This is the "simpler reshuffle" property of B' noted in the design doc.
 */

#ifndef OPG_LAYOUT_B_RESHUFFLE_HPP
#define OPG_LAYOUT_B_RESHUFFLE_HPP

#include <opg_common/permutation/permutation.hpp>
#include "particle_container.hpp"

namespace opg::layout_b {

using namespace opg::common;

/**
 * Apply a PH-key permutation to all common arrays.
 *
 * @param container container whose registry holds the common arrays
 * @param perm      forward permutation (perm[i] = old index)
 * @param scratch   scratch buffer ≥ registry.max_elem_bytes()
 *
 * Type-specific arrays are unchanged. PLinkage is itself a common array, so
 * it is permuted along with the others — its values (indices into Gas/Star/BH)
 * remain valid because those arrays haven't moved.
 */
template<PhysicsConfig Cfg>
void reshuffle_common(ParticleContainer<Cfg>& container,
                      const idx_t* perm,
                      void*        scratch) noexcept
{
    const count_t n = container.count();
    container.registry().reshuffle_all(perm, n, scratch);
}

/**
 * Compact the type-specific arrays and rewrite PLinkage values.
 *
 * After a batch of mutations (gas→star, BH mergers), type arrays develop
 * holes. This function compacts them and rewrites PLinkage[j].type_idx in
 * Common to point at the new slots. Stub for Phase 1.
 */
template<PhysicsConfig Cfg>
void compact_type_arrays(ParticleContainer<Cfg>& /*container*/) noexcept {
    // TODO: implement together with the B' mutation integration step.
}

} // namespace opg::layout_b

#endif // OPG_LAYOUT_B_RESHUFFLE_HPP
