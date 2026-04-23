/**
 * @file mutation.hpp  (layout B')
 * @brief Gas → Star mutation for B'.
 *
 * Under B', mutation is local: the type field in Common[j] flips from Gas
 * to Star, a new slot is taken from the Star arrays (the type-specific
 * region), and PLinkage[j].type_idx is rewritten to point at the new slot.
 * No in-leaf swap is needed because B' does NOT rely on positional alignment.
 *
 * Stub — full implementation deferred to Phase 2; the design is straightforward
 * because B' has no positional invariant to preserve.
 */

#ifndef OPG_LAYOUT_B_MUTATION_HPP
#define OPG_LAYOUT_B_MUTATION_HPP

#include "particle_container.hpp"

namespace opg::layout_b {

using namespace opg::common;

/**
 * Convert Common[j] from gas to star in place. Requires HasStellarEvolution.
 * Returns the slot index in Star[] that now holds the new star's data.
 */
template<PhysicsConfig Cfg>
idx_t mutate_gas_to_star(ParticleContainer<Cfg>& /*container*/,
                         idx_t /*common_idx*/) noexcept
    requires HasStellarEvolution<Cfg>
{
    // 1) Find free slot in Star arrays (n_star < n_max_star).
    // 2) Copy physically meaningful common-layer state is unchanged: it stays
    //    in Common[j] — B' keeps common data regardless of mutation.
    // 3) Flip Common[j].type to Star, increment star count.
    // 4) Rewrite PLinkage[j].type_idx to the new Star slot.
    //
    // Implementation pending Phase 2.
    return -1;
}

} // namespace opg::layout_b

#endif // OPG_LAYOUT_B_MUTATION_HPP
