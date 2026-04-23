/**
 * @file compact.hpp  (layout B, coarse)
 * @brief Periodic compaction of type-specific arrays under layout B.
 *
 * Mutations (gas→star, BH mergers) leave orphaned slots in the type-specific
 * arrays. Under B this is tolerable — the old slots are just unused memory,
 * not incorrect data — but left uncompacted they inflate array sizes over
 * time and degrade cache behaviour.
 *
 * compact_type_arrays() scans a container's common arrays, builds an
 * occupancy bitmap over each type-specific array (a slot is "live" iff at
 * least one common entry has linkage[j].type_idx == that slot and
 * core[j].type matches), and rewrites the type-specific arrays in place to
 * remove holes. Linkage entries are updated accordingly.
 *
 * Cost: O(n_part) scan + O(n_gas + n_star + n_bh) rewrite. Runs much less
 * often than PH-reshuffle — a reasonable cadence is "every N full rebuilds"
 * or "when orphan fraction exceeds X%".
 */

#ifndef OPG_LAYOUT_B_COARSE_COMPACT_HPP
#define OPG_LAYOUT_B_COARSE_COMPACT_HPP

#include <vector>
#include <cstring>

#include "particle_container.hpp"

namespace opg::layout_b_coarse {

using namespace opg::common;

// Internal helper: build a forward-map (old_idx -> new_idx) for one type.
// Returns new_count. Slots referenced by at least one live common entry
// get new_map[old_idx] = their position in the compacted array, in order.
// Orphans get new_map[old_idx] = -1.
inline idx_t build_forward_map_for_type(
    const PCoreB*   core,
    const PLinkageB* linkage,
    count_t         n_part,
    ParticleType    target_type,
    count_t         n_old,
    idx_t*          new_map /* size n_old */) noexcept
{
    // First pass: mark live slots with 1 (order-independent).
    for (count_t i = 0; i < n_old; ++i) new_map[i] = -1;
    const uint8_t bare = (uint8_t)target_type;
    for (count_t j = 0; j < n_part; ++j) {
        if ((core[j].type & 0x7Fu) != bare) continue;
        const uint32_t old = linkage[j].type_idx;
        if (old < n_old) new_map[old] = 1;   // mark as live
    }
    // Second pass: assign compacted indices in array order.
    idx_t next = 0;
    for (count_t i = 0; i < n_old; ++i) {
        if (new_map[i] == 1) new_map[i] = next++;
    }
    return next;
}

template<PhysicsConfig Cfg>
struct CompactStats {
    count_t gas_before  = 0, gas_after  = 0;
    count_t star_before = 0, star_after = 0;
    count_t bh_before   = 0, bh_after   = 0;
};

template<PhysicsConfig Cfg>
CompactStats<Cfg> compact_type_arrays(ParticleContainer<Cfg>& c)
{
    CompactStats<Cfg> stats;

    const count_t  N       = c.n_part();
    const PCoreB*  core    = c.core();
    PLinkageB*     linkage = c.linkage();

    // ---- Gas ----
    {
        const count_t n_old = c.n_gas();
        stats.gas_before = n_old;
        std::vector<idx_t> fmap(n_old);
        const idx_t n_new = build_forward_map_for_type(
            core, linkage, N, ParticleType::Gas, n_old, fmap.data());

        // Rewrite gas_all in place: move live entries to their new positions.
        // Safe because new_idx <= old_idx (compacted array is a left-pack).
        auto* gas_all = c.gas_all();
        for (count_t old = 0; old < n_old; ++old) {
            const idx_t nw = fmap[old];
            if (nw < 0) continue;
            if ((idx_t)old != nw) {
                std::memcpy(&gas_all[nw], &gas_all[old], sizeof(gas_all[0]));
            }
        }

        // Rewrite linkage for all Gas common entries.
        for (count_t j = 0; j < N; ++j) {
            if ((core[j].type & 0x7Fu) != (uint8_t)ParticleType::Gas) continue;
            const uint32_t oldi = linkage[j].type_idx;
            linkage[j].type_idx = (uint32_t)fmap[oldi];
        }

        c.set_count_gas((count_t)n_new);
        stats.gas_after = (count_t)n_new;
    }

    // ---- Star ----
    if constexpr (HasStellarEvolution<Cfg>) {
        const count_t n_old = c.n_star();
        stats.star_before = n_old;
        std::vector<idx_t> fmap(n_old);
        const idx_t n_new = build_forward_map_for_type(
            core, linkage, N, ParticleType::Star, n_old, fmap.data());

        auto* star_all = c.star_all();
        for (count_t old = 0; old < n_old; ++old) {
            const idx_t nw = fmap[old];
            if (nw < 0) continue;
            if ((idx_t)old != nw) {
                std::memcpy(&star_all[nw], &star_all[old], sizeof(star_all[0]));
            }
        }
        for (count_t j = 0; j < N; ++j) {
            if ((core[j].type & 0x7Fu) != (uint8_t)ParticleType::Star) continue;
            const uint32_t oldi = linkage[j].type_idx;
            linkage[j].type_idx = (uint32_t)fmap[oldi];
        }
        c.set_count_star((count_t)n_new);
        stats.star_after = (count_t)n_new;
    }

    // ---- BH ----
    if constexpr (HasBlackHoles<Cfg>) {
        const count_t n_old = c.n_bh();
        stats.bh_before = n_old;
        std::vector<idx_t> fmap(n_old);
        const idx_t n_new = build_forward_map_for_type(
            core, linkage, N, ParticleType::BH, n_old, fmap.data());

        auto* bh_all = c.bh_all();
        for (count_t old = 0; old < n_old; ++old) {
            const idx_t nw = fmap[old];
            if (nw < 0) continue;
            if ((idx_t)old != nw) {
                std::memcpy(&bh_all[nw], &bh_all[old], sizeof(bh_all[0]));
            }
        }
        for (count_t j = 0; j < N; ++j) {
            if ((core[j].type & 0x7Fu) != (uint8_t)ParticleType::BH) continue;
            const uint32_t oldi = linkage[j].type_idx;
            linkage[j].type_idx = (uint32_t)fmap[oldi];
        }
        c.set_count_bh((count_t)n_new);
        stats.bh_after = (count_t)n_new;
    }

    return stats;
}

} // namespace opg::layout_b_coarse

#endif // OPG_LAYOUT_B_COARSE_COMPACT_HPP
