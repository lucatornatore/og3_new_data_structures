/**
 * @file example_layout_b.cpp
 * @brief Minimal demonstration of the B' prototype.
 *
 * Shows:
 *   - Arena setup
 *   - Container construction (NTTP PhysicsConfig)
 *   - Accessor guarantees (compile-time-available vs. rejected)
 *   - Registry contents
 *   - A trivial reshuffle on a small synthetic dataset
 */

#include <opg_layout_b/particle_container.hpp>
#include <opg_layout_b/reshuffle.hpp>

#include <cstdio>

using namespace opg::common;
using namespace opg::layout_b;

int main() {
    std::printf("=== B' prototype demonstration ===\n");
    std::printf("Layout: %s\n\n",
                ParticleContainer<physics_configs::StandardSPH>::layout_name());

    // ---- arena + capacity --------------------------------------------------
    ArenaConfig acfg{
        .total_particles  = 10'000,
        .n_tasks          = 1,
        .part_alloc_factor = 2.0,
        .alignment        = 64
    };
    std::printf("n_average = %llu,  n_max = %llu\n",
                (unsigned long long)acfg.n_average(),
                (unsigned long long)acfg.n_max());

    MemoryArena arena(128 * 1024 * 1024, acfg);

    Capacity cap{.n_max=20'000, .n_max_gas=10'000, .n_max_star=5'000, .n_max_bh=100};

    // ---- StandardSPH container --------------------------------------------
    ParticleContainer<physics_configs::StandardSPH> container(arena, cap);

    std::printf("\nStandardSPH container:\n");
    std::printf("  allocated bytes : %zu (%.1f MB)\n",
                container.total_bytes_allocated(),
                container.total_bytes_allocated() / (1024.0*1024.0));
    std::printf("  registry size   : %zu  (common arrays; type arrays don't move under B')\n",
                container.registry().size());

    // Demonstrate: optional accessors that ARE present at compile time.
    PLeap*       leap  = container.leap();        // HasLeapfrog → present
    auto*        gmet  = container.gas_metal();   // HasMetals   → present
    auto*        scor  = container.star_core();   // HasStellarEvolution → present
    auto*        bhc   = container.bh_core();     // HasBlackHoles → present
    (void)leap; (void)gmet; (void)scor; (void)bhc;
    std::printf("  leap/gas_metal/star_core/bh_core accessors compiled (as expected)\n");

    // ---- GravityOnly container: disabled features are absent ----------------
    ParticleContainer<physics_configs::GravityOnly> g_only(arena, cap);
    std::printf("\nGravityOnly container:\n");
    std::printf("  allocated bytes : %zu (%.1f MB)\n",
                g_only.total_bytes_allocated(),
                g_only.total_bytes_allocated() / (1024.0*1024.0));
    std::printf("  registry size   : %zu\n", g_only.registry().size());
    // g_only.gas_mag(), g_only.leap(), etc. would FAIL TO COMPILE — not runtime null.
    std::printf("  (gas_mag/leap/star_core are absent at compile time for GravityOnly)\n");

    // ---- trivial reshuffle on 5 particles ----------------------------------
    container.set_counts(5, 3, 2, 0);
    auto* core = container.core();
    for (int i = 0; i < 5; ++i) {
        core[i].key    = (pkey_t)(100 - 10*i);  // descending — will be sorted ascending
        core[i].mass   = 1.0f + i;
        core[i].type   = (uint8)((i < 3) ? ParticleType::Gas : ParticleType::Star);
    }

    // Build sort helpers, sort, extract perm.
    SortHelper helpers[5];
    uint8 types[5];
    pkey_t keys[5];
    for (int i=0;i<5;++i){ keys[i] = core[i].key; types[i] = core[i].type; }
    PHKeySorter::build_helpers(keys, types, 5, helpers);
    PHKeySorter::sort_by_key (helpers, 5);
    idx_t perm[5];
    PHKeySorter::extract_permutation(helpers, 5, perm, nullptr);

    std::printf("\nBefore reshuffle: keys = ");
    for (int i=0;i<5;++i) std::printf("%llu ", (unsigned long long)core[i].key);
    std::printf("\nPermutation     : ");
    for (int i=0;i<5;++i) std::printf("%d ", (int)perm[i]);

    // Scratch buffer from the arena temp stack: largest element size × active particles.
    ArenaCheckpoint scratch_scope(arena);
    const size_t scratch_bytes = container.registry().max_elem_size() * 5;
    void* scratch = arena.allocate_temp_bytes(scratch_bytes, arena.config().alignment).ptr;

    reshuffle_common(container, perm, scratch);

    std::printf("\nAfter reshuffle : keys = ");
    for (int i=0;i<5;++i) std::printf("%llu ", (unsigned long long)core[i].key);
    std::printf("\n");

    return 0;
}
