/**
 * @file example_layout_c.cpp
 * @brief Minimal demonstration of the C prototype.
 *
 * Shows:
 *   - Arena setup
 *   - Container construction (NTTP PhysicsConfig)
 *   - Four-registry architecture for lockstep reshuffle
 *   - The staged-mutation protocol: stage_gas_to_star preserves positional
 *     alignment for surviving gas particles WITHOUT a per-particle branch.
 */

#include <opg_layout_c/particle_container.hpp>
#include <opg_layout_c/reshuffle.hpp>
#include <opg_layout_c/mutation.hpp>

#include <cstdio>

using namespace opg::common;
using namespace opg::layout_c;

int main() {
    std::printf("=== C prototype demonstration ===\n");
    std::printf("Layout: %s\n\n",
                ParticleContainer<physics_configs::StandardSPH>::layout_name());

    ArenaConfig acfg{
        .total_particles  = 1'000,
        .n_tasks          = 1,
        .part_alloc_factor = 2.0,
        .alignment        = 64
    };
    MemoryArena arena(64 * 1024 * 1024, acfg);

    Capacity cap{.n_max=2'000, .n_max_gas=1'000, .n_max_star=500, .n_max_bh=10};

    ParticleContainer<physics_configs::StandardSPH> container(arena, cap);

    std::printf("Container total bytes: %zu (%.1f MB)\n",
                container.total_bytes_allocated(),
                container.total_bytes_allocated() / (1024.0*1024.0));
    std::printf("Registries            : common=%zu gas=%zu star=%zu bh=%zu\n\n",
                container.registry_common().size(),
                container.registry_gas().size(),
                container.registry_star().size(),
                container.registry_bh().size());

    // ---- set up a toy leaf with 4 gas + 1 star -----------------------------
    container.set_counts(/*np*/7, /*ng*/4, /*ns*/1, /*nb*/0);

    BoxLeaf leaf{};
    leaf.common_begin = 100;
    leaf.common_end   = 105;
    leaf.type_offset[(int)ParticleType::Gas]  = 0;
    leaf.type_count [(int)ParticleType::Gas]  = 4;
    leaf.type_offset[(int)ParticleType::Star] = 4;
    leaf.type_count [(int)ParticleType::Star] = 1;
    leaf.gas_array_begin  = 0;
    leaf.star_array_begin = 0;

    auto* core = container.core();
    auto* gas  = container.gas_core();
    for (int j=0; j<4; ++j) {
        core[100+j].type = (uint8)ParticleType::Gas;
        core[100+j].mass = 1.0f + j;
        gas[j].hsml      = 0.5f;
        gas[j].density   = 10.0f + j;
    }
    core[104].type = (uint8)ParticleType::Star;
    core[104].mass = 7.0f;

    std::printf("Before staging:\n");
    for (int j=0;j<leaf.gas_count();++j) {
        int cidx = 100 + leaf.type_offset[(int)ParticleType::Gas] + j;
        int gidx = leaf.gas_array_begin + j;
        std::printf("  leaf.gas[%d] -> common[%d] mass=%.1f; gas[%d] density=%.1f\n",
                    j, cidx, core[cidx].mass, gidx, gas[gidx].density);
    }
    std::printf("  count_star = %llu, leaf.staged_star_count = %d\n\n",
                (unsigned long long)container.count_star(),
                (int)leaf.staged_star_count);

    // ---- stage gas #1 for conversion to star -------------------------------
    idx_t new_star_slot = stage_gas_to_star(container, leaf, /*j_within_leaf*/1);

    std::printf("Staged leaf.gas[1] -> new star staging slot %d\n", (int)new_star_slot);
    std::printf("  count_star = %llu, leaf.staged_star_count = %d, leaf.gas_count=%d\n\n",
                (unsigned long long)container.count_star(),
                (int)leaf.staged_star_count,
                (int)leaf.gas_count());

    std::printf("After staging (surviving gas — positional alignment must hold):\n");
    for (int j=0; j<leaf.gas_count(); ++j) {
        int cidx = 100 + leaf.type_offset[(int)ParticleType::Gas] + j;
        int gidx = leaf.gas_array_begin + j;
        std::printf("  leaf.gas[%d] -> common[%d] mass=%.1f; gas[%d] density=%.1f\n",
                    j, cidx, core[cidx].mass, gidx, gas[gidx].density);
    }

    std::printf("\nStaged particle parked at common[103] with type = 0x%02x (high bit set = staged)\n",
                core[103].type);
    std::printf("SPH loops see only leaf.gas_count() = %d particles — no per-particle branch.\n",
                (int)leaf.gas_count());

    return 0;
}
