/**
 * @file example_layout_b_coarse.cpp
 * @brief End-to-end smoke test for layout B (coarse) with the tight MPI bridge.
 */

#include "opg_layout_b_coarse/particle_container.hpp"
#include "opg_layout_b_coarse/reshuffle.hpp"
#include "opg_layout_b_coarse/pack_unpack.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

using namespace opg::common;
using namespace opg::layout_b_coarse;
using namespace opg::common::physics_configs;

constexpr auto CFG = StandardSPH;
using Container = ParticleContainer<CFG>;

static void print_sizes() {
    std::printf("=== Layout B (coarse) sizes under StandardSPH ===\n");
    std::printf("  PCoreB                     = %zu B\n", sizeof(PCoreB));
    std::printf("  PDynB<Cfg>                 = %zu B\n", sizeof(PDynB<CFG>));
    std::printf("  PAuxB<Cfg>                 = %zu B\n", sizeof(PAuxB<CFG>));
    std::printf("  PLinkageB                  = %zu B\n", sizeof(PLinkageB));
    std::printf("  common wire bytes/particle = %zu B\n", common_record_bytes<CFG>());
    std::printf("  GasAllB<Cfg>               = %zu B\n", sizeof(GasAllB<CFG>));
    std::printf("  StarAllB<Cfg>              = %zu B\n", sizeof(StarAllB<CFG>));
    std::printf("  BHAllB<Cfg>                = %zu B\n", sizeof(BHAllB<CFG>));
    std::printf("  tight header               = %zu B\n", packed_header_bytes<CFG>());
    std::printf("  legacy fixed-record bytes  = %zu B\n", legacy_fixed_record_bytes<CFG>());
    std::printf("\n");
}

static void seed(Container& c, std::mt19937_64& rng) {
    const count_t N_gas  = 10;
    const count_t N_dm   = 20;
    const count_t N_star = 5;
    const count_t N_bh   = 2;
    const count_t N_tot  = N_gas + N_dm + N_star + N_bh;

    c.set_count(N_tot);
    c.set_count_gas(N_gas);
    c.set_count_star(N_star);
    c.set_count_bh(N_bh);

    std::uniform_real_distribution<double> upos(0.0, 1.0);
    std::uniform_int_distribution<uint64_t> ukey(0, (1ULL << 63) - 1);

    idx_t j = 0;
    for (count_t i = 0; i < N_gas; ++i, ++j) {
        auto& co = c.core()[j];
        co.pos = {upos(rng), upos(rng), upos(rng)};
        co.mass = 0.01f;
        co.key  = ukey(rng);
        co.type = static_cast<uint8_t>(ParticleType::Gas);
        co.flags = 0;
        co.time_bin = 10;
        co.leaf_idx = -1;

        auto& dy = c.dyn()[j];
        dy.vel = {static_cast<float>(i), static_cast<float>(i + 1), static_cast<float>(i + 2)};
        dy.grav_acc = {0.1 * i, 0.2 * i, 0.3 * i};
        dy.old_acc = 0.01f * i;

        auto& au = c.aux()[j];
        au.id = 1000 + i;
        au.ti_begin = 0;
        au.ti_current = 100;
        au.dt_step = 50;
        au.grav_cost = 1.0f;
        au.num_ngb = 32.0f;
        au.true_ngb = 31;

        auto& g = c.gas_all()[i];
        std::memset(&g, 0, sizeof(g));
        g.core.hsml = 0.05f + 0.001f * i;
        g.core.density = 1.0f + 0.1f * i;
        g.core.pressure = 0.5f * i;
        g.core.sound_speed = 1.2f;

        c.linkage()[j].type_idx = static_cast<uint32_t>(i);
    }

    for (count_t i = 0; i < N_dm; ++i, ++j) {
        auto& co = c.core()[j];
        co.pos = {upos(rng), upos(rng), upos(rng)};
        co.mass = 0.1f;
        co.key  = ukey(rng);
        co.type = static_cast<uint8_t>(ParticleType::DM1);
        co.flags = 0;
        co.time_bin = 8;
        co.leaf_idx = -1;
        c.dyn()[j] = {};
        c.dyn()[j].vel = {static_cast<float>(j), static_cast<float>(j), static_cast<float>(j)};
        c.aux()[j] = {};
        c.aux()[j].id = 2000 + i;
        c.linkage()[j].type_idx = 0;
    }

    for (count_t i = 0; i < N_star; ++i, ++j) {
        auto& co = c.core()[j];
        co.pos = {upos(rng), upos(rng), upos(rng)};
        co.mass = 0.005f;
        co.key  = ukey(rng);
        co.type = static_cast<uint8_t>(ParticleType::Star);
        co.flags = 0;
        co.time_bin = 12;
        co.leaf_idx = -1;
        c.dyn()[j] = {};
        c.aux()[j] = {};
        c.aux()[j].id = 3000 + i;

        auto& s = c.star_all()[i];
        std::memset(&s, 0, sizeof(s));
        s.core.stellar_age = 1e8f + 1e7f * i;
        s.meta.z_age = 1.0f + i;
        s.meta.init_z = 0.02f;

        c.linkage()[j].type_idx = static_cast<uint32_t>(i);
    }

    for (count_t i = 0; i < N_bh; ++i, ++j) {
        auto& co = c.core()[j];
        co.pos = {upos(rng), upos(rng), upos(rng)};
        co.mass = 1e5f;
        co.key  = ukey(rng);
        co.type = static_cast<uint8_t>(ParticleType::BH);
        co.flags = 0;
        co.time_bin = 14;
        co.leaf_idx = -1;
        c.dyn()[j] = {};
        c.aux()[j] = {};
        c.aux()[j].id = 4000 + i;

        auto& b = c.bh_all()[i];
        std::memset(&b, 0, sizeof(b));
        b.core.bh_mass = 1e5f + 1e4f * i;
        b.core.bh_mdot = 1e-3f;

        c.linkage()[j].type_idx = static_cast<uint32_t>(i);
    }
}

static int verify_round_trip(const Container& src, const Container& dst) {
    const count_t N = src.n_part();
    int fails = 0;
    int fails_core = 0;
    int fails_dyn = 0;
    int fails_aux = 0;
    int fails_typed = 0;

    for (count_t i = 0; i < N; ++i) {
        if (std::memcmp(&src.core()[i], &dst.core()[i], sizeof(PCoreB)) != 0) {
            ++fails; ++fails_core;
        }
        if (std::memcmp(&src.dyn()[i], &dst.dyn()[i], sizeof(src.dyn()[i])) != 0) {
            ++fails; ++fails_dyn;
        }
        if (std::memcmp(&src.aux()[i], &dst.aux()[i], sizeof(src.aux()[i])) != 0) {
            ++fails; ++fails_aux;
        }

        const uint8_t t = src.core()[i].type & 0x7Fu;
        const idx_t s_idx = static_cast<idx_t>(src.linkage()[i].type_idx);
        const idx_t d_idx = static_cast<idx_t>(dst.linkage()[i].type_idx);

        switch (static_cast<ParticleType>(t)) {
            case ParticleType::Gas:
                if (std::memcmp(&src.gas_all()[s_idx], &dst.gas_all()[d_idx], sizeof(src.gas_all()[0])) != 0) {
                    ++fails; ++fails_typed;
                }
                break;
            case ParticleType::Star:
                if (std::memcmp(&src.star_all()[s_idx], &dst.star_all()[d_idx], sizeof(src.star_all()[0])) != 0) {
                    ++fails; ++fails_typed;
                }
                break;
            case ParticleType::BH:
                if (std::memcmp(&src.bh_all()[s_idx], &dst.bh_all()[d_idx], sizeof(src.bh_all()[0])) != 0) {
                    ++fails; ++fails_typed;
                }
                break;
            default:
                break;
        }
    }

    if (fails != 0) {
        std::fprintf(stderr, "  fails: core=%d dyn=%d aux=%d typed=%d\n",
                     fails_core, fails_dyn, fails_aux, fails_typed);
    }
    return fails;
}

int main() {
    print_sizes();

    ArenaConfig acfg{};
    MemoryArena arena_src(16 * 1024 * 1024, acfg);
    MemoryArena arena_dst(16 * 1024 * 1024, acfg);

    Capacity cap{.n_max = 64, .n_max_gas = 32, .n_max_star = 16, .n_max_bh = 8};
    Container src(arena_src, cap);
    Container dst(arena_dst, cap);

    std::mt19937_64 rng(42);
    seed(src, rng);

    std::printf("=== Seeded source container ===\n"
                "  n_part = %lu (gas=%lu, dm=%lu, star=%lu, bh=%lu)\n\n",
                static_cast<unsigned long>(src.n_part()),
                static_cast<unsigned long>(src.n_gas()),
                static_cast<unsigned long>(src.n_part() - src.n_gas() - src.n_star() - src.n_bh()),
                static_cast<unsigned long>(src.n_star()),
                static_cast<unsigned long>(src.n_bh()));

    std::printf("=== PH-key reshuffle ===\n");
    std::printf("  common registry size = %zu (expect 4: Core, Dyn, Aux, Linkage)\n",
                src.common_registry().size());

    std::vector<SortHelper> helpers(src.n_part());
    std::vector<idx_t> perm(src.n_part());
    ArenaTempCheckpoint scratch_scope(arena_src);
    const size_t scratch_bytes = src.common_registry().max_elem_size() * src.n_part();
    auto scratch = arena_src.allocate_temp_bytes(scratch_bytes, arena_src.config().alignment, "reshuffle scratch");

    std::vector<opg::common::pid_t> ids_before(src.n_part());
    for (count_t i = 0; i < src.n_part(); ++i) ids_before[i] = src.aux()[i].id;

    do_reshuffle_by_key(src, helpers.data(), perm.data(), scratch.ptr);

    bool ok_permutation = true;
    for (count_t i = 0; i < src.n_part(); ++i) {
        if (src.aux()[i].id != ids_before[perm[i]]) {
            ok_permutation = false;
            break;
        }
    }

    bool ok_sorted = true;
    for (count_t i = 1; i < src.n_part(); ++i) {
        if (src.core()[i].key < src.core()[i - 1].key) {
            ok_sorted = false;
            break;
        }
    }

    std::printf("  permutation consistent across common arrays: %s\n", ok_permutation ? "OK" : "FAIL");
    std::printf("  keys are non-decreasing after reshuffle: %s\n\n", ok_sorted ? "OK" : "FAIL");

    std::printf("=== Tight MPI-byte bridge round-trip ===\n");

    const count_t N = src.n_part();
    const PackedExchangeCounts counts = count_range_exchange(src, 0, N);
    const size_t tight_bytes = exchange_bytes<CFG>(counts);
    const size_t legacy_bytes = static_cast<size_t>(N) * legacy_fixed_record_bytes<CFG>();

    ArenaTempCheckpoint wire_scope(arena_src);
    auto wire_block = arena_src.allocate_temp_bytes(tight_bytes, arena_src.config().alignment, "tight wire buffer");
    auto wire = make_exchange_view<CFG>(wire_block.ptr, wire_block.bytes);

    const bool packed_ok = pack_range_exchange(src, 0, N, wire);
    const PackedExchangeHeader header = load_exchange_header(wire);
    const bool header_ok = validate_exchange_header<CFG>(wire, header);
    const bool unpack_ok = unpack_exchange(dst, 0, wire.data, wire.bytes);

    std::printf("  header counts            = total=%lu gas=%lu star=%lu bh=%lu\n",
                static_cast<unsigned long>(header.n_particles),
                static_cast<unsigned long>(header.n_gas),
                static_cast<unsigned long>(header.n_star),
                static_cast<unsigned long>(header.n_bh));
    std::printf("  tight wire bytes         = %zu\n", tight_bytes);
    std::printf("  legacy wire bytes        = %zu\n", legacy_bytes);
    if (legacy_bytes > 0) {
        const double saved = 100.0 * static_cast<double>(legacy_bytes - tight_bytes) / static_cast<double>(legacy_bytes);
        std::printf("  wire-volume reduction    = %.2f %%\n", saved);
    }
    std::printf("  unpacked into destination: %s\n\n", (packed_ok && header_ok && unpack_ok) ? "OK" : "FAIL");

    if (!packed_ok || !header_ok || !unpack_ok) return 1;

    const int fails = verify_round_trip(src, dst);
    std::printf("  destination counts       = gas=%lu star=%lu bh=%lu\n",
                static_cast<unsigned long>(dst.n_gas()),
                static_cast<unsigned long>(dst.n_star()),
                static_cast<unsigned long>(dst.n_bh()));

    if (fails == 0) {
        std::printf("  round-trip fidelity      = OK (bit-for-bit match across %lu particles)\n\n",
                    static_cast<unsigned long>(N));
    } else {
        std::fprintf(stderr, "  round-trip fidelity: FAIL (%d mismatches)\n", fails);
        return 1;
    }

    std::printf("All checks passed.\n");
    return 0;
}
