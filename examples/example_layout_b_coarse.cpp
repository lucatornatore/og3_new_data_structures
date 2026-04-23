/**
 * @file example_layout_b_coarse.cpp
 * @brief End-to-end smoke test for layout B (coarse).
 *
 * Exercises:
 *   1. Container construction and struct-size report.
 *   2. Seeding: mixed gas/DM/star/BH particles with synthetic PH keys.
 *   3. PH-key reshuffle — demonstrates that under B only the 4 common arrays
 *      participate in the permutation (registry size = 4).
 *   4. Pack/unpack round-trip:
 *         a. Pack particles 0..N-1 of a SOURCE container into a linear
 *            ParticleRecord[] buffer (the distribute.c `data_type*` analogue).
 *         b. Unpack the buffer into a DESTINATION container and verify every
 *            field round-trips bit-for-bit.
 *   5. Print the sizes that the distribute.c harness needs:
 *         sizeof(ParticleRecord<Cfg>)    — the DATA_SIZE to compile with.
 *         sizeof(PCoreB / PDynB / PAuxB) — the common payload breakdown.
 */

#include "opg_layout_b_coarse/particle_container.hpp"
#include "opg_layout_b_coarse/reshuffle.hpp"
#include "opg_layout_b_coarse/pack_unpack.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <random>
#include <vector>

using namespace opg::common;
using namespace opg::layout_b_coarse;
using namespace opg::common::physics_configs;

// Pick a representative config. Adjust to experiment with sizes.
constexpr auto CFG = StandardSPH;

using Container = ParticleContainer<CFG>;
using Record    = ParticleRecord<CFG>;

static void print_sizes() {
    std::printf("=== Layout B (coarse) sizes under %s ===\n",
                "StandardSPH (SPH + metals + BHs, no MHD)");
    std::printf("  PCoreB                 = %zu B\n", sizeof(PCoreB));
    std::printf("  PDynB<Cfg>             = %zu B\n", sizeof(PDynB<CFG>));
    std::printf("  PAuxB<Cfg>             = %zu B\n", sizeof(PAuxB<CFG>));
    std::printf("  PLinkageB              = %zu B\n", sizeof(PLinkageB));
    std::printf("  GasAllB<Cfg>           = %zu B\n", sizeof(GasAllB<CFG>));
    std::printf("  StarAllB<Cfg>          = %zu B\n", sizeof(StarAllB<CFG>));
    std::printf("  BHAllB<Cfg>            = %zu B\n", sizeof(BHAllB<CFG>));
    std::printf("  ParticleRecord<Cfg>    = %zu B   (the exchange unit)\n",
                sizeof(Record));
    std::printf("  typed-payload bytes    = %zu B\n",
                record_typed_payload_bytes<CFG>());
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

    // Layout of the common arrays: first gas, then DM (type=1), then star, then BH.
    // (Not strictly required under B — but mimics initial-file convention.)
    idx_t j = 0;

    // gas
    for (count_t i = 0; i < N_gas; ++i, ++j) {
        auto& co = c.core()[j];
        co.pos = {upos(rng), upos(rng), upos(rng)};
        co.mass = 0.01f;
        co.key  = ukey(rng);
        co.type = (uint8_t)ParticleType::Gas;
        co.flags = 0;
        co.time_bin = 10;
        co.leaf_idx = -1;

        // seed dyn/aux with distinguishable patterns
        auto& dy = c.dyn()[j];
        dy.vel = {(float)i, (float)(i+1), (float)(i+2)};
        dy.grav_acc = {0.1*i, 0.2*i, 0.3*i};
        dy.old_acc  = 0.01f * i;

        auto& au = c.aux()[j];
        au.id = 1000 + i;
        au.ti_begin = 0;
        au.ti_current = 100;
        au.dt_step = 50;
        au.grav_cost = 1.0f;
        au.num_ngb = 32.0f;
        au.true_ngb = 31;

        // typed payload: write GasCore fields
        auto& g = c.gas_all()[i];
        std::memset(&g, 0, sizeof(g));
        g.core.hsml = 0.05f + 0.001f * i;
        g.core.density = 1.0f + 0.1f * i;
        g.core.pressure = 0.5f * i;
        g.core.sound_speed = 1.2f;

        // linkage: j-th common slot → i-th gas slot
        c.linkage()[j].type_idx = (uint32_t)i;
    }

    // DM (type=1)
    for (count_t i = 0; i < N_dm; ++i, ++j) {
        auto& co = c.core()[j];
        co.pos = {upos(rng), upos(rng), upos(rng)};
        co.mass = 0.1f;
        co.key  = ukey(rng);
        co.type = 1;
        co.flags = 0;
        co.time_bin = 8;
        co.leaf_idx = -1;
        c.dyn()[j] = {};
        c.dyn()[j].vel = {(float)j, (float)j, (float)j};
        c.aux()[j] = {};
        c.aux()[j].id = 2000 + i;
        c.linkage()[j].type_idx = 0;  // unused for DM
    }

    // star
    for (count_t i = 0; i < N_star; ++i, ++j) {
        auto& co = c.core()[j];
        co.pos = {upos(rng), upos(rng), upos(rng)};
        co.mass = 0.005f;
        co.key  = ukey(rng);
        co.type = (uint8_t)ParticleType::Star;
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

        c.linkage()[j].type_idx = (uint32_t)i;
    }

    // BH
    for (count_t i = 0; i < N_bh; ++i, ++j) {
        auto& co = c.core()[j];
        co.pos = {upos(rng), upos(rng), upos(rng)};
        co.mass = 1e5f;
        co.key  = ukey(rng);
        co.type = (uint8_t)ParticleType::BH;
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

        c.linkage()[j].type_idx = (uint32_t)i;
    }
}

int main() {
    print_sizes();

    // --- two containers that will round-trip particles between them -----------
    ArenaConfig acfg{};
    MemoryArena arena_src(16 * 1024 * 1024, acfg);
    MemoryArena arena_dst(16 * 1024 * 1024, acfg);

    Capacity cap{
        .n_max      = 64,
        .n_max_gas  = 32,
        .n_max_star = 16,
        .n_max_bh   = 8,
    };

    Container src(arena_src, cap);
    Container dst(arena_dst, cap);

    std::mt19937_64 rng(42);
    seed(src, rng);

    std::printf("=== Seeded source container ===\n"
                "  n_part = %lu (gas=%lu, dm=%lu, star=%lu, bh=%lu)\n\n",
                (unsigned long)src.n_part(),
                (unsigned long)src.n_gas(),
                (unsigned long)(src.n_part() - src.n_gas() - src.n_star() - src.n_bh()),
                (unsigned long)src.n_star(),
                (unsigned long)src.n_bh());

    // --- 1. test reshuffle ---------------------------------------------------
    std::printf("=== PH-key reshuffle ===\n");
    std::printf("  common registry size = %zu (expect 4: Core, Dyn, Aux, Linkage)\n",
                src.common_registry().size());

    std::vector<SortHelper> helpers(src.n_part());
    std::vector<idx_t>      perm(src.n_part());
    // scratch: need at least max_elem_size * n bytes
    const size_t scratch_bytes = 64 * src.n_part();
    std::vector<std::byte>  scratch(scratch_bytes);

    // Snapshot of common state BEFORE the shuffle (to verify values re-appear
    // in their new positions after permutation).
    std::vector<opg::common::pid_t> ids_before(src.n_part());
    for (count_t i = 0; i < src.n_part(); ++i)
        ids_before[i] = src.aux()[i].id;

    do_reshuffle_by_key(src, helpers.data(), perm.data(), scratch.data());

    // Verify: src.aux()[i].id should now equal ids_before[perm[i]], for every i.
    bool ok_permutation = true;
    for (count_t i = 0; i < src.n_part(); ++i) {
        if (src.aux()[i].id != ids_before[perm[i]]) {
            ok_permutation = false;
            break;
        }
    }
    std::printf("  permutation consistent across common arrays: %s\n",
                ok_permutation ? "OK" : "FAIL");

    // Verify keys are now sorted.
    bool ok_sorted = true;
    for (count_t i = 1; i < src.n_part(); ++i) {
        if (src.core()[i].key < src.core()[i-1].key) {
            ok_sorted = false;
            break;
        }
    }
    std::printf("  keys are non-decreasing after reshuffle: %s\n\n",
                ok_sorted ? "OK" : "FAIL");

    // --- 2. test pack/unpack round-trip --------------------------------------
    std::printf("=== Pack/unpack round-trip ===\n");

    const count_t N = src.n_part();
    std::vector<Record> wire(N);  // the "MPI send/recv buffer"

    pack_range(src, 0, N, wire.data());

    const count_t unpacked = unpack_range(dst, 0, N, wire.data());
    std::printf("  packed %lu records of %zu B each (= %zu kB total wire)\n",
                (unsigned long)N, sizeof(Record), (N * sizeof(Record)) / 1024);
    std::printf("  unpacked %lu records into destination\n",
                (unsigned long)unpacked);

    if (unpacked != N) {
        std::fprintf(stderr, "  FAIL: partial unpack\n");
        return 1;
    }

    // --- 3. verify round-trip fidelity --------------------------------------
    int fails = 0;
    int fails_core = 0, fails_dyn = 0, fails_aux = 0, fails_typed = 0;

    for (count_t i = 0; i < N; ++i) {
        const auto& src_core = src.core()[i];
        const auto& dst_core = dst.core()[i];
        if (std::memcmp(&src_core, &dst_core, sizeof(PCoreB)) != 0) {
            ++fails; ++fails_core;
            if (fails_core <= 2) {
                std::fprintf(stderr,
                    "  core mismatch at i=%lu: src.type=%u dst.type=%u src.key=%lx dst.key=%lx\n",
                    (unsigned long)i,
                    (unsigned)src_core.type, (unsigned)dst_core.type,
                    (unsigned long)src_core.key, (unsigned long)dst_core.key);
            }
        }
        if (std::memcmp(&src.dyn()[i], &dst.dyn()[i], sizeof(src.dyn()[i])) != 0) {
            ++fails; ++fails_dyn;
        }
        if (std::memcmp(&src.aux()[i], &dst.aux()[i], sizeof(src.aux()[i])) != 0) {
            ++fails; ++fails_aux;
            if (fails_aux <= 2) {
                std::fprintf(stderr,
                    "  aux mismatch at i=%lu: src.id=%lu dst.id=%lu\n",
                    (unsigned long)i,
                    (unsigned long)src.aux()[i].id,
                    (unsigned long)dst.aux()[i].id);
            }
        }

        // Type-specific: the destination has reassigned linkage, so we look up
        // the typed entry on each side via its own linkage.
        const uint8_t t = src_core.type & 0x7Fu;
        const idx_t   s_idx = (idx_t)src.linkage()[i].type_idx;
        const idx_t   d_idx = (idx_t)dst.linkage()[i].type_idx;

        switch (static_cast<ParticleType>(t)) {
            case ParticleType::Gas:
                if (std::memcmp(&src.gas_all()[s_idx], &dst.gas_all()[d_idx],
                                sizeof(src.gas_all()[0])) != 0) {
                    ++fails; ++fails_typed;
                    if (fails_typed <= 2) {
                        std::fprintf(stderr,
                            "  gas-typed mismatch at i=%lu (s_idx=%d, d_idx=%d)\n",
                            (unsigned long)i, (int)s_idx, (int)d_idx);
                    }
                }
                break;
            case ParticleType::Star:
                if (std::memcmp(&src.star_all()[s_idx], &dst.star_all()[d_idx],
                                sizeof(src.star_all()[0])) != 0) {
                    ++fails; ++fails_typed;
                }
                break;
            case ParticleType::BH:
                if (std::memcmp(&src.bh_all()[s_idx], &dst.bh_all()[d_idx],
                                sizeof(src.bh_all()[0])) != 0) {
                    ++fails; ++fails_typed;
                }
                break;
            default:
                break;  // DM: no typed payload to compare
        }
    }
    std::fprintf(stderr, "  fails: core=%d dyn=%d aux=%d typed=%d\n",
                 fails_core, fails_dyn, fails_aux, fails_typed);

    // Also verify type-specific counts came out right on the destination.
    std::printf("  dst counts: gas=%lu, star=%lu, bh=%lu (expect %lu, %lu, %lu)\n",
                (unsigned long)dst.n_gas(),
                (unsigned long)dst.n_star(),
                (unsigned long)dst.n_bh(),
                (unsigned long)src.n_gas(),
                (unsigned long)src.n_star(),
                (unsigned long)src.n_bh());

    if (fails == 0) {
        std::printf("  round-trip fidelity: OK (bit-for-bit match across %lu particles)\n\n",
                    (unsigned long)N);
    } else {
        std::fprintf(stderr, "  round-trip fidelity: FAIL (%d mismatches)\n", fails);
        return 1;
    }

    std::printf("All checks passed.\n");
    return 0;
}
