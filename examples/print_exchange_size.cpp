/**
 * @file print_exchange_size.cpp
 * @brief Utility: print exact coarse-B packed-message sizes for a particle mix.
 *
 * Usage:
 *   ./print_exchange_size <preset> <n_total> <n_gas> <n_star> <n_bh>
 */

#include <cstdio>
#include <cstdlib>
#include <string>

#include "opg_layout_b_coarse/pack_unpack.hpp"

using namespace opg::common;
using namespace opg::layout_b_coarse;
using namespace opg::common::physics_configs;

template<PhysicsConfig Cfg>
static int emit(const char* preset_name,
                const char* total_s,
                const char* gas_s,
                const char* star_s,
                const char* bh_s) {
    const auto n_total = static_cast<count_t>(std::strtoull(total_s, nullptr, 10));
    const auto n_gas   = static_cast<count_t>(std::strtoull(gas_s,   nullptr, 10));
    const auto n_star  = static_cast<count_t>(std::strtoull(star_s,  nullptr, 10));
    const auto n_bh    = static_cast<count_t>(std::strtoull(bh_s,    nullptr, 10));

    if (n_gas + n_star + n_bh > n_total) {
        std::fprintf(stderr, "ERR: typed counts exceed n_total\n");
        return 2;
    }

    PackedExchangeCounts counts{};
    counts.n_particles = n_total;
    counts.n_gas = n_gas;
    counts.n_star = n_star;
    counts.n_bh = n_bh;

    const auto hdr = make_exchange_header<Cfg>(counts);
    const size_t legacy_bytes = static_cast<size_t>(n_total) * legacy_fixed_record_bytes<Cfg>();
    const double saving = legacy_bytes > 0
        ? 100.0 * (1.0 - static_cast<double>(hdr.total_bytes) / static_cast<double>(legacy_bytes))
        : 0.0;

    std::printf("preset              : %s\n", preset_name);
    std::printf("header_bytes        : %zu\n", sizeof(PackedExchangeHeader));
    std::printf("common_record_bytes : %zu\n", common_record_bytes<Cfg>());
    std::printf("gas_record_bytes    : %zu\n", gas_record_bytes<Cfg>());
    std::printf("star_record_bytes   : %zu\n", star_record_bytes<Cfg>());
    std::printf("bh_record_bytes     : %zu\n", bh_record_bytes<Cfg>());
    std::printf("message_bytes       : %zu\n", static_cast<size_t>(hdr.total_bytes));
    std::printf("  common_block      : %zu\n", static_cast<size_t>(hdr.common_bytes));
    std::printf("  gas_block         : %zu\n", static_cast<size_t>(hdr.gas_bytes));
    std::printf("  star_block        : %zu\n", static_cast<size_t>(hdr.star_bytes));
    std::printf("  bh_block          : %zu\n", static_cast<size_t>(hdr.bh_bytes));
    std::printf("legacy_fixed_bytes  : %zu\n", legacy_bytes);
    std::printf("wire_saving_percent : %.3f\n", saving);
    return 0;
}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::fprintf(stderr,
            "usage: %s <preset> <n_total> <n_gas> <n_star> <n_bh>\n"
            "  <preset> : GravityOnly | StandardSPH | FullMHD | MFMHydro\n",
            argv[0]);
        return 1;
    }

    const std::string preset = argv[1];
    if      (preset == "GravityOnly") return emit<GravityOnly>(argv[1], argv[2], argv[3], argv[4], argv[5]);
    else if (preset == "StandardSPH") return emit<StandardSPH>(argv[1], argv[2], argv[3], argv[4], argv[5]);
    else if (preset == "FullMHD")     return emit<FullMHD>(argv[1], argv[2], argv[3], argv[4], argv[5]);
    else if (preset == "MFMHydro")    return emit<MFMHydro>(argv[1], argv[2], argv[3], argv[4], argv[5]);

    std::fprintf(stderr, "ERR: unknown preset '%s'\n", preset.c_str());
    return 1;
}
