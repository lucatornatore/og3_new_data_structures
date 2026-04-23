/**
 * @file print_record_size.cpp
 * @brief Utility: print sizeof(ParticleRecord<Cfg>) for distribute.c wiring.
 *
 * The distribute.c exchange benchmark takes a -DDATA_SIZE=N compile flag
 * where N is the number of ints per cell. To run distribute.c with a
 * "payload size" equal to one layout-B particle record, build this utility
 * and pipe its output into the distribute.c build:
 *
 *     g++ -std=c++20 -I .../opg_common/include -I .../opg_layout_b_coarse/include \
 *         print_record_size.cpp -o print_record_size
 *
 *     N_INTS=$(./print_record_size StandardSPH ints)
 *     mpicc -DDATA_SIZE=$N_INTS distribute.c -o distribute_sph
 *
 * Or, under -DUSE_SoA, the bytes figure may be useful too:
 *
 *     N_BYTES=$(./print_record_size StandardSPH bytes)
 *
 * Supported presets: GravityOnly, StandardSPH, FullMHD, MFMHydro.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "opg_layout_b_coarse/pack_unpack.hpp"

using namespace opg::common;
using namespace opg::layout_b_coarse;
using namespace opg::common::physics_configs;

template<PhysicsConfig Cfg>
static void emit(const char* mode) {
    const size_t bytes = record_bytes<Cfg>();
    if (std::strcmp(mode, "bytes") == 0) {
        std::printf("%zu\n", bytes);
    } else if (std::strcmp(mode, "ints") == 0) {
        // distribute.c uses `int data[DATA_SIZE]` -> DATA_SIZE is ints-per-cell.
        std::printf("%zu\n", bytes / sizeof(int));
    } else {
        std::fprintf(stderr, "ERR: unknown mode '%s'; use 'bytes' or 'ints'\n", mode);
        std::exit(2);
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
            "usage: %s <preset> [bytes|ints]\n"
            "  <preset>   : GravityOnly | StandardSPH | FullMHD | MFMHydro\n"
            "  [mode]     : 'bytes' (default) or 'ints' (for distribute.c -DDATA_SIZE)\n",
            argv[0]);
        return 1;
    }
    const std::string preset = argv[1];
    const char* mode = (argc >= 3) ? argv[2] : "bytes";

    if      (preset == "GravityOnly")  emit<GravityOnly>(mode);
    else if (preset == "StandardSPH")  emit<StandardSPH>(mode);
    else if (preset == "FullMHD")      emit<FullMHD>(mode);
    else if (preset == "MFMHydro")     emit<MFMHydro>(mode);
    else {
        std::fprintf(stderr, "ERR: unknown preset '%s'\n", preset.c_str());
        return 1;
    }
    return 0;
}
