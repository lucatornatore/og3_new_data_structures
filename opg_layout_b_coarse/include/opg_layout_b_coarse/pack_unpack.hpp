/**
 * @file pack_unpack.hpp  (layout B, coarse)
 * @brief MPI-ready pack/unpack API for layout B's particle exchange test.
 *
 * ---------------------------------------------------------------------------
 * RELATIONSHIP TO distribute.c
 * ---------------------------------------------------------------------------
 * The existing distribute.c test harness exchanges fixed-size `data_type`
 * records across MPI tasks according to a 3D region map. To reuse that
 * harness for particle exchange, we need a single type that plays the role
 * of `data_type`: a fixed-size POD describing one entire particle.
 *
 * `ParticleRecord<Cfg>` is that type. It is a flat POD containing:
 *
 *   +----------------------------------+
 *   |  PCoreB     (64 B)               | <- always
 *   +----------------------------------+
 *   |  PDynB<Cfg> (64 B)               | <- always
 *   +----------------------------------+
 *   |  PAuxB<Cfg> (64 B)               | <- always
 *   +----------------------------------+
 *   |  typed payload: max of           |
 *   |    sizeof(GasAllB<Cfg>),         |
 *   |    sizeof(StarAllB<Cfg>),        |
 *   |    sizeof(BHAllB<Cfg>)           |
 *   |  (The actual type is recorded in |
 *   |   record.core.type.)             |
 *   +----------------------------------+
 *
 * PLinkageB is deliberately NOT serialised. On the sender it indexes into
 * the sender's type-specific arrays; on the receiver those indices are
 * meaningless. The receiver reconstructs linkage from the record's type
 * field as it places the typed payload into its own arrays.
 *
 * ---------------------------------------------------------------------------
 * TWO PACKING MODES
 * ---------------------------------------------------------------------------
 * MODE 1 (this header): one fixed-size record per particle.
 *   - Pro: trivial mapping to distribute.c's data_type[]. One contiguous
 *          send buffer per target task.
 *   - Con: DM particles waste the typed-payload bytes (padded, not sent
 *          content meaningful to the receiver).
 *
 * MODE 2 (future): grouped by type.
 *   - Sender partitions the outgoing particles by type, sends common data
 *     for all particles, then separate gas/star/bh payloads only for the
 *     matching particles.
 *   - Optimal bandwidth; more bookkeeping. Will be written when we move
 *     beyond the first test.
 *
 * For the first test (reproducing distribute.c's AoS-vs-SoA measurement
 * on particle data), MODE 1 is the direct analogue. That is what this
 * file provides.
 *
 * ---------------------------------------------------------------------------
 * THE DATA_TYPE BRIDGE TO distribute.c
 * ---------------------------------------------------------------------------
 * In distribute.c, `data_type` is
 *     typedef struct { int data[DATA_SIZE]; } data_type;
 * where DATA_SIZE comes from a compile-time `-DDATA_SIZE=...` flag.
 *
 * The C++ ParticleRecord<Cfg> has a fixed, C++-derived size. If the test
 * harness requires DATA_SIZE to match sizeof(ParticleRecord), the user
 * builds distribute.c with
 *     -DDATA_SIZE=$(printf '%d' $(./print_record_size StandardSPH))
 * where `print_record_size` is a tiny program that prints
 * `sizeof(ParticleRecord<Cfg>) / sizeof(int)`. See examples/.
 */

#ifndef OPG_LAYOUT_B_COARSE_PACK_UNPACK_HPP
#define OPG_LAYOUT_B_COARSE_PACK_UNPACK_HPP

#include <cstddef>
#include <cstring>
#include <algorithm>

#include "particle_container.hpp"

namespace opg::layout_b_coarse {

using namespace opg::common;

// =============================================================================
// Helpers for compile-time payload size
// =============================================================================

namespace detail {
    template<PhysicsConfig Cfg>
    consteval size_t typed_payload_bytes() noexcept {
        size_t s = sizeof(GasAllB<Cfg>);
        if constexpr (HasStellarEvolution<Cfg>)
            s = s > sizeof(StarAllB<Cfg>) ? s : sizeof(StarAllB<Cfg>);
        if constexpr (HasBlackHoles<Cfg>)
            s = s > sizeof(BHAllB<Cfg>)   ? s : sizeof(BHAllB<Cfg>);
        return s;
    }
} // namespace detail

// =============================================================================
// ParticleRecord<Cfg> — the exchange unit
// =============================================================================
// Trivially-copyable, standard-layout, fixed-size. Suitable for raw MPI_BYTE
// transfer. The typed payload is `alignas(64)`-padded bytes; the actual
// interpretation is determined at unpack time from record.core.type.

template<PhysicsConfig Cfg>
struct alignas(64) ParticleRecord {
    PCoreB       core;
    PDynB<Cfg>   dyn;
    PAuxB<Cfg>   aux;
    alignas(64) std::byte typed_payload[detail::typed_payload_bytes<Cfg>()];
};

// Validation: one record must be trivially copyable for MPI_BYTE, and
// standard-layout so its binary representation is portable across TUs.
namespace _validate_record {
    using namespace opg::common::physics_configs;
    static_assert(std::is_trivially_copyable_v<ParticleRecord<StandardSPH>>);
    static_assert(std::is_standard_layout_v<ParticleRecord<StandardSPH>>);
    static_assert(std::is_trivially_copyable_v<ParticleRecord<FullMHD>>);
    static_assert(std::is_standard_layout_v<ParticleRecord<FullMHD>>);
}

// =============================================================================
// Size query helpers — for the distribute.c bridge and for MPI buffer sizing
// =============================================================================

template<PhysicsConfig Cfg>
consteval size_t record_bytes() noexcept {
    return sizeof(ParticleRecord<Cfg>);
}

template<PhysicsConfig Cfg>
consteval size_t record_typed_payload_bytes() noexcept {
    return detail::typed_payload_bytes<Cfg>();
}

// =============================================================================
// pack_record — serialise one particle into a ParticleRecord
// =============================================================================
// Reads from the container's arrays at common index j. Looks up the
// particle's type-specific slot through linkage_[j].type_idx and copies the
// appropriate GasAllB / StarAllB / BHAllB entry into record.typed_payload.
//
// For DM particles the typed payload is left zeroed (the receiver will
// ignore it anyway based on record.core.type).

template<PhysicsConfig Cfg>
void pack_record(const ParticleContainer<Cfg>& container,
                 idx_t j,
                 ParticleRecord<Cfg>& record) noexcept
{
    // Use explicit memcpy (not struct assignment) to guarantee that all
    // bytes — including alignas-driven trailing padding — flow from source
    // arrays through the wire buffer. Struct assignment on trivially-copyable
    // aggregates is permitted by the standard to skip padding; explicit
    // memcpy is the byte-for-byte semantics MPI transfer requires.
    std::memcpy(&record.core, &container.core()[j], sizeof(PCoreB));
    std::memcpy(&record.dyn,  &container.dyn() [j], sizeof(record.dyn));
    std::memcpy(&record.aux,  &container.aux() [j], sizeof(record.aux));

    // Zero the typed-payload region first so the record is deterministic
    // regardless of particle type (DM has no typed content; its bytes must
    // still be reproducible for any downstream checksum / MPI CRC).
    std::memset(record.typed_payload, 0, sizeof(record.typed_payload));

    const uint8_t bare_type = static_cast<uint8_t>(record.core.type) & 0x7Fu;
    const idx_t   type_idx  = static_cast<idx_t>(
                                  container.linkage()[j].type_idx);

    switch (static_cast<ParticleType>(bare_type)) {
        case ParticleType::Gas: {
            const auto& src = container.gas_all()[type_idx];
            std::memcpy(record.typed_payload, &src, sizeof(src));
            break;
        }
        case ParticleType::Star: {
            if constexpr (HasStellarEvolution<Cfg>) {
                const auto& src = container.star_all()[type_idx];
                std::memcpy(record.typed_payload, &src, sizeof(src));
            }
            break;
        }
        case ParticleType::BH: {
            if constexpr (HasBlackHoles<Cfg>) {
                const auto& src = container.bh_all()[type_idx];
                std::memcpy(record.typed_payload, &src, sizeof(src));
            }
            break;
        }
        default:
            // DM (types 1..3) and any unknown type: no typed payload.
            break;
    }
}

// =============================================================================
// unpack_record — materialise a received record into the container
// =============================================================================
// The receiver supplies the common index j_dst at which the new particle
// will live (typically the next free slot) and the container updates its
// arrays. For typed particles a new slot is claimed in the matching
// type-specific array, the payload is copied there, and linkage is written.
//
// Returns true on success, false if a type-specific array is out of
// capacity or the type is unsupported by the current Cfg.

template<PhysicsConfig Cfg>
bool unpack_record(ParticleContainer<Cfg>& container,
                   idx_t j_dst,
                   const ParticleRecord<Cfg>& record) noexcept
{
    // Byte-for-byte copy (see pack_record comment).
    std::memcpy(&container.core()[j_dst], &record.core, sizeof(PCoreB));
    std::memcpy(&container.dyn() [j_dst], &record.dyn,  sizeof(record.dyn));
    std::memcpy(&container.aux() [j_dst], &record.aux,  sizeof(record.aux));

    const uint8_t bare_type = static_cast<uint8_t>(record.core.type) & 0x7Fu;

    switch (static_cast<ParticleType>(bare_type)) {
        case ParticleType::Gas: {
            const count_t n = container.n_gas();
            if (n >= container.capacity().n_max_gas) return false;
            std::memcpy(&container.gas_all()[n],
                        record.typed_payload,
                        sizeof(container.gas_all()[0]));
            container.linkage()[j_dst].type_idx = static_cast<uint32_t>(n);
            container.set_count_gas(n + 1);
            break;
        }
        case ParticleType::Star: {
            if constexpr (HasStellarEvolution<Cfg>) {
                const count_t n = container.n_star();
                if (n >= container.capacity().n_max_star) return false;
                std::memcpy(&container.star_all()[n],
                            record.typed_payload,
                            sizeof(container.star_all()[0]));
                container.linkage()[j_dst].type_idx = static_cast<uint32_t>(n);
                container.set_count_star(n + 1);
            } else {
                return false;  // received a star but this build has no stars
            }
            break;
        }
        case ParticleType::BH: {
            if constexpr (HasBlackHoles<Cfg>) {
                const count_t n = container.n_bh();
                if (n >= container.capacity().n_max_bh) return false;
                std::memcpy(&container.bh_all()[n],
                            record.typed_payload,
                            sizeof(container.bh_all()[0]));
                container.linkage()[j_dst].type_idx = static_cast<uint32_t>(n);
                container.set_count_bh(n + 1);
            } else {
                return false;
            }
            break;
        }
        default:
            // DM: no typed payload, linkage is unused.
            container.linkage()[j_dst].type_idx = 0;
            break;
    }
    return true;
}

// =============================================================================
// Bulk helpers — pack/unpack contiguous ranges from/to a record buffer
// =============================================================================

template<PhysicsConfig Cfg>
void pack_range(const ParticleContainer<Cfg>& container,
                idx_t start,
                count_t count,
                ParticleRecord<Cfg>* out) noexcept
{
    for (count_t i = 0; i < count; ++i) {
        pack_record(container, start + static_cast<idx_t>(i), out[i]);
    }
}

// Pack particles at non-contiguous indices (the common case when the sender
// has already binned particles by target task but the bins refer to original
// indices that are not contiguous in the array).
template<PhysicsConfig Cfg>
void pack_scatter(const ParticleContainer<Cfg>& container,
                  const idx_t* indices,
                  count_t count,
                  ParticleRecord<Cfg>* out) noexcept
{
    for (count_t i = 0; i < count; ++i) {
        pack_record(container, indices[i], out[i]);
    }
}

// Unpack a contiguous range of records into the container starting at j_dst.
// Returns the number of records successfully unpacked (may be less than
// `count` if a type-specific array fills up mid-batch).
template<PhysicsConfig Cfg>
count_t unpack_range(ParticleContainer<Cfg>& container,
                     idx_t j_dst,
                     count_t count,
                     const ParticleRecord<Cfg>* in) noexcept
{
    for (count_t i = 0; i < count; ++i) {
        if (!unpack_record(container,
                           j_dst + static_cast<idx_t>(i), in[i])) {
            container.set_count(j_dst + i);
            return i;
        }
    }
    container.set_count(j_dst + count);
    return count;
}

} // namespace opg::layout_b_coarse

#endif // OPG_LAYOUT_B_COARSE_PACK_UNPACK_HPP
