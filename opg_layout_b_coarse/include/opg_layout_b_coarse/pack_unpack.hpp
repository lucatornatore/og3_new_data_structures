/**
 * @file pack_unpack.hpp
 * @brief Tight MPI bridge for layout B-coarse.
 *
 * Resident storage stays split as:
 *   - common arrays: PCoreB[], PDynB<Cfg>[], PAuxB<Cfg>[], PLinkageB[]
 *   - typed arrays : GasAllB<Cfg>[], StarAllB<Cfg>[], BHAllB<Cfg>[]
 *
 * The MPI bridge packs one contiguous MPI_BYTE message per neighbour:
 *
 *   [ PackedExchangeHeader ]
 *   [ common block = core,dyn,aux for every outgoing particle ]
 *   [ gas block    = GasAllB entries only for outgoing gas ]
 *   [ star block   = StarAllB entries only for outgoing stars ]
 *   [ bh block     = BHAllB entries only for outgoing BHs ]
 *
 * DM particles carry only common data on the wire. Linkage is rebuilt on the
 * receiver from the type tags stored in PCoreB.
 */

#ifndef OPG_LAYOUT_B_COARSE_PACK_UNPACK_HPP
#define OPG_LAYOUT_B_COARSE_PACK_UNPACK_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>

#include "particle_container.hpp"

namespace opg::layout_b_coarse {

using namespace opg::common;

namespace detail {

inline constexpr uint32_t packed_exchange_magic = 0x4F504742u; // 'OPGB'
inline constexpr uint16_t packed_exchange_version = 2u;

constexpr uint8_t bare_type(uint8_t t) noexcept {
    return static_cast<uint8_t>(t & 0x7Fu);
}

template<typename T>
inline T load_trivial(const std::byte* src) noexcept {
    T value{};
    std::memcpy(&value, src, sizeof(T));
    return value;
}

} // namespace detail

struct PackedExchangeCounts {
    count_t n_particles = 0;
    count_t n_gas = 0;
    count_t n_star = 0;
    count_t n_bh = 0;
};

struct PackedExchangeHeader {
    uint32_t magic = detail::packed_exchange_magic;
    uint16_t version = detail::packed_exchange_version;
    uint16_t reserved = 0;
    count_t n_particles = 0;
    count_t n_gas = 0;
    count_t n_star = 0;
    count_t n_bh = 0;
    uint64_t common_bytes = 0;
    uint64_t gas_bytes = 0;
    uint64_t star_bytes = 0;
    uint64_t bh_bytes = 0;
    uint64_t total_bytes = 0;
};

static_assert(std::is_trivially_copyable_v<PackedExchangeHeader>);
static_assert(std::is_standard_layout_v<PackedExchangeHeader>);

template<PhysicsConfig Cfg>
consteval size_t packed_header_bytes() noexcept {
    return sizeof(PackedExchangeHeader);
}

template<PhysicsConfig Cfg>
consteval size_t common_record_bytes() noexcept {
    return sizeof(PCoreB) + sizeof(PDynB<Cfg>) + sizeof(PAuxB<Cfg>);
}

template<PhysicsConfig Cfg>
consteval size_t gas_record_bytes() noexcept {
    if constexpr (HasHydro<Cfg>) return sizeof(GasAllB<Cfg>);
    return 0;
}

template<PhysicsConfig Cfg>
consteval size_t star_record_bytes() noexcept {
    if constexpr (HasStellarEvolution<Cfg>) return sizeof(StarAllB<Cfg>);
    return 0;
}

template<PhysicsConfig Cfg>
consteval size_t bh_record_bytes() noexcept {
    if constexpr (HasBlackHoles<Cfg>) return sizeof(BHAllB<Cfg>);
    return 0;
}

template<PhysicsConfig Cfg>
consteval size_t legacy_fixed_record_bytes() noexcept {
    size_t typed = gas_record_bytes<Cfg>();
    typed = typed > star_record_bytes<Cfg>() ? typed : star_record_bytes<Cfg>();
    typed = typed > bh_record_bytes<Cfg>() ? typed : bh_record_bytes<Cfg>();
    return common_record_bytes<Cfg>() + typed;
}

template<PhysicsConfig Cfg>
constexpr size_t exchange_bytes(const PackedExchangeCounts& counts) noexcept {
    return sizeof(PackedExchangeHeader)
         + static_cast<size_t>(counts.n_particles) * common_record_bytes<Cfg>()
         + static_cast<size_t>(counts.n_gas) * gas_record_bytes<Cfg>()
         + static_cast<size_t>(counts.n_star) * star_record_bytes<Cfg>()
         + static_cast<size_t>(counts.n_bh) * bh_record_bytes<Cfg>();
}

template<PhysicsConfig Cfg>
constexpr PackedExchangeHeader make_exchange_header(const PackedExchangeCounts& counts) noexcept {
    PackedExchangeHeader hdr{};
    hdr.n_particles = counts.n_particles;
    hdr.n_gas = counts.n_gas;
    hdr.n_star = counts.n_star;
    hdr.n_bh = counts.n_bh;
    hdr.common_bytes = static_cast<uint64_t>(counts.n_particles) * common_record_bytes<Cfg>();
    hdr.gas_bytes = static_cast<uint64_t>(counts.n_gas) * gas_record_bytes<Cfg>();
    hdr.star_bytes = static_cast<uint64_t>(counts.n_star) * star_record_bytes<Cfg>();
    hdr.bh_bytes = static_cast<uint64_t>(counts.n_bh) * bh_record_bytes<Cfg>();
    hdr.total_bytes = sizeof(PackedExchangeHeader)
                    + hdr.common_bytes + hdr.gas_bytes + hdr.star_bytes + hdr.bh_bytes;
    return hdr;
}

template<PhysicsConfig Cfg>
struct PackedExchangeView {
    std::byte* data = nullptr;
    size_t bytes = 0;
    constexpr explicit operator bool() const noexcept {
        return data != nullptr && bytes >= sizeof(PackedExchangeHeader);
    }
};

template<PhysicsConfig Cfg>
struct PackedConstExchangeView {
    const std::byte* data = nullptr;
    size_t bytes = 0;
    constexpr explicit operator bool() const noexcept {
        return data != nullptr && bytes >= sizeof(PackedExchangeHeader);
    }
};

template<PhysicsConfig Cfg>
inline PackedExchangeView<Cfg> make_exchange_view(void* data, size_t bytes) noexcept {
    return PackedExchangeView<Cfg>{static_cast<std::byte*>(data), bytes};
}

template<PhysicsConfig Cfg>
inline PackedConstExchangeView<Cfg> make_const_exchange_view(const void* data, size_t bytes) noexcept {
    return PackedConstExchangeView<Cfg>{static_cast<const std::byte*>(data), bytes};
}

template<PhysicsConfig Cfg>
inline PackedExchangeHeader load_exchange_header(const PackedConstExchangeView<Cfg>& view) noexcept {
    PackedExchangeHeader hdr{};
    if (!view || view.bytes < sizeof(hdr)) return hdr;
    std::memcpy(&hdr, view.data, sizeof(hdr));
    return hdr;
}

template<PhysicsConfig Cfg>
inline PackedExchangeHeader load_exchange_header(const PackedExchangeView<Cfg>& view) noexcept {
    return load_exchange_header<Cfg>(make_const_exchange_view<Cfg>(view.data, view.bytes));
}

template<PhysicsConfig Cfg>
inline void store_exchange_header(const PackedExchangeView<Cfg>& view,
                                  const PackedExchangeHeader& hdr) noexcept {
    assert(view);
    assert(view.bytes >= sizeof(hdr));
    std::memcpy(view.data, &hdr, sizeof(hdr));
}

template<PhysicsConfig Cfg>
inline bool validate_exchange_header(const PackedConstExchangeView<Cfg>& view,
                                     const PackedExchangeHeader& hdr) noexcept {
    if (!view) return false;
    if (hdr.magic != detail::packed_exchange_magic) return false;
    if (hdr.version != detail::packed_exchange_version) return false;
    if (hdr.common_bytes != static_cast<uint64_t>(hdr.n_particles) * common_record_bytes<Cfg>()) return false;
    if (hdr.gas_bytes != static_cast<uint64_t>(hdr.n_gas) * gas_record_bytes<Cfg>()) return false;
    if (hdr.star_bytes != static_cast<uint64_t>(hdr.n_star) * star_record_bytes<Cfg>()) return false;
    if (hdr.bh_bytes != static_cast<uint64_t>(hdr.n_bh) * bh_record_bytes<Cfg>()) return false;
    if (hdr.total_bytes != sizeof(PackedExchangeHeader)
                         + hdr.common_bytes + hdr.gas_bytes + hdr.star_bytes + hdr.bh_bytes) return false;
    if (view.bytes < hdr.total_bytes) return false;
    return true;
}

template<PhysicsConfig Cfg>
inline bool validate_exchange_header(const PackedExchangeView<Cfg>& view,
                                     const PackedExchangeHeader& hdr) noexcept {
    return validate_exchange_header<Cfg>(make_const_exchange_view<Cfg>(view.data, view.bytes), hdr);
}

template<PhysicsConfig Cfg>
inline const std::byte* common_block_ptr(const PackedConstExchangeView<Cfg>& view) noexcept {
    return view.data + sizeof(PackedExchangeHeader);
}

template<PhysicsConfig Cfg>
inline std::byte* common_block_ptr(const PackedExchangeView<Cfg>& view) noexcept {
    return view.data + sizeof(PackedExchangeHeader);
}

template<PhysicsConfig Cfg>
inline const std::byte* gas_block_ptr(const PackedConstExchangeView<Cfg>& view,
                                      const PackedExchangeHeader& hdr) noexcept {
    return common_block_ptr<Cfg>(view) + hdr.common_bytes;
}

template<PhysicsConfig Cfg>
inline std::byte* gas_block_ptr(const PackedExchangeView<Cfg>& view,
                                const PackedExchangeHeader& hdr) noexcept {
    return common_block_ptr<Cfg>(view) + hdr.common_bytes;
}

template<PhysicsConfig Cfg>
inline const std::byte* star_block_ptr(const PackedConstExchangeView<Cfg>& view,
                                       const PackedExchangeHeader& hdr) noexcept {
    return gas_block_ptr<Cfg>(view, hdr) + hdr.gas_bytes;
}

template<PhysicsConfig Cfg>
inline std::byte* star_block_ptr(const PackedExchangeView<Cfg>& view,
                                 const PackedExchangeHeader& hdr) noexcept {
    return gas_block_ptr<Cfg>(view, hdr) + hdr.gas_bytes;
}

template<PhysicsConfig Cfg>
inline const std::byte* bh_block_ptr(const PackedConstExchangeView<Cfg>& view,
                                     const PackedExchangeHeader& hdr) noexcept {
    return star_block_ptr<Cfg>(view, hdr) + hdr.star_bytes;
}

template<PhysicsConfig Cfg>
inline std::byte* bh_block_ptr(const PackedExchangeView<Cfg>& view,
                               const PackedExchangeHeader& hdr) noexcept {
    return star_block_ptr<Cfg>(view, hdr) + hdr.star_bytes;
}

template<PhysicsConfig Cfg>
PackedExchangeCounts count_range_exchange(const ParticleContainer<Cfg>& container,
                                          idx_t start,
                                          count_t count) noexcept {
    PackedExchangeCounts counts{};
    counts.n_particles = count;
    for (count_t i = 0; i < count; ++i) {
        const idx_t j = start + static_cast<idx_t>(i);
        switch (static_cast<ParticleType>(detail::bare_type(static_cast<uint8_t>(container.core()[j].type)))) {
            case ParticleType::Gas: ++counts.n_gas; break;
            case ParticleType::Star: ++counts.n_star; break;
            case ParticleType::BH: ++counts.n_bh; break;
            default: break;
        }
    }
    return counts;
}

template<PhysicsConfig Cfg>
PackedExchangeCounts count_scatter_exchange(const ParticleContainer<Cfg>& container,
                                            const idx_t* indices,
                                            count_t count) noexcept {
    PackedExchangeCounts counts{};
    counts.n_particles = count;
    for (count_t i = 0; i < count; ++i) {
        const idx_t j = indices[i];
        switch (static_cast<ParticleType>(detail::bare_type(static_cast<uint8_t>(container.core()[j].type)))) {
            case ParticleType::Gas: ++counts.n_gas; break;
            case ParticleType::Star: ++counts.n_star; break;
            case ParticleType::BH: ++counts.n_bh; break;
            default: break;
        }
    }
    return counts;
}

template<PhysicsConfig Cfg>
PackedExchangeView<Cfg> allocate_exchange_buffer(MemoryArena& arena,
                                                 const PackedExchangeCounts& counts,
                                                 const char* name = "B/tight_exchange") {
    const size_t bytes = exchange_bytes<Cfg>(counts);
    const ArenaBlock block = arena.allocate_temp_bytes(bytes, 64, name);
    return make_exchange_view<Cfg>(block.ptr, block.bytes);
}

template<PhysicsConfig Cfg>
bool pack_range_exchange(const ParticleContainer<Cfg>& container,
                         idx_t start,
                         count_t count,
                         const PackedExchangeView<Cfg>& view) noexcept {
    const PackedExchangeCounts counts = count_range_exchange(container, start, count);
    const PackedExchangeHeader hdr = make_exchange_header<Cfg>(counts);
    if (!view || view.bytes < hdr.total_bytes) return false;

    store_exchange_header<Cfg>(view, hdr);

    std::byte* common_out = common_block_ptr<Cfg>(view);
    std::byte* gas_out = gas_block_ptr<Cfg>(view, hdr);
    std::byte* star_out = star_block_ptr<Cfg>(view, hdr);
    std::byte* bh_out = bh_block_ptr<Cfg>(view, hdr);

    for (count_t i = 0; i < count; ++i) {
        const idx_t j = start + static_cast<idx_t>(i);
        std::memcpy(common_out, &container.core()[j], sizeof(PCoreB));
        common_out += sizeof(PCoreB);
        std::memcpy(common_out, &container.dyn()[j], sizeof(typename ParticleContainer<Cfg>::PDynT));
        common_out += sizeof(typename ParticleContainer<Cfg>::PDynT);
        std::memcpy(common_out, &container.aux()[j], sizeof(typename ParticleContainer<Cfg>::PAuxT));
        common_out += sizeof(typename ParticleContainer<Cfg>::PAuxT);

        const uint8_t bare = detail::bare_type(static_cast<uint8_t>(container.core()[j].type));
        const idx_t type_idx = static_cast<idx_t>(container.linkage()[j].type_idx);
        switch (static_cast<ParticleType>(bare)) {
            case ParticleType::Gas:
                if constexpr (HasHydro<Cfg>) {
                    std::memcpy(gas_out, &container.gas_all()[type_idx], sizeof(typename ParticleContainer<Cfg>::GasAllT));
                    gas_out += sizeof(typename ParticleContainer<Cfg>::GasAllT);
                } else {
                    return false;
                }
                break;
            case ParticleType::Star:
                if constexpr (HasStellarEvolution<Cfg>) {
                    std::memcpy(star_out, &container.star_all()[type_idx], sizeof(typename ParticleContainer<Cfg>::StarAllT));
                    star_out += sizeof(typename ParticleContainer<Cfg>::StarAllT);
                } else {
                    return false;
                }
                break;
            case ParticleType::BH:
                if constexpr (HasBlackHoles<Cfg>) {
                    std::memcpy(bh_out, &container.bh_all()[type_idx], sizeof(typename ParticleContainer<Cfg>::BHAllT));
                    bh_out += sizeof(typename ParticleContainer<Cfg>::BHAllT);
                } else {
                    return false;
                }
                break;
            default:
                break;
        }
    }

    return true;
}

template<PhysicsConfig Cfg>
bool pack_scatter_exchange(const ParticleContainer<Cfg>& container,
                           const idx_t* indices,
                           count_t count,
                           const PackedExchangeView<Cfg>& view) noexcept {
    const PackedExchangeCounts counts = count_scatter_exchange(container, indices, count);
    const PackedExchangeHeader hdr = make_exchange_header<Cfg>(counts);
    if (!view || view.bytes < hdr.total_bytes) return false;

    store_exchange_header<Cfg>(view, hdr);

    std::byte* common_out = common_block_ptr<Cfg>(view);
    std::byte* gas_out = gas_block_ptr<Cfg>(view, hdr);
    std::byte* star_out = star_block_ptr<Cfg>(view, hdr);
    std::byte* bh_out = bh_block_ptr<Cfg>(view, hdr);

    for (count_t i = 0; i < count; ++i) {
        const idx_t j = indices[i];
        std::memcpy(common_out, &container.core()[j], sizeof(PCoreB));
        common_out += sizeof(PCoreB);
        std::memcpy(common_out, &container.dyn()[j], sizeof(typename ParticleContainer<Cfg>::PDynT));
        common_out += sizeof(typename ParticleContainer<Cfg>::PDynT);
        std::memcpy(common_out, &container.aux()[j], sizeof(typename ParticleContainer<Cfg>::PAuxT));
        common_out += sizeof(typename ParticleContainer<Cfg>::PAuxT);

        const uint8_t bare = detail::bare_type(static_cast<uint8_t>(container.core()[j].type));
        const idx_t type_idx = static_cast<idx_t>(container.linkage()[j].type_idx);
        switch (static_cast<ParticleType>(bare)) {
            case ParticleType::Gas:
                if constexpr (HasHydro<Cfg>) {
                    std::memcpy(gas_out, &container.gas_all()[type_idx], sizeof(typename ParticleContainer<Cfg>::GasAllT));
                    gas_out += sizeof(typename ParticleContainer<Cfg>::GasAllT);
                } else {
                    return false;
                }
                break;
            case ParticleType::Star:
                if constexpr (HasStellarEvolution<Cfg>) {
                    std::memcpy(star_out, &container.star_all()[type_idx], sizeof(typename ParticleContainer<Cfg>::StarAllT));
                    star_out += sizeof(typename ParticleContainer<Cfg>::StarAllT);
                } else {
                    return false;
                }
                break;
            case ParticleType::BH:
                if constexpr (HasBlackHoles<Cfg>) {
                    std::memcpy(bh_out, &container.bh_all()[type_idx], sizeof(typename ParticleContainer<Cfg>::BHAllT));
                    bh_out += sizeof(typename ParticleContainer<Cfg>::BHAllT);
                } else {
                    return false;
                }
                break;
            default:
                break;
        }
    }

    return true;
}

template<PhysicsConfig Cfg>
PackedExchangeView<Cfg> pack_range_exchange(MemoryArena& arena,
                                            const ParticleContainer<Cfg>& container,
                                            idx_t start,
                                            count_t count,
                                            const char* name = "B/tight_exchange") {
    const PackedExchangeCounts counts = count_range_exchange(container, start, count);
    PackedExchangeView<Cfg> view = allocate_exchange_buffer<Cfg>(arena, counts, name);
    const bool ok = pack_range_exchange(container, start, count, view);
    assert(ok);
    (void)ok;
    return view;
}

template<PhysicsConfig Cfg>
PackedExchangeView<Cfg> pack_scatter_exchange(MemoryArena& arena,
                                              const ParticleContainer<Cfg>& container,
                                              const idx_t* indices,
                                              count_t count,
                                              const char* name = "B/tight_exchange") {
    const PackedExchangeCounts counts = count_scatter_exchange(container, indices, count);
    PackedExchangeView<Cfg> view = allocate_exchange_buffer<Cfg>(arena, counts, name);
    const bool ok = pack_scatter_exchange(container, indices, count, view);
    assert(ok);
    (void)ok;
    return view;
}

template<PhysicsConfig Cfg>
bool unpack_exchange(ParticleContainer<Cfg>& container,
                     idx_t j_dst,
                     const PackedConstExchangeView<Cfg>& view) noexcept {
    const PackedExchangeHeader hdr = load_exchange_header<Cfg>(view);
    if (!validate_exchange_header<Cfg>(view, hdr)) return false;
    if (j_dst < 0) return false;
    if (static_cast<count_t>(j_dst) + hdr.n_particles > container.capacity().n_max) return false;
    if (container.n_gas() + hdr.n_gas > container.capacity().n_max_gas) return false;
    if (container.n_star() + hdr.n_star > container.capacity().n_max_star) return false;
    if (container.n_bh() + hdr.n_bh > container.capacity().n_max_bh) return false;

    const std::byte* common_in = common_block_ptr<Cfg>(view);
    const std::byte* gas_in = gas_block_ptr<Cfg>(view, hdr);
    const std::byte* star_in = star_block_ptr<Cfg>(view, hdr);
    const std::byte* bh_in = bh_block_ptr<Cfg>(view, hdr);

    count_t next_gas = container.n_gas();
    count_t next_star = container.n_star();
    count_t next_bh = container.n_bh();

    for (count_t i = 0; i < hdr.n_particles; ++i) {
        const idx_t j = j_dst + static_cast<idx_t>(i);
        const PCoreB core = detail::load_trivial<PCoreB>(common_in);
        common_in += sizeof(PCoreB);
        const typename ParticleContainer<Cfg>::PDynT dyn =
            detail::load_trivial<typename ParticleContainer<Cfg>::PDynT>(common_in);
        common_in += sizeof(typename ParticleContainer<Cfg>::PDynT);
        const typename ParticleContainer<Cfg>::PAuxT aux =
            detail::load_trivial<typename ParticleContainer<Cfg>::PAuxT>(common_in);
        common_in += sizeof(typename ParticleContainer<Cfg>::PAuxT);

        std::memcpy(&container.core()[j], &core, sizeof(core));
        std::memcpy(&container.dyn()[j], &dyn, sizeof(dyn));
        std::memcpy(&container.aux()[j], &aux, sizeof(aux));

        switch (static_cast<ParticleType>(detail::bare_type(static_cast<uint8_t>(core.type)))) {
            case ParticleType::Gas:
                if constexpr (HasHydro<Cfg>) {
                    typename ParticleContainer<Cfg>::GasAllT gas{};
                    std::memcpy(&gas, gas_in, sizeof(gas));
                    gas_in += sizeof(gas);
                    std::memcpy(&container.gas_all()[next_gas], &gas, sizeof(gas));
                    container.linkage()[j].type_idx = static_cast<uint32_t>(next_gas);
                    ++next_gas;
                } else {
                    return false;
                }
                break;
            case ParticleType::Star:
                if constexpr (HasStellarEvolution<Cfg>) {
                    typename ParticleContainer<Cfg>::StarAllT star{};
                    std::memcpy(&star, star_in, sizeof(star));
                    star_in += sizeof(star);
                    std::memcpy(&container.star_all()[next_star], &star, sizeof(star));
                    container.linkage()[j].type_idx = static_cast<uint32_t>(next_star);
                    ++next_star;
                } else {
                    return false;
                }
                break;
            case ParticleType::BH:
                if constexpr (HasBlackHoles<Cfg>) {
                    typename ParticleContainer<Cfg>::BHAllT bh{};
                    std::memcpy(&bh, bh_in, sizeof(bh));
                    bh_in += sizeof(bh);
                    std::memcpy(&container.bh_all()[next_bh], &bh, sizeof(bh));
                    container.linkage()[j].type_idx = static_cast<uint32_t>(next_bh);
                    ++next_bh;
                } else {
                    return false;
                }
                break;
            default:
                container.linkage()[j].type_idx = 0;
                break;
        }
    }

    container.set_count(static_cast<count_t>(j_dst) + hdr.n_particles);
    container.set_count_gas(next_gas);
    container.set_count_star(next_star);
    container.set_count_bh(next_bh);
    return true;
}

template<PhysicsConfig Cfg>
bool unpack_exchange(ParticleContainer<Cfg>& container,
                     idx_t j_dst,
                     const void* data,
                     size_t bytes) noexcept {
    return unpack_exchange<Cfg>(container, j_dst, make_const_exchange_view<Cfg>(data, bytes));
}

} // namespace opg::layout_b_coarse

#endif // OPG_LAYOUT_B_COARSE_PACK_UNPACK_HPP
