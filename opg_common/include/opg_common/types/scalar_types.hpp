#ifndef OPG_COMMON_TYPES_SCALAR_TYPES_HPP
#define OPG_COMMON_TYPES_SCALAR_TYPES_HPP

#include <cstdint>
#include <cstddef>
#include <type_traits>

namespace opg::common {

using int8   = std::int8_t;
using int16  = std::int16_t;
using int32  = std::int32_t;
using int64  = std::int64_t;
using uint8  = std::uint8_t;
using uint16 = std::uint16_t;
using uint32 = std::uint32_t;
using uint64 = std::uint64_t;

#ifndef OPG_PRECISION
#define OPG_PRECISION 1
#endif

#ifndef OPG_ID_PRECISION
#define OPG_ID_PRECISION 1
#endif

#if OPG_PRECISION == 0
using real_t     = float;
using real_acc_t = float;
using real_pos_t = float;
#elif OPG_PRECISION == 1
using real_t     = float;
using real_acc_t = double;
using real_pos_t = double;
#else
using real_t     = double;
using real_acc_t = double;
using real_pos_t = double;
#endif

using pos_t        = real_pos_t;
using mass_t       = real_t;
using vel_t        = real_t;
using acc_t        = real_acc_t;

using time_int_t   = uint32;
using time_float_t = real_t;

using pkey_t       = uint64;

#if OPG_ID_PRECISION == 0
using pid_t        = uint32;
#else
using pid_t        = uint64;
#endif

// FIX 1: Unified idx_t and count_t
#ifdef OPG_LARGE_DATASETS
    using count_t = std::uint64_t;
#else
    using count_t = std::uint32_t;
#endif

using local_count_t = uint32;
using idx_t         = count_t; 
using global_idx_t  = int64;

using dens_t       = real_t;
using pres_t       = real_t;
using enth_t       = real_t;
using temp_t       = real_t;
using hsml_t       = real_t;
using bfield_t     = real_t;

template<typename T, int N>
struct vec {
    T data[N];
    constexpr T&       operator[](int i)       noexcept { return data[i]; }
    constexpr const T& operator[](int i) const noexcept { return data[i]; }
};

using pos3_t    = vec<pos_t,    3>;
using vel3_t    = vec<vel_t,    3>;
using acc3_t    = vec<acc_t,    3>;
using bfield3_t = vec<bfield_t, 3>;

static_assert(std::is_trivially_copyable_v<pos3_t>);
static_assert(std::is_trivially_copyable_v<vel3_t>);
static_assert(std::is_trivially_copyable_v<acc3_t>);
static_assert(std::is_standard_layout_v<pos3_t>);

template<typename T>
concept ParticleComponent =
    std::is_trivially_copyable_v<T> &&
    std::is_standard_layout_v<T>;

template<typename T>
concept Scalar = std::is_arithmetic_v<T>;

#define OPG_ASSERT_PARTICLE_COMPONENT(...)                                     \
    static_assert(std::is_trivially_copyable_v<__VA_ARGS__>,                   \
                  #__VA_ARGS__ " must be trivially copyable (required for "    \
                  "memcpy and MPI_BYTE transfer in opg_common)");                \
    static_assert(std::is_standard_layout_v<__VA_ARGS__>,                      \
                  #__VA_ARGS__ " must be standard-layout (required for "       \
                  "portable field offsets across translation units)")

} // namespace opg::common

#endif // OPG_COMMON_TYPES_SCALAR_TYPES_HPP
