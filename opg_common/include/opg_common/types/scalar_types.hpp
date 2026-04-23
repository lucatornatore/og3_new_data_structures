/**
 * @file scalar_types.hpp
 * @brief Internal scalar types for precision control.
 *
 * All field types throughout the codebase should use these aliases so that
 * single/double/mixed precision can be selected at build time via
 * OPG_PRECISION = 0 (single), 1 (mixed), 2 (double).
 *
 * Default is mixed: single-precision field storage, double-precision
 * accumulators and positions.
 */

#ifndef OPG_COMMON_TYPES_SCALAR_TYPES_HPP
#define OPG_COMMON_TYPES_SCALAR_TYPES_HPP

#include <cstdint>
#include <cstddef>
#include <type_traits>

namespace opg::common {

// =============================================================================
// Fixed-width integer types
// =============================================================================

using int8   = std::int8_t;
using int16  = std::int16_t;
using int32  = std::int32_t;
using int64  = std::int64_t;
using uint8  = std::uint8_t;
using uint16 = std::uint16_t;
using uint32 = std::uint32_t;
using uint64 = std::uint64_t;

// =============================================================================
// Precision selectors
// =============================================================================

#ifndef OPG_PRECISION
#define OPG_PRECISION 1  // 0=single, 1=mixed, 2=double
#endif

#ifndef OPG_ID_PRECISION
#define OPG_ID_PRECISION 1  // 0=32-bit IDs, 1=64-bit IDs
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

// =============================================================================
// Semantic aliases
// =============================================================================

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

using count_t       = uint64;
using local_count_t = uint32;
using idx_t         = int32;
using global_idx_t  = int64;

using dens_t       = real_t;
using pres_t       = real_t;
using enth_t       = real_t;
using temp_t       = real_t;
using hsml_t       = real_t;
using bfield_t     = real_t;

// =============================================================================
// POD vector (no constructors, keeps trivial-copyability)
// =============================================================================

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

// =============================================================================
// Concepts (C++20)
// =============================================================================

template<typename T>
concept ParticleComponent =
    std::is_trivially_copyable_v<T> &&
    std::is_standard_layout_v<T>;

template<typename T>
concept Scalar = std::is_arithmetic_v<T>;

// =============================================================================
// ASSERT_PARTICLE_COMPONENT — single tripwire covering both traits
// =============================================================================
//
// `static_assert(std::is_trivially_copyable_v<T>);`
// `static_assert(std::is_standard_layout_v<T>);`
// were spelled out separately after every particle struct definition. That
// was noisy, easy to forget one half of, and did in fact get half-forgotten
// in particle_gas.hpp and particle_star_bh.hpp. This macro asserts both
// traits in one line, and emits distinct diagnostics so a failure tells you
// which property broke.
//
// Why a macro and not a `static_assert(ParticleComponent<T>)`? Because a
// single combined concept assert would fail with "ParticleComponent not
// satisfied" without saying which half broke. The two sub-asserts below
// give a pointed diagnostic: "not trivially copyable" vs "not standard
// layout" — which is usually what you need to diagnose the fix quickly.
//
// Usage:
//
//   struct alignas(32) GasCore { ... };
//   OPG_ASSERT_PARTICLE_COMPONENT(GasCore);
//
//   // Templated types with multiple template args are OK — the macro is
//   // variadic so `GasCR<1, 1>` does not trip the preprocessor comma split:
//   OPG_ASSERT_PARTICLE_COMPONENT(GasCR<1, 1>);

#define OPG_ASSERT_PARTICLE_COMPONENT(...)                                     \
    static_assert(std::is_trivially_copyable_v<__VA_ARGS__>,                   \
                  #__VA_ARGS__ " must be trivially copyable (required for "    \
                  "memcpy and MPI_BYTE transfer in opg_common)");                \
    static_assert(std::is_standard_layout_v<__VA_ARGS__>,                      \
                  #__VA_ARGS__ " must be standard-layout (required for "       \
                  "portable field offsets across translation units)")

} // namespace opg::common

#endif // OPG_COMMON_TYPES_SCALAR_TYPES_HPP
