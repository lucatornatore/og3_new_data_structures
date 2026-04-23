/**
 * @file physics_config.hpp
 * @brief Particle-type enum, PhysicsConfig NTTP, physics-preset configs, concepts.
 *
 * `PhysicsConfig` is a structural type suitable for use as a C++20 non-type
 * template parameter of class type. All its predicates are `consteval`
 * (immediate functions): evaluating them on a non-constant-evaluated instance
 * is a compile error — the strongest enforcement of design note v4 § 7's
 * zero-overhead rule.
 *
 * `ArenaConfig` (elsewhere) stays `constexpr` because it legitimately
 * straddles compile and runtime.
 */

#ifndef OPG_COMMON_TYPES_PARTICLE_TYPES_HPP
#define OPG_COMMON_TYPES_PARTICLE_TYPES_HPP

#include "scalar_types.hpp"
#include <type_traits>
#include <utility>   // std::forward for opt_if

// =============================================================================
// Build-time configuration macros — single source of truth
// =============================================================================
//
// The C codebase uses compile flags like -DNMET=... to set the metal-species
// count once at build time. The C++ code mirrors that convention: the ONE
// place to set the metal-species count is the build flag. All PhysicsConfig
// presets then read from OPG_NMET_DEFAULT — they do not each carry their own
// independent 11.
//
// Override at build time:    g++ ... -DOPG_NMET=7 ...
// Default (matches C code):  11
//
// A preset that legitimately wants a different count (e.g. GravityOnly has
// no metals at all) explicitly sets .n_metal_species = 0; all other presets
// inherit from the macro below.

#ifndef OPG_NMET
#define OPG_NMET 11
#endif

#ifndef OPG_CR_PROTON_BINS
#define OPG_CR_PROTON_BINS 0
#endif

#ifndef OPG_CR_ELECTRON_BINS
#define OPG_CR_ELECTRON_BINS 0
#endif

namespace opg::common {

// The macro is surfaced as a namespaced constexpr so code that wants to read
// it can do so without consulting the preprocessor.
inline constexpr int OPG_NMET_DEFAULT          = OPG_NMET;
inline constexpr int OPG_CR_PROTON_BINS_DEFAULT   = OPG_CR_PROTON_BINS;
inline constexpr int OPG_CR_ELECTRON_BINS_DEFAULT = OPG_CR_ELECTRON_BINS;

// =============================================================================
// Particle type enumeration
// =============================================================================

enum class ParticleType : uint8 {
    Gas  = 0,
    DM1  = 1,
    DM2  = 2,
    DM3  = 3,
    Star = 4,
    BH   = 5,
    NumTypes = 6,
    Common   = 6,     // pseudo-type for "all types"
    Invalid  = 255
};

inline constexpr int NTYPES = static_cast<int>(ParticleType::NumTypes);

template<ParticleType PT>
inline constexpr bool is_gas_v = (PT == ParticleType::Gas);

template<ParticleType PT>
inline constexpr bool is_dm_v =
    (PT == ParticleType::DM1 || PT == ParticleType::DM2 || PT == ParticleType::DM3);

template<ParticleType PT>
inline constexpr bool is_star_v = (PT == ParticleType::Star);

template<ParticleType PT>
inline constexpr bool is_bh_v = (PT == ParticleType::BH);

template<ParticleType PT>
// Baryons: gas and stars. Black holes are NOT baryons — they are compact
// objects whose matter content (beyond the horizon) is not meaningfully
// described as a baryon count in this code's accounting. Use is_luminous_v
// or is_collisional_v if "things that aren't dark matter" is what you mean.
inline constexpr bool is_baryon_v =
    is_gas_v<PT> || is_star_v<PT>;

// =============================================================================
// PhysicsConfig — structural type for C++20 NTTP of class type
// =============================================================================
//
// Structural-type requirements (C++20 [temp.param]): all members public, no
// user-declared copy/move/dtor, every member's type itself structural.
// The `consteval` member functions do not disqualify structural status.

struct PhysicsConfig {
    // -------- hydro --------
    bool hydro_sph         = true;
    bool hydro_mfm         = false;

    // -------- gravity --------
    bool gravity_enabled   = true;
    bool gravity_pm        = true;
    bool gravity_adaptive  = false;

    // -------- sub-grid --------
    bool cooling           = true;
    bool star_formation    = true;
    bool stellar_winds     = false;

    // -------- MHD --------
    bool magnetic          = false;
    bool divb_cleaning     = false;

    // -------- chemistry / stellar evolution --------
    bool stellar_evolution = true;

    // Metal species count — single source of truth is the OPG_NMET build
    // flag (see top of file), surfaced here as OPG_NMET_DEFAULT. Presets
    // that legitimately differ (GravityOnly has none) override to 0.
    int  n_metal_species   = OPG_NMET_DEFAULT;

    // -------- black holes --------
    bool black_holes       = true;
    bool bh_spin           = false;
    bool bh_kinetic_fb     = false;

    // -------- cosmic rays --------
    bool cosmic_rays       = false;
    int  cr_proton_bins    = OPG_CR_PROTON_BINS_DEFAULT;
    int  cr_electron_bins  = OPG_CR_ELECTRON_BINS_DEFAULT;

    // -------- output / debug --------
    bool output_potential  = false;
    bool output_accel      = true;
    bool force_test        = false;

    // -------- derived queries: consteval — forbids runtime use --------

    consteval bool has_metals() const noexcept {
        return stellar_evolution || (star_formation && n_metal_species > 0);
    }

    consteval bool has_type_specific_data(ParticleType pt) const noexcept {
        switch (pt) {
            case ParticleType::Gas:  return true;
            case ParticleType::Star: return stellar_evolution;
            case ParticleType::BH:   return black_holes;
            default: return false;
        }
    }

    consteval bool has_leapfrog_predictor() const noexcept {
        return hydro_sph || hydro_mfm;
    }

    // NOTE: The hard-coded struct_size_estimate() helpers that lived here
    // have been removed. They were dead code (nothing called them) and,
    // as Luca rightly noted, would be unmaintainable if physics fields
    // grew: every new field would have required an ad-hoc number update
    // inside a comment-free consteval function.
    //
    // The right way to size an arena is from real types:
    //
    //    sizeof(GasCore) + (HasSPH<Cfg>      ? sizeof(GasGrad)              : 0)
    //                    + (HasMagnetic<Cfg> ? sizeof(GasMag)               : 0)
    //                    + (HasMetals<Cfg>   ? sizeof(GasMetal<Cfg.n_metal_species>) : 0)
    //                    + ...
    //
    // That helper belongs on the container (where the struct types are
    // already known) and should be added as part of the Phase-2 arena
    // pre-computation work (see README TODOs). It is deliberately NOT
    // placed on PhysicsConfig, because PhysicsConfig must not depend on
    // particle struct definitions (they include it).
};

// =============================================================================
// Preset configurations
// =============================================================================

namespace physics_configs {

inline constexpr PhysicsConfig GravityOnly{
    .hydro_sph = false, .hydro_mfm = false,
    .gravity_enabled = true, .gravity_pm = true, .gravity_adaptive = false,
    .cooling = false, .star_formation = false, .stellar_winds = false,
    .magnetic = false, .divb_cleaning = false,
    .stellar_evolution = false, .n_metal_species = 0,
    .black_holes = false, .bh_spin = false, .bh_kinetic_fb = false,
    .cosmic_rays = false, .cr_proton_bins = 0, .cr_electron_bins = 0,
    .output_potential = false, .output_accel = true, .force_test = false
};

inline constexpr PhysicsConfig StandardSPH{
    .hydro_sph = true, .hydro_mfm = false,
    .gravity_enabled = true, .gravity_pm = true, .gravity_adaptive = false,
    .cooling = true, .star_formation = true, .stellar_winds = true,
    .magnetic = false, .divb_cleaning = false,
    .stellar_evolution = true, .n_metal_species = OPG_NMET_DEFAULT,
    .black_holes = true, .bh_spin = false, .bh_kinetic_fb = false,
    .cosmic_rays = false, .cr_proton_bins = 0, .cr_electron_bins = 0,
    .output_potential = false, .output_accel = true, .force_test = false
};

inline constexpr PhysicsConfig FullMHD{
    .hydro_sph = true, .hydro_mfm = false,
    .gravity_enabled = true, .gravity_pm = true, .gravity_adaptive = false,
    .cooling = true, .star_formation = true, .stellar_winds = true,
    .magnetic = true, .divb_cleaning = true,
    .stellar_evolution = true, .n_metal_species = OPG_NMET_DEFAULT,
    .black_holes = true, .bh_spin = true, .bh_kinetic_fb = true,
    .cosmic_rays = false, .cr_proton_bins = 0, .cr_electron_bins = 0,
    .output_potential = true, .output_accel = true, .force_test = false
};

inline constexpr PhysicsConfig MFMHydro{
    .hydro_sph = false, .hydro_mfm = true,
    .gravity_enabled = true, .gravity_pm = true, .gravity_adaptive = false,
    .cooling = true, .star_formation = true, .stellar_winds = true,
    .magnetic = false, .divb_cleaning = false,
    .stellar_evolution = true, .n_metal_species = OPG_NMET_DEFAULT,
    .black_holes = true, .bh_spin = false, .bh_kinetic_fb = false,
    .cosmic_rays = false, .cr_proton_bins = 0, .cr_electron_bins = 0,
    .output_potential = false, .output_accel = true, .force_test = false
};

} // namespace physics_configs

// =============================================================================
// Physics concepts — the C++20-native "feature-is-enabled" checks
// =============================================================================

template<PhysicsConfig Cfg>
concept HasHydro = Cfg.hydro_sph || Cfg.hydro_mfm;

template<PhysicsConfig Cfg>
concept HasSPH = Cfg.hydro_sph;

template<PhysicsConfig Cfg>
concept HasMFM = Cfg.hydro_mfm;

template<PhysicsConfig Cfg>
concept HasMagnetic = Cfg.magnetic;

template<PhysicsConfig Cfg>
concept HasMagnetizedGas = Cfg.hydro_sph && Cfg.magnetic;

template<PhysicsConfig Cfg>
concept HasMetals = Cfg.has_metals();

template<PhysicsConfig Cfg>
concept HasStellarEvolution = Cfg.stellar_evolution;

template<PhysicsConfig Cfg>
concept HasStarFormation = Cfg.star_formation;

template<PhysicsConfig Cfg>
concept HasBlackHoles = Cfg.black_holes;

template<PhysicsConfig Cfg>
concept HasBHSpin = Cfg.black_holes && Cfg.bh_spin;

template<PhysicsConfig Cfg>
concept HasBHKineticFB = Cfg.black_holes && Cfg.bh_kinetic_fb;

template<PhysicsConfig Cfg>
concept HasCosmicRays = Cfg.cosmic_rays;

template<PhysicsConfig Cfg>
concept HasPotentialOutput = Cfg.output_potential;

template<PhysicsConfig Cfg>
concept HasLeapfrog = Cfg.has_leapfrog_predictor();

// =============================================================================
// optional_ptr — used with [[no_unique_address]] for zero-cost disabled fields
// =============================================================================

template<bool Enabled, typename T>
struct optional_ptr { T* p = nullptr; };

template<typename T>
struct optional_ptr<false, T> {};   // empty → 0 bytes with [[no_unique_address]]

// =============================================================================
// optional_field — value-holding counterpart for nested-struct optional members
// =============================================================================
//
// optional_ptr<Enabled,T>:   storage is a POINTER to T (used by container classes
//                            that own arrays of T). Disabled = empty = 0 bytes.
//
// optional_field<Enabled,T>: storage is a VALUE of T (used inside a particle
//                            struct to carry an optional sub-block of fields).
//                            Disabled = empty = 0 bytes.
//
// Access under the enabled specialisation is direct: `x.data.member`.
// User code must gate access with `if constexpr (EnabledConcept<Cfg>)` or
// with a `requires`-constrained accessor, exactly as for optional_ptr.
//
// Both specialisations are trivially-copyable and standard-layout iff T is,
// so types using optional_field remain valid ParticleComponents.

template<bool Enabled, typename T>
struct optional_field { T data; };

template<typename T>
struct optional_field<false, T> {};  // empty → 0 bytes with [[no_unique_address]]

// =============================================================================
// Helpers for accessing optional_field members ergonomically
// =============================================================================
//
// `optional_field<Enabled,T>` lets a struct carry a conditionally-present
// sub-block without paying storage when disabled. However, accessing it
// naïvely is clumsy at call sites:
//
//     if constexpr (HasMagnetic<Cfg>) {
//         record.mag.data.B = ...;                          // (a)
//     }
//
// Three layers of verbosity:
//   1. The `if constexpr (HasMagnetic<Cfg>)` gate is mandatory when the
//      feature may be disabled — without it, `record.mag.data.B` would fail
//      to compile because the Enabled=false specialisation has no `.data`.
//   2. The `.data.` hop is language plumbing, not intent.
//   3. The gate must cite the right concept name, which lives elsewhere.
//
// The helpers below address (2) and (3). (1) cannot go away without more
// invasive machinery; but when the mutation scope is large, the gate can be
// folded outside several accesses so it's not written per-line.
//
// ------ Read: value of an optional_field, with a fallback when disabled ---
//
//   opt_read(record.mag, bfield3_t{})  ==  record.mag.data when enabled,
//                                      ==  bfield3_t{}     when disabled.
//
// ------ Write / mutate: call a lambda only when enabled ----------------
//
//   opt_if(record.mag, [&](auto& m){ m.B = new_B; m.div_B = 0.0f; });
//
// The lambda is only instantiated (and only called) on the enabled branch.
// Under the disabled specialisation this expands to nothing. The style
// replaces the typical
//
//     if constexpr (HasMagnetic<Cfg>) {
//         record.mag.data.B     = new_B;
//         record.mag.data.div_B = 0.0f;
//     }
//
// with the visibly lighter
//
//     opt_if(record.mag, [&](auto& m){ m.B = new_B; m.div_B = 0.0f; });
//
// ------ On the surviving gate at the call site --------------------------
//
// When the feature is BOTH disabled AND the caller is in a context where
// the lambda body can compile even for the disabled type, `opt_if` alone
// suffices — the compiler instantiates the non-op version. When the lambda
// body refers to types that only exist when Enabled=true (e.g. a member of
// T only present under the enabled specialisation — but that is not what
// `optional_field<false,T>` does here, it is simply absent), then an outer
// `if constexpr` is still required as a type-existence gate. In the common
// case — reading/writing fields of T — `opt_if` alone is sufficient.

template<bool Enabled, typename T, typename U>
constexpr T opt_read(const optional_field<Enabled, T>& f, U&& fallback) noexcept {
    if constexpr (Enabled) return f.data;
    else                   return static_cast<T>(std::forward<U>(fallback));
}

// The noexcept-spec cannot refer to `f.data` directly because the disabled
// specialisation has no `data` member — the expression is evaluated before
// the `if constexpr` filters it. Forward-by-value of the lambda and leave
// noexcept unspecified; in practice the lambdas are trivial and the
// compiler infers it.

template<bool Enabled, typename T, typename F>
constexpr void opt_if(optional_field<Enabled, T>& f, F&& fn) {
    if constexpr (Enabled) fn(f.data);
    // else: no-op; the lambda is never instantiated for Enabled=false.
    (void)f;  // suppress -Wunused-parameter under the disabled spec
}

template<bool Enabled, typename T, typename F>
constexpr void opt_if(const optional_field<Enabled, T>& f, F&& fn) {
    if constexpr (Enabled) fn(f.data);
    (void)f;
}

} // namespace opg::common

#endif // OPG_COMMON_TYPES_PARTICLE_TYPES_HPP
