/**
 * @file particles.hpp  (layout B, coarse)
 * @brief All particle structs for layout B (coarse) — self-contained in this layout.
 *
 * Layout B (coarse): coarser-granularity counterpart of layout_b (B'). The
 * connection scheme is the same as B' (explicit cross-index via PLinkageB),
 * but the data is packed into fewer, wider structs per container:
 *
 *   Common side:
 *     PCoreB     — hot core (pos, mass, key, type, flags, time_bin, leaf_idx)
 *     PDynB<Cfg> — velocity + gravitational acceleration + old_acc
 *                  + optional leapfrog predictor
 *     PAuxB<Cfg> — ID + time-integration bookkeeping + neighbour counts
 *                  + optional potential output
 *     PLinkageB  — cross-index into the type-specific arrays
 *
 *   Type-specific side (one struct per type):
 *     GasAllB<Cfg>  — SPH hot loop + optional gradients, MHD, metals, SFR
 *     StarAllB<Cfg> — StarCore + StarMeta rolled into one
 *     BHAllB<Cfg>   — BHCore + BHEnv + BHRepos + optional spin, kinetic FB
 *
 * Total: 7 arrays for FullMHD vs ~17 under B'/C.
 *
 * Everything in this file lives in namespace opg::layout_b_coarse. There is
 * NO sharing of particle types with layout_b or layout_c: each layout owns
 * its own data model. The only shared dependencies are the generic
 * opg_common plumbing (PhysicsConfig, scalar types, optional_field, etc.).
 *
 * This file is intentionally a single-shot definition of every particle
 * struct the layout needs. If you are reading the layout in isolation,
 * everything particle-side lives here; the container, reshuffle, pack/unpack,
 * mutation, and compaction headers just plug these types into arrays.
 */

#ifndef OPG_LAYOUT_B_COARSE_PARTICLES_HPP
#define OPG_LAYOUT_B_COARSE_PARTICLES_HPP

#include <opg_common/types/scalar_types.hpp>
#include <opg_common/types/physics_config.hpp>
#include <type_traits>

namespace opg::layout_b_coarse {

using namespace opg::common;

// Zero-overhead storage block for compile-time-sized optional arrays.
template<int N>
struct FieldArray {
    real_t data[N > 0 ? N : 1];
    constexpr real_t& operator[](int i) noexcept { return data[i]; }
    constexpr const real_t& operator[](int i) const noexcept { return data[i]; }
    constexpr int size() const noexcept { return N; }
};

template<>
struct FieldArray<0> {
    constexpr int size() const noexcept { return 0; }
};

// =============================================================================
// PCoreB — hot core fields (always present, not templated)
// =============================================================================
// Sized to exactly one cache line under mixed precision: 24+4+8+1+1+2+4 = 44 B
// of data + 20 B alignas padding.
struct alignas(64) PCoreB {
    pos3_t  pos;         // 24 B  position (real_pos_t)
    mass_t  mass;        //  4 B  particle mass
    pkey_t  key;         //  8 B  Peano-Hilbert key
    uint8   type;        //  1 B  ParticleType; bit 7 = staged sentinel
    uint8   flags;       //  1 B  particle_flags bitfield
    int16   time_bin;    //  2 B  current time bin
    idx_t   leaf_idx;    //  4 B  containing box-leaf index
};
OPG_ASSERT_PARTICLE_COMPONENT(PCoreB);
static_assert(sizeof(PCoreB) == 64, "PCoreB must fit exactly one cache line");

namespace particle_flags {
    inline constexpr uint8 Active          = 1u << 0;
    inline constexpr uint8 Boundary        = 1u << 1;
    inline constexpr uint8 NeedsWakeup     = 1u << 2;
    inline constexpr uint8 ToBeDeleted     = 1u << 3;
    inline constexpr uint8 NewlyCreated    = 1u << 4;
    inline constexpr uint8 PendingMutation = 1u << 5;
}

// =============================================================================
// PDynB<Cfg> — kick/drift state + optional leapfrog predictor
// =============================================================================
template<PhysicsConfig Cfg>
struct alignas(64) PDynB {
    vel3_t  vel;           // 12 B  velocity
    acc3_t  grav_acc;      // 24 B  gravitational acceleration
    real_t  old_acc;       //  4 B  magnitude of previous acc; opening criterion
    real_t  _pad0;         //  4 B  explicit pad for 8-byte alignment of leap

    struct LeapState { vel3_t vel_pred; };
    [[no_unique_address]] optional_field<HasLeapfrog<Cfg>, LeapState> leap;
};

// =============================================================================
// PAuxB<Cfg> — bookkeeping: ID, time integration, neighbour counts, potential
// =============================================================================
template<PhysicsConfig Cfg>
struct alignas(64) PAuxB {
    pid_t       id;            // 8 B   unique particle ID
    time_int_t  ti_begin;      // 4 B   timestep start
    time_int_t  ti_current;    // 4 B   timestep current
    int32       dt_step;       // 4 B   timestep length
    real_t      grav_cost;     // 4 B   load-balancing cost metric
    real_t      num_ngb;       // 4 B   kernel estimate of neighbour count
    idx_t       true_ngb;      // 4 B   actual count from last loop

    [[no_unique_address]] optional_field<HasPotentialOutput<Cfg>, real_t> potential;
};

// =============================================================================
// PLinkageB — cross-index from common slot j to a slot in GasAllB/StarAllB/BHAllB
// =============================================================================
// 8 B per particle. NOT serialised during MPI exchange (per-task-local mapping;
// the receiver rebuilds linkage as it inserts particles into its own arrays).
struct alignas(8) PLinkageB {
    uint32_t type_idx;
    uint32_t _pad;
};
OPG_ASSERT_PARTICLE_COMPONENT(PLinkageB);
static_assert(sizeof(PLinkageB) == 8);

// =============================================================================
// Gas sub-structs (nested inside GasAllB)
// =============================================================================

struct alignas(32) GasCore {
    hsml_t  hsml;
    dens_t  density;
    pres_t  pressure;
    enth_t  entropy;
    enth_t  entropy_pred;
    acc3_t  hydro_acc;
    real_t  dt_entropy;
    real_t  max_signal_vel;
    real_t  sound_speed;
};
OPG_ASSERT_PARTICLE_COMPONENT(GasCore);

struct alignas(32) GasGrad {
    real_t  div_vel;
    real_t  curl_vel;
    real_t  dvel_matrix[3][3];
    real_t  dhsml_density_factor;
    real_t  alpha;
    real_t  f_balsara;
    real_t  dt_alpha;
};
OPG_ASSERT_PARTICLE_COMPONENT(GasGrad);

struct alignas(32) GasMag {
    bfield3_t B;
    bfield3_t B_pred;
    bfield3_t dt_B;
    real_t    div_B;
    bfield3_t mag_acc;
    bfield3_t mag_corr;
    real_t    phi;
    real_t    phi_pred;
    real_t    dt_phi;
    real_t    clean_vel;
};
OPG_ASSERT_PARTICLE_COMPONENT(GasMag);

template<int NMet>
struct alignas(32) GasMetal {
    static_assert(NMet > 0);
    real_t  metals[NMet];
    real_t  temperature;
    real_t  mass_res;
    real_t  egy_res;
    real_t  x_cold_cloud;
    real_t  egy_step;
    real_t  z_smooth;
    constexpr int num_species() const noexcept { return NMet; }
};

template<>
struct alignas(32) GasMetal<0> {
    real_t  temperature;
    real_t  mass_res;
    real_t  egy_res;
    real_t  x_cold_cloud;
    real_t  egy_step;
    real_t  z_smooth;
    constexpr int num_species() const noexcept { return 0; }
};
OPG_ASSERT_PARTICLE_COMPONENT(GasMetal<0>);
OPG_ASSERT_PARTICLE_COMPONENT(GasMetal<11>);

struct alignas(32) GasSF {
    real_t  sfr;
    real_t  delay_time;
    real_t  egy_sn;
    int32   multiphase;
    int32   out_part;
};
OPG_ASSERT_PARTICLE_COMPONENT(GasSF);

template<int D>
struct alignas(32) GasMFM {
    real_t  num_dens;
    real_t  dhsml_numdens_factor;
    real_t  dist_ngb_sqd_max;
    real_t  condition_number;
    static constexpr int NumVar = D + 2;
    real_t  alpha_slope[NumVar];
    real_t  internal_energy;
    real_t  internal_energy_pred;
    real_t  dt_internal_energy;
    real_t  max_vel_sqr_ngb;
    acc3_t  hydro_acc_extra;
    real_t  dt_entropy_extra;
};
OPG_ASSERT_PARTICLE_COMPONENT(GasMFM<3>);

template<int NProton, int NElectron>
struct alignas(32) GasCR {
    real_t  cr_p_pressure;
    real_t  cr_e_pressure;
    real_t  density_old;
    [[no_unique_address]] FieldArray<NProton>  cr_p_norm;
    [[no_unique_address]] FieldArray<NProton>  cr_p_slope;
    real_t  cr_p_cut;
    [[no_unique_address]] FieldArray<NElectron> cr_e_norm;
    [[no_unique_address]] FieldArray<NElectron> cr_e_slope;
    real_t  cr_e_cut;
};
OPG_ASSERT_PARTICLE_COMPONENT(GasCR<1,1>);

// =============================================================================
// GasAllB<Cfg> — all gas fields in one struct
// =============================================================================
template<PhysicsConfig Cfg>
struct alignas(64) GasAllB {
    GasCore core;
    [[no_unique_address]] optional_field<HasSPH<Cfg>,             GasGrad>                             grad;
    [[no_unique_address]] optional_field<HasMagnetic<Cfg>,        GasMag>                              mag;
    [[no_unique_address]] optional_field<HasThermoChemistry<Cfg>, GasMetal<Cfg.n_metal_species>> metals;
    [[no_unique_address]] optional_field<HasStarFormation<Cfg>,   GasSF>                               sfr;
    [[no_unique_address]] optional_field<HasMFM<Cfg>,             GasMFM<3>>                           mfm;
    [[no_unique_address]] optional_field<HasCosmicRays<Cfg>,      GasCR<Cfg.cr_proton_bins, Cfg.cr_electron_bins>> cr;
};

// =============================================================================
// Star sub-structs (nested inside StarAllB)
// =============================================================================

template<int NMet>
struct alignas(32) StarCore {
    static_assert(NMet > 0);
    real_t  stellar_age;
    real_t  last_chem_time;
    real_t  initial_mass;
    real_t  metals[NMet];
    real_t  weight;
    idx_t   pid;
    int32   chem_time_bin;
    real_t  mean_hsml;
    real_t  mean_rho;
    constexpr int num_species() const noexcept { return NMet; }
};

template<>
struct alignas(32) StarCore<0> {
    real_t  stellar_age;
    real_t  last_chem_time;
    real_t  initial_mass;
    real_t  weight;
    idx_t   pid;
    int32   chem_time_bin;
    real_t  mean_hsml;
    real_t  mean_rho;
    constexpr int num_species() const noexcept { return 0; }
};
OPG_ASSERT_PARTICLE_COMPONENT(StarCore<0>);
OPG_ASSERT_PARTICLE_COMPONENT(StarCore<11>);

struct alignas(32) StarMeta {
    real_t  z_age;
    real_t  z_age_llv;
    real_t  z_smooth;
    real_t  contrib_snii;
    real_t  contrib_snia;
    real_t  contrib_agb;
    int32   stellar_phase;
    real_t  init_z;
};
OPG_ASSERT_PARTICLE_COMPONENT(StarMeta);

template<PhysicsConfig Cfg>
struct alignas(64) StarAllB {
    StarCore<Cfg.n_metal_species> core;
    StarMeta                      meta;
};

// =============================================================================
// BH sub-structs (nested inside BHAllB)
// =============================================================================

struct alignas(32) BHCore {
    real_acc_t  bh_mass;
    real_t      bh_mdot;
    real_t      stellar_age;
    real_t      bh_density;
    real_t      bh_entropy;
    vel3_t      surrounding_vel;
    real_t      total_fb_eff;
    real_t      accreted_mass;
    real_acc_t  accreted_bh_mass;
    vel3_t      accreted_mom;
    pid_t       swallow_id;
    real_t      swallow_pot;
    int32       count_progs;
    int32       flag_merged;
    idx_t       pid;
    int32       time_bin_gas_ngb;
};
OPG_ASSERT_PARTICLE_COMPONENT(BHCore);

struct alignas(32) BHEnv {
    real_t  cold_density;
    real_t  cold_entropy;
    vel3_t  cold_gas_vel;
    real_t  hot_density;
    real_t  hot_entropy;
    vel3_t  hot_gas_vel;
    vel3_t  df_force;
    int32   df_ngb;
    real_t  sigma;
    real_t  bmax;
    real_t  mean_hsml;
    real_t  mean_rho;
};
OPG_ASSERT_PARTICLE_COMPONENT(BHEnv);

struct alignas(32) BHRepos {
    pos3_t  min_pot_pos;
    real_t  min_pot;
    pos3_t  star_dens_pos;
    real_t  star_dens;
    vel3_t  star_dens_vel;
};
OPG_ASSERT_PARTICLE_COMPONENT(BHRepos);

struct alignas(32) BHSpin {
    acc3_t      spin;
    real_t      eta_disk;
    real_t      alpha_disk;
    real_acc_t  mass_sum;
    real_acc_t  accr_mass_sum;
    real_acc_t  dt_sum;
    acc3_t      j_gas_sum;
    int32       accr_counter;
    int32       mode;
};
OPG_ASSERT_PARTICLE_COMPONENT(BHSpin);

struct alignas(32) BHKinFB {
    real_t  gas_mass_cone_up;
    real_t  gas_mass_cone_dn;
    int32   np_cone_up;
    int32   np_cone_dn;
    vel3_t  wind_direction;
    int16   flag_agn_on;
    int16   pad_;
    real_t  time_agn_last;
};
OPG_ASSERT_PARTICLE_COMPONENT(BHKinFB);

template<PhysicsConfig Cfg>
struct alignas(64) BHAllB {
    BHCore  core;
    BHEnv   env;
    BHRepos repos;
    [[no_unique_address]] optional_field<HasBHSpin<Cfg>,      BHSpin>  spin;
    [[no_unique_address]] optional_field<HasBHKineticFB<Cfg>, BHKinFB> kinfb;
};

// =============================================================================
// Tripwires: instantiate templated types against every preset
// =============================================================================
namespace _validate {
    using namespace opg::common::physics_configs;
    OPG_ASSERT_PARTICLE_COMPONENT(PDynB<GravityOnly>);
    OPG_ASSERT_PARTICLE_COMPONENT(PDynB<StandardSPH>);
    OPG_ASSERT_PARTICLE_COMPONENT(PDynB<FullMHD>);
    OPG_ASSERT_PARTICLE_COMPONENT(PAuxB<GravityOnly>);
    OPG_ASSERT_PARTICLE_COMPONENT(PAuxB<StandardSPH>);
    OPG_ASSERT_PARTICLE_COMPONENT(PAuxB<FullMHD>);
    OPG_ASSERT_PARTICLE_COMPONENT(GasAllB<StandardSPH>);
    OPG_ASSERT_PARTICLE_COMPONENT(GasAllB<FullMHD>);
    OPG_ASSERT_PARTICLE_COMPONENT(StarAllB<StandardSPH>);
    OPG_ASSERT_PARTICLE_COMPONENT(BHAllB<StandardSPH>);
    OPG_ASSERT_PARTICLE_COMPONENT(BHAllB<FullMHD>);
}

} // namespace opg::layout_b_coarse

#endif // OPG_LAYOUT_B_COARSE_PARTICLES_HPP
