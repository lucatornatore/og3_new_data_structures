/**
 * @file particles.hpp  (layout C, fine)
 * @brief All particle structs for layout C — self-contained in this layout.
 *
 * Layout C (fine-grained, positional mapping): each access-pattern cluster
 * is a separate small struct, and each small struct lives in its own flat
 * array owned by the container. Common-side and type-specific-side arrays
 * are linked by POSITIONAL MAPPING through the box-leaf descriptor (no
 * stored cross-index): the j-th gas particle inside a leaf lives at
 * CommonCore[common_begin+j] and GasCore[gas_begin+j].
 *
 * Everything in this file lives in namespace opg::layout_c. There is NO
 * sharing of particle types with layout_b or layout_b_coarse: each layout
 * owns its own data model. The only shared dependencies are the generic
 * opg_common plumbing (PhysicsConfig, scalar types, optional_field, etc.).
 *
 * This file is intentionally a single-shot definition of every particle
 * struct the layout needs. If you are reading the layout in isolation,
 * everything particle-side lives here; the container, reshuffle, and
 * mutation headers just plug these types into arrays.
 */

#ifndef OPG_LAYOUT_C_PARTICLES_HPP
#define OPG_LAYOUT_C_PARTICLES_HPP

#include <opg_common/types/scalar_types.hpp>
#include <opg_common/types/physics_config.hpp>
#include <type_traits>

namespace opg::layout_c {

using namespace opg::common;

// =============================================================================
// Common-side fine-grained structs
// =============================================================================

struct alignas(64) PCore {
    pos3_t  pos;
    mass_t  mass;
    pkey_t  key;
    uint8   type;
    uint8   flags;
    int16   time_bin;
    idx_t   leaf_idx;
};
OPG_ASSERT_PARTICLE_COMPONENT(PCore);

namespace particle_flags {
    inline constexpr uint8 Active          = 1u << 0;
    inline constexpr uint8 Boundary        = 1u << 1;
    inline constexpr uint8 NeedsWakeup     = 1u << 2;
    inline constexpr uint8 ToBeDeleted     = 1u << 3;
    inline constexpr uint8 NewlyCreated    = 1u << 4;
    inline constexpr uint8 PendingMutation = 1u << 5;
}

struct alignas(32) PDyn {
    vel3_t  vel;
    acc3_t  grav_acc;
    real_t  old_acc;
};
OPG_ASSERT_PARTICLE_COMPONENT(PDyn);

struct alignas(32) PLeap {
    vel3_t  dp;
    vel3_t  vel_pred;
};
OPG_ASSERT_PARTICLE_COMPONENT(PLeap);

struct alignas(32) PTime {
    time_int_t  ti_begin;
    time_int_t  ti_current;
    int32       dt_step;
    real_t      grav_cost;
};
OPG_ASSERT_PARTICLE_COMPONENT(PTime);

struct alignas(32) PMeta {
    pid_t   id;
    real_t  num_ngb;
    idx_t   true_ngb;
};
OPG_ASSERT_PARTICLE_COMPONENT(PMeta);

struct PPotential {
    real_t  potential;
    real_t  pm_potential;
};
OPG_ASSERT_PARTICLE_COMPONENT(PPotential);

// =============================================================================
// Gas fine-grained structs
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
    real_t  metals[NMet > 0 ? NMet : 1];
    real_t  temperature;
    real_t  mass_res;
    real_t  egy_res;
    real_t  x_cold_cloud;
    real_t  egy_step;
    real_t  z_smooth;
    constexpr int num_species() const noexcept { return NMet; }
};
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
    real_t  cr_p_norm [NProton  > 0 ? NProton  : 1];
    real_t  cr_p_slope[NProton  > 0 ? NProton  : 1];
    real_t  cr_p_cut;
    real_t  cr_e_norm [NElectron > 0 ? NElectron : 1];
    real_t  cr_e_slope[NElectron > 0 ? NElectron : 1];
    real_t  cr_e_cut;
};
OPG_ASSERT_PARTICLE_COMPONENT(GasCR<1,1>);

// =============================================================================
// Star fine-grained structs
// =============================================================================

template<int NMet>
struct alignas(32) StarCore {
    real_t  stellar_age;
    real_t  last_chem_time;
    real_t  initial_mass;
    real_t  metals[NMet > 0 ? NMet : 1];
    real_t  weight;
    idx_t   pid;
    int32   chem_time_bin;
    real_t  mean_hsml;
    real_t  mean_rho;
    constexpr int num_species() const noexcept { return NMet; }
};
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

// =============================================================================
// BH fine-grained structs
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

} // namespace opg::layout_c

#endif // OPG_LAYOUT_C_PARTICLES_HPP
