/* =========================================================================
 * fields.h
 *
 * "Atom" struct definitions shared by layouts B', B_coarse, and C. Layout
 * A does NOT use these atoms — it declares fat per-type structs that
 * inline equivalent fields.
 *
 * Physics configuration (OpenGadget3 flags):
 *   HYDRO=SPH, PMGRID, STAR_FORMATION (SFR), COOLING,
 *   STELLAR_EVOLUTION (LT_STELLAREVOLUTION), BLACK_HOLES, STELLARAGE,
 *   WINDS, CONDUCTION, LT_METAL_COOLING_WAL, LT_NMet=16.
 *
 * Alignment policy:
 *   Every atom uses  `alignas(OG3_ATOM_ALIGN)`  (set in bench_common.h,
 *   default 64). NO explicit padding arrays are used — the compiler
 *   rounds sizeof(T) up to a multiple of OG3_ATOM_ALIGN. The user can
 *   change the macro to trade memory for cache-line alignment.
 *
 * Raw (pre-alignment) byte counts are quoted in each atom's comment.
 * Actual sizeof(T) at OG3_ATOM_ALIGN = 64 is indicated too.
 *
 * Precision model (matches OG3 prototype):
 *   MyLongDouble, MyAtLeastDouble, MyDouble  -> double (8 B)
 *   MyFloat                                    -> float  (4 B)
 *   MyIDType                                   -> uint64_t
 *   integertime, peanokey                      -> 64-bit
 *   LT_NMetP                                   -> 16
 * ========================================================================= */

#ifndef FIELDS_H
#define FIELDS_H

#include "bench_common.h"

#include <stdalign.h>
#include <stdint.h>

/* OG3_ATOM is defined in bench_common.h and expands to alignas(OG3_ATOM_ALIGN). */

/* ----- fine common atoms ------------------------------------------------ */

/* PCore: hot fields read by almost every kernel.
 * raw 48 B, final sizeof depends on OG3_ATOM_ALIGN. */
typedef struct OG3_ATOM {
    pos_t      pos[3];       /* 24 B  Pos                              */
    real_t     mass;         /*  4 B  Mass                             */
    pkey_t     key;          /*  8 B  PH key                           */
    uint8_t    type;         /*  1 B  Type                             */
    uint8_t    flags;        /*  1 B                                   */
    int16_t    time_bin;     /*  2 B  TimeBin                          */
    int32_t    leaf_idx;     /*  4 B  box-leaf index                   */
} PCore;

/* PDyn: dynamical state (dp + GravAccel + GravPM + drift scalars).
 * raw 84 B, final sizeof depends on OG3_ATOM_ALIGN. */
typedef struct OG3_ATOM {
    real_t     vel[3];         /* 12 B  Vel                            */
    real_t     dp[3];          /* 12 B  dp                             */
    acc_t      grav_acc[3];    /* 24 B  GravAccel                      */
    acc_t      grav_pm[3];     /* 24 B  GravPM (PMGRID)                */
    real_t     old_acc;        /*  4 B  OldAcc                         */
    real_t     num_ngb;        /*  4 B  NumNgb                         */
    int32_t    true_ngb;       /*  4 B  TrueNGB                        */
} PDyn;

/* PTime: time integrators + gravity-cost vector.
 * raw 44 B, final sizeof depends on OG3_ATOM_ALIGN. */
typedef struct OG3_ATOM {
    time_int_t ti_begin;       /*  8 B  Ti_begstep                     */
    time_int_t ti_current;     /*  8 B  Ti_current                     */
    int32_t    dt_step;        /*  4 B  dt_step                        */
    real_t     grav_cost[6];   /* 24 B  GravCost[Gravcostlevels]       */
} PTime;

/* PMeta: identity + auxiliary flags.
 * raw 16 B, final sizeof depends on OG3_ATOM_ALIGN. */
typedef struct OG3_ATOM {
    ppid_t     id;             /*  8 B  ID                             */
    uint32_t   hsml_index;     /*  4 B                                 */
    uint32_t   flags_ext;      /*  4 B                                 */
} PMeta;

/* PLinkage: cross-index used by B' and B_coarse.
 * raw sizeof(count_t), final sizeof depends on OG3_ATOM_ALIGN. */
typedef struct OG3_ATOM {
    count_t    type_idx;       /* particle-array index                 */
} PLinkage;

/* ----- fine gas atoms (SphP under SPH + SFR + COOLING + LT_SE + WINDS + BH)
 * ------------------------------------------------------------------------ */

/* GasCore: primary hydro scalars.
 * raw 60 B, final sizeof depends on OG3_ATOM_ALIGN. */
typedef struct OG3_ATOM {
    real_t     density;            /*  4 B  Density                        */
    real_t     density_old;        /*  4 B  DensityOld                     */
    real_t     pressure;           /*  4 B  Pressure                       */
    real_t     entropy;            /*  4 B  Entropy                        */
    real_t     entropy_pred;       /*  4 B  EntropyPred                    */
    real_t     dt_entropy;         /*  4 B  DtEntropy                      */
    acc_t      dhsml_factor;       /*  8 B  DhsmlDensityFactor (HYDRO_SPH) */
    real_t     alpha;              /*  4 B  alpha                          */
    real_t     f_balsara;          /*  4 B  F_Balsara                      */
    real_t     temperature;        /*  4 B  Temperature (LT_SE)            */
    real_t     x_cold_cloud;       /*  4 B  XColdCloud (LT_SE)             */
    real_t     elec;               /*  4 B  elec (COOLING)                 */
    real_t     injected_bh_energy; /*  4 B  Injected_BH_Energy (BH)        */
    real_t     max_signal_vel;     /*  4 B  MaxSignalVel                   */
} GasCore;

/* GasGrad: vector quantities + integrators touched during force/drift.
 * raw 64 B */
typedef struct OG3_ATOM {
    real_t     vel_pred[3];        /* 12 B  VelPred                        */
    acc_t      hydro_acc[3];       /* 24 B  HydroAccel (MyAtLeastDouble)   */
    real_t     div_v;              /*  4 B  DivVel                         */
    real_t     curl_v;             /*  4 B  CurlVel (scalar form)          */
    acc_t      egy_step;           /*  8 B  EgyStep (MyAtLeastDouble)      */
    double     mstar;              /*  8 B  mstar                          */
} GasGrad;

/* GasMetal: metallicity vector (LT_NMetP=16).
 * raw 64 B */
typedef struct OG3_ATOM {
    real_t     metals[16];         /* 64 B  Metals[LT_NMetP] */
} GasMetal;

/* GasSF: SF / WINDS / BH-swallow fields.
 * raw 24 B, final sizeof depends on OG3_ATOM_ALIGN. */
typedef struct OG3_ATOM {
    ppid_t     swallow_id;         /*  8 B  SwallowID (BH)                 */
    real_t     sfr;                /*  4 B  Sfr                            */
    real_t     delay_time;         /*  4 B  DelayTime (WINDS)              */
    real_t     mass_res;           /*  4 B  MassRes (LT_SE)                */
    real_t     egy_res;            /*  4 B  EgyRes (LT_SE)                 */
} GasSF;

/* ----- fine star atoms (MetP under LT_SE + STELLARAGE) ------------------ */

/* StarCore: non-metal MetP scalars.
 * raw 28 B, final sizeof depends on OG3_ATOM_ALIGN. */
typedef struct OG3_ATOM {
    double     weight;             /*  8 B  weight                         */
    real_t     last_chem_time;     /*  4 B  LastChemTime                   */
    real_t     i_mass;             /*  4 B  iMass                          */
    uint32_t   pid;                /*  4 B  PID                            */
    real_t     stellar_age;        /*  4 B  StellarAge                     */
    int32_t    chem_time_bin;      /*  4 B  ChemTimeBin                    */
} StarCore;

/* StarMeta: Metals[16] for stars.
 * raw 64 B */
typedef struct OG3_ATOM {
    real_t     metals[16];         /* 64 B  Metals[LT_NMetP] */
} StarMeta;

/* ----- fine BH atoms (BHP under BLACK_HOLES) ---------------------------- */

/* BHCore: identity + mass + accretion scalars.
 * raw 48 B, final sizeof depends on OG3_ATOM_ALIGN. */
typedef struct OG3_ATOM {
    uint32_t   pid;                /*  4 B  PID                            */
    real_t     stellar_age;        /*  4 B  StellarAge                     */
    ppid_t     swallow_id;         /*  8 B  SwallowID                      */
    real_t     swallow_pot;        /*  4 B  SwallowPot                     */
    int32_t    count_progs;        /*  4 B  BH_CountProgs                  */
    double     bh_mass;            /*  8 B  BH_Mass (MyAtLeastDouble)      */
    real_t     bh_mdot;            /*  4 B  BH_Mdot                        */
    int32_t    time_bin_gas_ngb;   /*  4 B  BH_TimeBinGasNeighbor          */
    real_t     bh_density;         /*  4 B  BH_Density                     */
    real_t     bh_entropy;         /*  4 B  BH_Entropy                     */
} BHCore;

/* BHEnv: surrounding-gas environment (three phases).
 * raw 52 B, final sizeof depends on OG3_ATOM_ALIGN. */
typedef struct OG3_ATOM {
    real_t     surr_gas_vel[3];    /* 12 B  BH_SurroundingGasVel           */
    real_t     cold_density;       /*  4 B  BH_ColdDensity                 */
    real_t     cold_entropy;       /*  4 B  BH_ColdEntropy                 */
    real_t     surr_cold_vel[3];   /* 12 B  BH_SurroundingColdGasVel       */
    real_t     hot_density;        /*  4 B  BH_HotDensity                  */
    real_t     hot_entropy;        /*  4 B  BH_HotEntropy                  */
    real_t     surr_hot_vel[3];    /* 12 B  BH_SurroundingHotGasVel        */
} BHEnv;

/* BHRepos: accreted mass + repositioning + min-potential tracking.
 * raw 116 B, final sizeof depends on OG3_ATOM_ALIGN. */
typedef struct OG3_ATOM {
    real_t     accreted_mass;      /*  4 B  BH_accreted_Mass               */
    real_t     accreted_bh_mass;   /*  4 B  BH_accreted_BHMass             */
    real_t     accreted_mom[3];    /* 12 B  BH_accreted_momentum           */
    int32_t    flag_merged;        /*  4 B  BH_flag_merged                 */
    double     final_spin[3];      /* 24 B  BH_final_spin                  */
    real_t     surr_vel[3];        /* 12 B  BH_SurroundingVel              */
    real_t     surr_density;       /*  4 B  BH_SurroundingDensity          */
    real_t     sigma;              /*  4 B  BH_sigma                       */
    real_t     bmax;               /*  4 B  BH_bmax                        */
    real_t     tot_fb_eff;         /*  4 B  BH_TotalFeedbackEfficiency     */
    real_t     mean_hsml;          /*  4 B  mean_hsml                      */
    real_t     mean_rho;           /*  4 B  mean_rho                       */
    double     min_pot_pos[3];     /* 24 B  BH_MinPotPos (MyLongDouble)    */
    real_t     min_pot;            /*  4 B  BH_MinPot                      */
} BHRepos;

/* ------------------------------------------------------------------------ *
 * Size helpers.
 *
 * These query actual sizeof at runtime so callers report truthful bytes
 * regardless of OG3_ATOM_ALIGN.
 * ------------------------------------------------------------------------ */

static inline size_t gas_fine_bytes(void)
{
    return sizeof(GasCore) + sizeof(GasGrad) + sizeof(GasMetal) + sizeof(GasSF);
}
static inline size_t star_fine_bytes(void)
{
    return sizeof(StarCore) + sizeof(StarMeta);
}
static inline size_t bh_fine_bytes(void)
{
    return sizeof(BHCore) + sizeof(BHEnv) + sizeof(BHRepos);
}
static inline size_t common_fine_bytes(int with_linkage)
{
    return sizeof(PCore) + sizeof(PDyn) + sizeof(PTime) + sizeof(PMeta)
         + (with_linkage ? sizeof(PLinkage) : 0u);
}

#endif /* FIELDS_H */
