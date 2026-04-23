/**
 * @file particle_container.hpp  (layout B')
 * @brief B' particle container: Common + clustered type-specific arrays.
 *
 * Connected by explicit PLinkage cross-indices. Under B', the type-specific
 * arrays do NOT participate in the PH-key reshuffle (they are re-sorted
 * independently when needed). Therefore the container carries ONE registry
 * (for common arrays only), whereas C carries four.
 *
 * Optional physics is carried via `[[no_unique_address]] optional_ptr<E,T>`,
 * zero bytes when disabled. Accessors are `requires`-constrained concepts,
 * so disabled features are NOT present as members — a misuse is "no member
 * named …" at compile time, not a runtime null deref.
 */

#ifndef OPG_LAYOUT_B_PARTICLE_CONTAINER_HPP
#define OPG_LAYOUT_B_PARTICLE_CONTAINER_HPP

#include <opg_common/types/scalar_types.hpp>
#include <opg_common/types/physics_config.hpp>
#include <opg_common/memory/memory_arena.hpp>
#include "particles.hpp"
#include <opg_common/permutation/permutation.hpp>
#include "box_leaf.hpp"

namespace opg::layout_b {

using namespace opg::common;

// =============================================================================
// Capacity — per-container runtime counts
// =============================================================================

struct Capacity {
    count_t n_max      = 0;
    count_t n_max_gas  = 0;
    count_t n_max_star = 0;
    count_t n_max_bh   = 0;
};

// =============================================================================
// ParticleContainer<Cfg> — Cfg is NTTP of class type
// =============================================================================

template<PhysicsConfig Cfg>
class ParticleContainer {
public:
    // Compile-time instantiated templated structs (A.7 fix).
    using GasMetalT = GasMetal<Cfg.n_metal_species>;
    using GasMFMT   = GasMFM<3>;
    using GasCRT    = GasCR<Cfg.cr_proton_bins, Cfg.cr_electron_bins>;
    using StarCoreT = StarCore<Cfg.n_metal_species>;

    ParticleContainer() = default;

    ParticleContainer(MemoryArena& arena, const Capacity& cap)
        : arena_(&arena), capacity_(cap)
    {
        allocate_all();
        register_all_arrays();
    }

    ParticleContainer(const ParticleContainer&) = delete;
    ParticleContainer& operator=(const ParticleContainer&) = delete;

    // ========== required common-array accessors ===============================
    PCore*    core()    noexcept { return core_;    }
    PDyn*     dyn()     noexcept { return dyn_;     }
    PTime*    time()    noexcept { return time_;    }
    PMeta*    meta()    noexcept { return meta_;    }
    PLinkage* linkage() noexcept { return linkage_; }

    const PCore*    core()    const noexcept { return core_;    }
    const PDyn*     dyn()     const noexcept { return dyn_;     }
    const PTime*    time()    const noexcept { return time_;    }
    const PMeta*    meta()    const noexcept { return meta_;    }
    const PLinkage* linkage() const noexcept { return linkage_; }

    // ========== optional common-array accessors ===============================
    PLeap*      leap()      noexcept requires HasLeapfrog<Cfg>        { return leap_storage_.p; }
    PPotential* potential() noexcept requires HasPotentialOutput<Cfg> { return pot_storage_.p;  }

    // ========== gas arrays ====================================================
    GasCore* gas_core() noexcept                              { return gas_core_; }
    GasGrad* gas_grad() noexcept requires HasSPH<Cfg>         { return gas_grad_storage_.p; }
    GasMag*  gas_mag()  noexcept requires HasMagnetic<Cfg>    { return gas_mag_storage_.p;  }
    GasMetalT* gas_metal() noexcept requires HasMetals<Cfg>   { return gas_metal_storage_.p; }
    GasSF*   gas_sf()   noexcept requires HasStarFormation<Cfg> { return gas_sf_storage_.p; }
    GasMFMT* gas_mfm()  noexcept requires HasMFM<Cfg>         { return gas_mfm_storage_.p; }

    // ========== star arrays ===================================================
    StarCoreT* star_core() noexcept requires HasStellarEvolution<Cfg> { return star_core_storage_.p; }
    StarMeta*  star_meta() noexcept requires HasStellarEvolution<Cfg> { return star_meta_storage_.p; }

    // ========== BH arrays =====================================================
    BHCore*  bh_core()  noexcept requires HasBlackHoles<Cfg>  { return bh_core_storage_.p;  }
    BHEnv*   bh_env()   noexcept requires HasBlackHoles<Cfg>  { return bh_env_storage_.p;   }
    BHRepos* bh_repos() noexcept requires HasBlackHoles<Cfg>  { return bh_repos_storage_.p; }
    BHSpin*  bh_spin()  noexcept requires HasBHSpin<Cfg>      { return bh_spin_storage_.p;  }
    BHKinFB* bh_kin_fb() noexcept requires HasBHKineticFB<Cfg>{ return bh_kin_fb_storage_.p; }

    // ========== counters ======================================================
    count_t count()      const noexcept { return n_part_; }
    count_t count_gas()  const noexcept { return n_gas_;  }
    count_t count_star() const noexcept { return n_star_; }
    count_t count_bh()   const noexcept { return n_bh_;   }

    void set_counts(count_t np, count_t ng, count_t ns, count_t nb) noexcept {
        n_part_ = np; n_gas_ = ng; n_star_ = ns; n_bh_ = nb;
    }

    const Capacity& capacity() const noexcept { return capacity_; }

    ArrayRegistry<>& registry() noexcept { return registry_; }

    size_t total_bytes_allocated() const noexcept { return bytes_total_; }

    static consteval const char* layout_name() noexcept {
        return "B' (Common + clustered type-specific, explicit cross-index)";
    }

private:
    void allocate_all() {
        auto& a = *arena_;
        const count_t N  = capacity_.n_max;
        const count_t Ng = capacity_.n_max_gas;
        const count_t Ns = capacity_.n_max_star;
        const count_t Nb = capacity_.n_max_bh;

        // Common (always).
        core_    = a.allocate<PCore>   (N, "PCore");
        dyn_     = a.allocate<PDyn>    (N, "PDyn");
        time_    = a.allocate<PTime>   (N, "PTime");
        meta_    = a.allocate<PMeta>   (N, "PMeta");
        linkage_ = a.allocate<PLinkage>(N, "PLinkage");

        bytes_total_ += N * (sizeof(PCore) + sizeof(PDyn) + sizeof(PTime)
                           + sizeof(PMeta) + sizeof(PLinkage));

        // Optional common.
        if constexpr (HasLeapfrog<Cfg>) {
            leap_storage_.p = a.allocate<PLeap>(N, "PLeap");
            bytes_total_ += N * sizeof(PLeap);
        }
        if constexpr (HasPotentialOutput<Cfg>) {
            pot_storage_.p = a.allocate<PPotential>(N, "PPotential");
            bytes_total_ += N * sizeof(PPotential);
        }

        // Gas (core always present; assumes Ng > 0 when hydro is enabled).
        gas_core_ = a.allocate<GasCore>(Ng, "GasCore");
        bytes_total_ += Ng * sizeof(GasCore);

        if constexpr (HasSPH<Cfg>) {
            gas_grad_storage_.p = a.allocate<GasGrad>(Ng, "GasGrad");
            bytes_total_ += Ng * sizeof(GasGrad);
        }
        if constexpr (HasMagnetic<Cfg>) {
            gas_mag_storage_.p = a.allocate<GasMag>(Ng, "GasMag");
            bytes_total_ += Ng * sizeof(GasMag);
        }
        if constexpr (HasMetals<Cfg>) {
            gas_metal_storage_.p = a.allocate<GasMetalT>(Ng, "GasMetal");
            bytes_total_ += Ng * sizeof(GasMetalT);
        }
        if constexpr (HasStarFormation<Cfg>) {
            gas_sf_storage_.p = a.allocate<GasSF>(Ng, "GasSF");
            bytes_total_ += Ng * sizeof(GasSF);
        }
        if constexpr (HasMFM<Cfg>) {
            gas_mfm_storage_.p = a.allocate<GasMFMT>(Ng, "GasMFM");
            bytes_total_ += Ng * sizeof(GasMFMT);
        }

        // Star.
        if constexpr (HasStellarEvolution<Cfg>) {
            star_core_storage_.p = a.allocate<StarCoreT>(Ns, "StarCore");
            star_meta_storage_.p = a.allocate<StarMeta>(Ns, "StarMeta");
            bytes_total_ += Ns * (sizeof(StarCoreT) + sizeof(StarMeta));
        }

        // BH.
        if constexpr (HasBlackHoles<Cfg>) {
            bh_core_storage_.p  = a.allocate<BHCore >(Nb, "BHCore");
            bh_env_storage_.p   = a.allocate<BHEnv  >(Nb, "BHEnv");
            bh_repos_storage_.p = a.allocate<BHRepos>(Nb, "BHRepos");
            bytes_total_ += Nb * (sizeof(BHCore) + sizeof(BHEnv) + sizeof(BHRepos));
        }
        if constexpr (HasBHSpin<Cfg>) {
            bh_spin_storage_.p = a.allocate<BHSpin>(Nb, "BHSpin");
            bytes_total_ += Nb * sizeof(BHSpin);
        }
        if constexpr (HasBHKineticFB<Cfg>) {
            bh_kin_fb_storage_.p = a.allocate<BHKinFB>(Nb, "BHKinFB");
            bytes_total_ += Nb * sizeof(BHKinFB);
        }
    }

    /**
     * Register ONLY the common arrays in the B' registry.
     * Type-specific arrays do not move with the PH-key reshuffle under B'.
     */
    void register_all_arrays() noexcept {
        const count_t N = capacity_.n_max;
        registry_.register_array(core_,    N, "PCore");
        registry_.register_array(dyn_,     N, "PDyn");
        registry_.register_array(time_,    N, "PTime");
        registry_.register_array(meta_,    N, "PMeta");
        registry_.register_array(linkage_, N, "PLinkage");

        if constexpr (HasLeapfrog<Cfg>) {
            registry_.register_array(leap_storage_.p, N, "PLeap");
        }
        if constexpr (HasPotentialOutput<Cfg>) {
            registry_.register_array(pot_storage_.p, N, "PPotential");
        }
    }

    MemoryArena* arena_ = nullptr;
    Capacity     capacity_{};
    ArrayRegistry<> registry_{};

    count_t n_part_ = 0, n_gas_ = 0, n_star_ = 0, n_bh_ = 0;
    size_t  bytes_total_ = 0;

    // Common pointers.
    PCore*    core_    = nullptr;
    PDyn*     dyn_     = nullptr;
    PTime*    time_    = nullptr;
    PMeta*    meta_    = nullptr;
    PLinkage* linkage_ = nullptr;

    // Optional common.
    [[no_unique_address]] optional_ptr<HasLeapfrog<Cfg>,        PLeap>      leap_storage_;
    [[no_unique_address]] optional_ptr<HasPotentialOutput<Cfg>, PPotential> pot_storage_;

    // Gas core is mandatory; optional gas parts are empty when disabled.
    GasCore* gas_core_ = nullptr;
    [[no_unique_address]] optional_ptr<HasSPH<Cfg>,           GasGrad>   gas_grad_storage_;
    [[no_unique_address]] optional_ptr<HasMagnetic<Cfg>,      GasMag>    gas_mag_storage_;
    [[no_unique_address]] optional_ptr<HasMetals<Cfg>,        GasMetalT> gas_metal_storage_;
    [[no_unique_address]] optional_ptr<HasStarFormation<Cfg>, GasSF>     gas_sf_storage_;
    [[no_unique_address]] optional_ptr<HasMFM<Cfg>,           GasMFMT>   gas_mfm_storage_;

    // Star.
    [[no_unique_address]] optional_ptr<HasStellarEvolution<Cfg>, StarCoreT> star_core_storage_;
    [[no_unique_address]] optional_ptr<HasStellarEvolution<Cfg>, StarMeta>  star_meta_storage_;

    // BH.
    [[no_unique_address]] optional_ptr<HasBlackHoles<Cfg>, BHCore>  bh_core_storage_;
    [[no_unique_address]] optional_ptr<HasBlackHoles<Cfg>, BHEnv>   bh_env_storage_;
    [[no_unique_address]] optional_ptr<HasBlackHoles<Cfg>, BHRepos> bh_repos_storage_;
    [[no_unique_address]] optional_ptr<HasBHSpin<Cfg>,     BHSpin>  bh_spin_storage_;
    [[no_unique_address]] optional_ptr<HasBHKineticFB<Cfg>,BHKinFB> bh_kin_fb_storage_;
};

} // namespace opg::layout_b

#endif // OPG_LAYOUT_B_PARTICLE_CONTAINER_HPP
