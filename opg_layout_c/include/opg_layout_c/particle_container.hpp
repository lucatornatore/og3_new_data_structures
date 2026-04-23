/**
 * @file particle_container.hpp  (layout C)
 * @brief C particle container — Common + per-type arrays connected by
 *        implicit positional mapping. No cross-indices.
 *
 * The hard price of this elegance: the common, gas, star, and bh arrays
 * must be permuted atomically in lockstep during PH-key reshuffle, and
 * every mutation (gas→star, BH merge) must preserve the positional invariant.
 *
 * C maintains FOUR separate ArrayRegistries (one per array length), because
 * common, gas, star, and bh each have different lengths and different
 * permutations (all derived from the same sort, but applied to different
 * subsets). See reshuffle.hpp for the driver.
 */

#ifndef OPG_LAYOUT_C_PARTICLE_CONTAINER_HPP
#define OPG_LAYOUT_C_PARTICLE_CONTAINER_HPP

#include <opg_common/types/scalar_types.hpp>
#include <opg_common/types/physics_config.hpp>
#include <opg_common/memory/memory_arena.hpp>
#include "particles.hpp"
#include <opg_common/permutation/permutation.hpp>
#include "box_leaf.hpp"

namespace opg::layout_c {

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
// ParticleContainer<Cfg>  — C layout
// =============================================================================

template<PhysicsConfig Cfg>
class ParticleContainer {
public:
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
    // Note: NO PLinkage under C — the mapping is positional.
    PCore* core() noexcept { return core_; }
    PDyn*  dyn()  noexcept { return dyn_;  }
    PTime* time() noexcept { return time_; }
    PMeta* meta() noexcept { return meta_; }

    const PCore* core() const noexcept { return core_; }
    const PDyn*  dyn()  const noexcept { return dyn_;  }
    const PTime* time() const noexcept { return time_; }
    const PMeta* meta() const noexcept { return meta_; }

    // ========== optional common ==============================================
    PLeap*      leap()      noexcept requires HasLeapfrog<Cfg>        { return leap_storage_.p; }
    PPotential* potential() noexcept requires HasPotentialOutput<Cfg> { return pot_storage_.p;  }

    // ========== gas ==========================================================
    GasCore* gas_core() noexcept requires HasHydro<Cfg> { return gas_core_storage_.p; }
    GasGrad* gas_grad() noexcept requires HasSPH<Cfg>           { return gas_grad_storage_.p; }
    GasMag*  gas_mag()  noexcept requires HasMagnetic<Cfg>      { return gas_mag_storage_.p;  }
    GasMetalT* gas_metal() noexcept requires HasMetalSpecies<Cfg>     { return gas_metal_storage_.p; }
    GasSF*   gas_sf()   noexcept requires HasStarFormation<Cfg> { return gas_sf_storage_.p;  }
    GasMFMT* gas_mfm()  noexcept requires HasMFM<Cfg>           { return gas_mfm_storage_.p; }

    // ========== star =========================================================
    StarCoreT* star_core() noexcept requires HasStellarEvolution<Cfg> { return star_core_storage_.p; }
    StarMeta*  star_meta() noexcept requires HasStellarEvolution<Cfg> { return star_meta_storage_.p; }

    // ========== BH ===========================================================
    BHCore*  bh_core()  noexcept requires HasBlackHoles<Cfg>  { return bh_core_storage_.p;  }
    BHEnv*   bh_env()   noexcept requires HasBlackHoles<Cfg>  { return bh_env_storage_.p;   }
    BHRepos* bh_repos() noexcept requires HasBlackHoles<Cfg>  { return bh_repos_storage_.p; }
    BHSpin*  bh_spin()  noexcept requires HasBHSpin<Cfg>      { return bh_spin_storage_.p;  }
    BHKinFB* bh_kin_fb() noexcept requires HasBHKineticFB<Cfg>{ return bh_kin_fb_storage_.p; }

    // ========== counters =====================================================
    count_t count()      const noexcept { return n_part_; }
    count_t count_gas()  const noexcept { return n_gas_;  }
    count_t count_star() const noexcept { return n_star_; }
    count_t count_bh()   const noexcept { return n_bh_;   }

    void set_counts(count_t np, count_t ng, count_t ns, count_t nb) noexcept {
        n_part_ = np; n_gas_ = ng; n_star_ = ns; n_bh_ = nb;
    }
    void set_count_star(count_t ns) noexcept { n_star_ = ns; }
    void set_count_gas (count_t ng) noexcept { n_gas_  = ng; }

    const Capacity& capacity() const noexcept { return capacity_; }

    // Four independent registries — one per array-length class.
    ArrayRegistry<>& registry_common() noexcept { return registry_common_; }
    ArrayRegistry<>& registry_gas()    noexcept { return registry_gas_;    }
    ArrayRegistry<>& registry_star()   noexcept { return registry_star_;   }
    ArrayRegistry<>& registry_bh()     noexcept { return registry_bh_;     }

    size_t total_bytes_allocated() const noexcept { return bytes_total_; }

    static consteval const char* layout_name() noexcept {
        return "C (Common + per-type AoSS, implicit positional mapping)";
    }

private:
    void allocate_all() {
        auto& a = *arena_;
        const count_t N  = capacity_.n_max;
        const count_t Ng = capacity_.n_max_gas;
        const count_t Ns = capacity_.n_max_star;
        const count_t Nb = capacity_.n_max_bh;

        // Common (no PLinkage under C).
        core_ = a.allocate<PCore>(N, "PCore");
        dyn_  = a.allocate<PDyn> (N, "PDyn");
        time_ = a.allocate<PTime>(N, "PTime");
        meta_ = a.allocate<PMeta>(N, "PMeta");
        bytes_total_ += N * (sizeof(PCore) + sizeof(PDyn) + sizeof(PTime) + sizeof(PMeta));

        if constexpr (HasLeapfrog<Cfg>) {
            leap_storage_.p = a.allocate<PLeap>(N, "PLeap");
            bytes_total_ += N * sizeof(PLeap);
        }
        if constexpr (HasPotentialOutput<Cfg>) {
            pot_storage_.p = a.allocate<PPotential>(N, "PPotential");
            bytes_total_ += N * sizeof(PPotential);
        }

        // Gas arrays.
        if constexpr (HasHydro<Cfg>) {
            gas_core_storage_.p = a.allocate<GasCore>(Ng, "GasCore");
            bytes_total_ += Ng * sizeof(GasCore);
        }
        if constexpr (HasSPH<Cfg>)           { gas_grad_storage_.p  = a.allocate<GasGrad> (Ng, "GasGrad");  bytes_total_ += Ng * sizeof(GasGrad); }
        if constexpr (HasMagnetic<Cfg>)      { gas_mag_storage_.p   = a.allocate<GasMag>  (Ng, "GasMag");   bytes_total_ += Ng * sizeof(GasMag);  }
        if constexpr (HasMetalSpecies<Cfg>)        { gas_metal_storage_.p = a.allocate<GasMetalT>(Ng, "GasMetal"); bytes_total_ += Ng * sizeof(GasMetalT); }
        if constexpr (HasStarFormation<Cfg>) { gas_sf_storage_.p    = a.allocate<GasSF>   (Ng, "GasSF");    bytes_total_ += Ng * sizeof(GasSF);   }
        if constexpr (HasMFM<Cfg>)           { gas_mfm_storage_.p   = a.allocate<GasMFMT> (Ng, "GasMFM");   bytes_total_ += Ng * sizeof(GasMFMT); }

        // Star arrays.
        if constexpr (HasStellarEvolution<Cfg>) {
            star_core_storage_.p = a.allocate<StarCoreT>(Ns, "StarCore");
            star_meta_storage_.p = a.allocate<StarMeta>(Ns, "StarMeta");
            bytes_total_ += Ns * (sizeof(StarCoreT) + sizeof(StarMeta));
        }

        // BH arrays.
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
     * Register every array in its length-matched registry.
     *
     * This is the SINGLE PLACE where the reshuffle contract is declared.
     * Every component array that participates in the positional invariant
     * MUST be registered here. Missing one is the canonical N2 bug — but it
     * is now localised to this one function, not scattered across every
     * reshuffle call site.
     */
    void register_all_arrays() noexcept {
        const count_t N  = capacity_.n_max;
        const count_t Ng = capacity_.n_max_gas;
        const count_t Ns = capacity_.n_max_star;
        const count_t Nb = capacity_.n_max_bh;

        // Common.
        registry_common_.register_array(core_, N, "PCore");
        registry_common_.register_array(dyn_,  N, "PDyn");
        registry_common_.register_array(time_, N, "PTime");
        registry_common_.register_array(meta_, N, "PMeta");
        if constexpr (HasLeapfrog<Cfg>)        registry_common_.register_array(leap_storage_.p, N, "PLeap");
        if constexpr (HasPotentialOutput<Cfg>) registry_common_.register_array(pot_storage_.p,  N, "PPotential");

        // Gas.
        if constexpr (HasHydro<Cfg>) registry_gas_.register_array(gas_core_storage_.p, Ng, "GasCore");
        if constexpr (HasSPH<Cfg>)           registry_gas_.register_array(gas_grad_storage_.p,  Ng, "GasGrad");
        if constexpr (HasMagnetic<Cfg>)      registry_gas_.register_array(gas_mag_storage_.p,   Ng, "GasMag");
        if constexpr (HasMetalSpecies<Cfg>)        registry_gas_.register_array(gas_metal_storage_.p, Ng, "GasMetal");
        if constexpr (HasStarFormation<Cfg>) registry_gas_.register_array(gas_sf_storage_.p,    Ng, "GasSF");
        if constexpr (HasMFM<Cfg>)           registry_gas_.register_array(gas_mfm_storage_.p,   Ng, "GasMFM");

        // Star.
        if constexpr (HasStellarEvolution<Cfg>) {
            registry_star_.register_array(star_core_storage_.p, Ns, "StarCore");
            registry_star_.register_array(star_meta_storage_.p, Ns, "StarMeta");
        }

        // BH.
        if constexpr (HasBlackHoles<Cfg>) {
            registry_bh_.register_array(bh_core_storage_.p,  Nb, "BHCore");
            registry_bh_.register_array(bh_env_storage_.p,   Nb, "BHEnv");
            registry_bh_.register_array(bh_repos_storage_.p, Nb, "BHRepos");
        }
        if constexpr (HasBHSpin<Cfg>)       registry_bh_.register_array(bh_spin_storage_.p,   Nb, "BHSpin");
        if constexpr (HasBHKineticFB<Cfg>)  registry_bh_.register_array(bh_kin_fb_storage_.p, Nb, "BHKinFB");
    }

    // -------- state ---------------------------------------------------------
    MemoryArena* arena_ = nullptr;
    Capacity     capacity_{};

    ArrayRegistry<> registry_common_{};
    ArrayRegistry<> registry_gas_{};
    ArrayRegistry<> registry_star_{};
    ArrayRegistry<> registry_bh_{};

    count_t n_part_ = 0, n_gas_ = 0, n_star_ = 0, n_bh_ = 0;
    size_t  bytes_total_ = 0;

    PCore* core_ = nullptr;
    PDyn*  dyn_  = nullptr;
    PTime* time_ = nullptr;
    PMeta* meta_ = nullptr;

    [[no_unique_address]] optional_ptr<HasLeapfrog<Cfg>,        PLeap>      leap_storage_;
    [[no_unique_address]] optional_ptr<HasPotentialOutput<Cfg>, PPotential> pot_storage_;

    [[no_unique_address]] optional_ptr<HasHydro<Cfg>, GasCore> gas_core_storage_;
    [[no_unique_address]] optional_ptr<HasSPH<Cfg>,           GasGrad>   gas_grad_storage_;
    [[no_unique_address]] optional_ptr<HasMagnetic<Cfg>,      GasMag>    gas_mag_storage_;
    [[no_unique_address]] optional_ptr<HasMetalSpecies<Cfg>,        GasMetalT> gas_metal_storage_;
    [[no_unique_address]] optional_ptr<HasStarFormation<Cfg>, GasSF>     gas_sf_storage_;
    [[no_unique_address]] optional_ptr<HasMFM<Cfg>,           GasMFMT>   gas_mfm_storage_;

    [[no_unique_address]] optional_ptr<HasStellarEvolution<Cfg>, StarCoreT> star_core_storage_;
    [[no_unique_address]] optional_ptr<HasStellarEvolution<Cfg>, StarMeta>  star_meta_storage_;

    [[no_unique_address]] optional_ptr<HasBlackHoles<Cfg>, BHCore>  bh_core_storage_;
    [[no_unique_address]] optional_ptr<HasBlackHoles<Cfg>, BHEnv>   bh_env_storage_;
    [[no_unique_address]] optional_ptr<HasBlackHoles<Cfg>, BHRepos> bh_repos_storage_;
    [[no_unique_address]] optional_ptr<HasBHSpin<Cfg>,     BHSpin>  bh_spin_storage_;
    [[no_unique_address]] optional_ptr<HasBHKineticFB<Cfg>,BHKinFB> bh_kin_fb_storage_;
};

} // namespace opg::layout_c

#endif // OPG_LAYOUT_C_PARTICLE_CONTAINER_HPP
