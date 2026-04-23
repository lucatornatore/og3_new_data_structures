/**
 * @file particle_container.hpp  (layout B, coarse)
 * @brief B (coarse) particle container: 3-way common split + one type-specific
 *        struct per type, connected by explicit PLinkageB cross-indices.
 *
 * Under B, the type-specific arrays (GasAllB, StarAllB, BHAllB) do NOT
 * participate in the PH-key reshuffle — they are re-sorted independently
 * when desired (compact_type_arrays, not implemented yet). Therefore only
 * the four common arrays need a registry for reshuffle_all().
 *
 * Compared to v2 layout_b (B', fine):
 *   - 3 common structs instead of 5+ (Core, Dyn, Aux vs Core, Dyn, Leap,
 *     Time, Meta, Potential).
 *   - 1 type-specific struct per type instead of 5+ (GasAll vs GasCore,
 *     GasGrad, GasMag, GasMetal, GasSF, …).
 *   - Total ~7 arrays instead of ~17.
 */

#ifndef OPG_LAYOUT_B_COARSE_PARTICLE_CONTAINER_HPP
#define OPG_LAYOUT_B_COARSE_PARTICLE_CONTAINER_HPP

#include <cstddef>
#include <cstring>

#include "opg_common/types/scalar_types.hpp"
#include "opg_common/types/physics_config.hpp"
#include "opg_common/memory/memory_arena.hpp"
#include "opg_common/permutation/permutation.hpp"
#include "opg_common/tree/box_leaf_base.hpp"
#include "particles.hpp"

namespace opg::layout_b_coarse {

using namespace opg::common;

// =============================================================================
// Capacity — per-container runtime counts
// =============================================================================

struct Capacity {
    count_t n_max      = 0;   // maximum common-array capacity (all types)
    count_t n_max_gas  = 0;   // maximum GasAllB capacity
    count_t n_max_star = 0;   // maximum StarAllB capacity (0 if stellar evolution off)
    count_t n_max_bh   = 0;   // maximum BHAllB capacity (0 if BHs off)
};

// =============================================================================
// BoxLeaf — B's leaf descriptor
// =============================================================================
// Under B there is no per-type-array begin index inside the leaf. Type-specific
// iteration goes through `linkage[j].type_idx`. So BoxLeaf adds nothing over
// BoxLeafBase — a simple alias would work, but wrapping it preserves the
// layout-specific namespace so that a misuse against layout C is a type error.

struct BoxLeaf : public BoxLeafBase {
    // Intentionally empty: BoxLeafBase already carries common_begin / common_end
    // and the per-type offset/count arrays into the common range. That is all
    // B needs — type-specific arrays are not leaf-sorted under B.
};

static_assert(std::is_trivially_copyable_v<BoxLeaf>);

// =============================================================================
// ParticleContainer<Cfg>
// =============================================================================

template<PhysicsConfig Cfg>
class ParticleContainer {
public:
    // Compile-time type aliases for templated structs.
    using PDynT    = PDynB<Cfg>;
    using PAuxT    = PAuxB<Cfg>;
    using GasAllT  = GasAllB<Cfg>;
    using StarAllT = StarAllB<Cfg>;
    using BHAllT   = BHAllB<Cfg>;

    ParticleContainer() = default;

    ParticleContainer(MemoryArena& arena, const Capacity& cap)
        : arena_(&arena), capacity_(cap)
    {
        allocate_all();
        register_all_arrays();
    }

    ParticleContainer(const ParticleContainer&) = delete;
    ParticleContainer& operator=(const ParticleContainer&) = delete;

    // ========== common-array accessors (always present) ======================
    PCoreB*    core()    noexcept { return core_;    }
    PDynT*     dyn()     noexcept { return dyn_;     }
    PAuxT*     aux()     noexcept { return aux_;     }
    PLinkageB* linkage() noexcept { return linkage_; }

    const PCoreB*    core()    const noexcept { return core_;    }
    const PDynT*     dyn()     const noexcept { return dyn_;     }
    const PAuxT*     aux()     const noexcept { return aux_;     }
    const PLinkageB* linkage() const noexcept { return linkage_; }

    // ========== type-specific accessors =======================================
    // Gas is always allocated (possibly size 0 under GravityOnly). Star and BH
    // accessors are `requires`-constrained — absent from the interface when
    // their physics is disabled. A misuse is a compile-time "no member".
    GasAllT* gas_all() noexcept { return gas_all_; }
    const GasAllT* gas_all() const noexcept { return gas_all_; }

    StarAllT* star_all() noexcept requires HasStellarEvolution<Cfg>
        { return star_storage_.p; }
    const StarAllT* star_all() const noexcept requires HasStellarEvolution<Cfg>
        { return star_storage_.p; }

    BHAllT* bh_all() noexcept requires HasBlackHoles<Cfg>
        { return bh_storage_.p; }
    const BHAllT* bh_all() const noexcept requires HasBlackHoles<Cfg>
        { return bh_storage_.p; }

    // ========== counts =========================================================
    count_t n_part() const noexcept { return n_part_; }
    count_t n_gas()  const noexcept { return n_gas_;  }
    count_t n_star() const noexcept { return n_star_; }
    count_t n_bh()   const noexcept { return n_bh_;   }

    void set_count(count_t n)      noexcept { n_part_ = n; }
    void set_count_gas(count_t n)  noexcept { n_gas_  = n; }
    void set_count_star(count_t n) noexcept { n_star_ = n; }
    void set_count_bh(count_t n)   noexcept { n_bh_   = n; }

    const Capacity& capacity() const noexcept { return capacity_; }

    // ========== registry access ===============================================
    // Under B there is exactly ONE registry that participates in the PH-key
    // reshuffle: the common-array registry. Type-specific arrays are not
    // permuted by PH-reshuffle under B.
    ArrayRegistry<>&       common_registry()       noexcept { return reg_common_; }
    const ArrayRegistry<>& common_registry() const noexcept { return reg_common_; }

    // ========== layout name (diagnostics) =====================================
    consteval static const char* layout_name() noexcept { return "B-coarse"; }

private:
    void allocate_all() {
        const count_t N    = capacity_.n_max;
        const count_t Ng   = capacity_.n_max_gas;
        const count_t Nstr = capacity_.n_max_star;
        const count_t Nbh  = capacity_.n_max_bh;

        core_    = arena_->allocate<PCoreB>(N,    "B/core");
        dyn_     = arena_->allocate<PDynT>(N,     "B/dyn");
        aux_     = arena_->allocate<PAuxT>(N,     "B/aux");
        linkage_ = arena_->allocate<PLinkageB>(N, "B/linkage");

        gas_all_ = arena_->allocate<GasAllT>(Ng, "B/gas_all");

        if constexpr (HasStellarEvolution<Cfg>) {
            star_storage_.p = arena_->allocate<StarAllT>(Nstr, "B/star_all");
        }
        if constexpr (HasBlackHoles<Cfg>) {
            bh_storage_.p = arena_->allocate<BHAllT>(Nbh, "B/bh_all");
        }
    }

    void register_all_arrays() noexcept {
        // Only common arrays participate in PH-key reshuffle.
        const count_t N = capacity_.n_max;
        reg_common_.register_array(core_,    N, "B/core");
        reg_common_.register_array(dyn_,     N, "B/dyn");
        reg_common_.register_array(aux_,     N, "B/aux");
        reg_common_.register_array(linkage_, N, "B/linkage");

        // Type-specific arrays are deliberately NOT registered in reg_common_.
        // Under B they are reshuffled independently (by a separate
        // compact_type_arrays step, not implemented yet) or simply tolerated
        // with stale ordering between full rebuilds.
    }

    MemoryArena* arena_   = nullptr;
    Capacity     capacity_{};

    // Active particle counts (runtime state)
    count_t n_part_ = 0;
    count_t n_gas_  = 0;
    count_t n_star_ = 0;
    count_t n_bh_   = 0;

    // Always-allocated arrays
    PCoreB*    core_    = nullptr;
    PDynT*     dyn_     = nullptr;
    PAuxT*     aux_     = nullptr;
    PLinkageB* linkage_ = nullptr;
    GasAllT*   gas_all_ = nullptr;

    // Conditionally-allocated arrays — zero bytes when physics disabled.
    [[no_unique_address]]
    optional_ptr<HasStellarEvolution<Cfg>, StarAllT> star_storage_;
    [[no_unique_address]]
    optional_ptr<HasBlackHoles<Cfg>,       BHAllT>   bh_storage_;

    // The single reshuffle registry.
    ArrayRegistry<> reg_common_;
};

} // namespace opg::layout_b_coarse

#endif // OPG_LAYOUT_B_COARSE_PARTICLE_CONTAINER_HPP
