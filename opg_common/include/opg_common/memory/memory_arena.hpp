/**
 * @file memory_arena.hpp
 * @brief Single-allocation arena for particle data.
 *
 * Design:
 *   - One posix_memalign at MPI-task startup.
 *   - No fragmentation: all component arrays are carved from one buffer.
 *   - No std::vector: avoids allocator complexity and hidden reallocs.
 *   - Concept-constrained allocate<T>: only ParticleComponent types allowed.
 *   - No unconditional zero-fill. Caller chooses ArenaFillPolicy.
 *
 * ArenaConfig::n_max() and n_average() are `constexpr` (legitimately
 * straddle compile and runtime): an ArenaConfig may be a compile-time
 * constant in a test or built at runtime from MPI_Comm_size + params.
 */

#ifndef OPG_COMMON_MEMORY_ARENA_HPP
#define OPG_COMMON_MEMORY_ARENA_HPP

#include "../types/scalar_types.hpp"
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace opg::common {

// =============================================================================
// ArenaConfig — constexpr (straddles compile and runtime)
// =============================================================================

struct ArenaConfig {
    count_t total_particles   = 0;
    int     n_tasks           = 1;
    double  part_alloc_factor = 1.2;
    int     numa_node         = -1;
    bool    use_huge_pages    = false;
    bool    use_mpi_shared    = false;
    size_t  alignment         = 64;

    constexpr count_t n_average() const noexcept {
        return n_tasks > 0
            ? total_particles / static_cast<count_t>(n_tasks)
            : count_t{0};
    }

    constexpr count_t n_max() const noexcept {
        return static_cast<count_t>(n_average() * part_alloc_factor);
    }
};

// =============================================================================
// Fill-policy for the arena buffer at allocation time
// =============================================================================

enum class ArenaFillPolicy {
    None,     // no fill — fastest, caller must initialise
    Zero,     // zero-fill (legacy behaviour)
    Poison    // fill with a debug pattern
};

// =============================================================================
// MemoryArena
// =============================================================================

class MemoryArena {
public:
    MemoryArena() = default;

    MemoryArena(size_t total_bytes,
                const ArenaConfig& config,
                ArenaFillPolicy fill_policy = ArenaFillPolicy::None)
        : config_(config), fill_policy_(fill_policy)
    {
        allocate_arena(total_bytes, config);
    }

    ~MemoryArena() { deallocate_arena(); }

    MemoryArena(const MemoryArena&) = delete;
    MemoryArena& operator=(const MemoryArena&) = delete;

    MemoryArena(MemoryArena&& other) noexcept { move_from(std::move(other)); }
    MemoryArena& operator=(MemoryArena&& other) noexcept {
        if (this != &other) {
            deallocate_arena();
            move_from(std::move(other));
        }
        return *this;
    }

    /**
     * Allocate n contiguous T objects from the arena.
     * Only ParticleComponent (trivially-copyable + standard-layout) allowed.
     */
    template<ParticleComponent T>
    T* allocate(count_t n, const char* /*name*/ = "") {
        const size_t bytes = n * sizeof(T);
        const size_t align = alignof(T) > config_.alignment
                             ? alignof(T) : config_.alignment;

        // Align current_offset_ up to `align`.
        size_t misaligned = current_offset_ % align;
        if (misaligned) current_offset_ += (align - misaligned);

        if (current_offset_ + bytes > total_size_) {
            throw std::bad_alloc();
        }

        void* raw = static_cast<char*>(arena_base_) + current_offset_;
        current_offset_ += bytes;
        if (current_offset_ > peak_usage_) peak_usage_ = current_offset_;

        switch (fill_policy_) {
            case ArenaFillPolicy::Zero:
                std::memset(raw, 0, bytes);
                break;
            case ArenaFillPolicy::Poison:
                std::memset(raw, 0xCD, bytes);
                break;
            default:
                break;
        }

        return static_cast<T*>(raw);
    }

    // Scoped checkpoint support
    size_t save_point() const noexcept { return current_offset_; }
    void   restore(size_t p)  noexcept { current_offset_ = p; }

    size_t total_size()     const noexcept { return total_size_; }
    size_t used()           const noexcept { return current_offset_; }
    size_t peak_usage()     const noexcept { return peak_usage_; }
    size_t remaining()      const noexcept { return total_size_ - current_offset_; }
    const ArenaConfig& config() const noexcept { return config_; }

private:
    void allocate_arena(size_t bytes, const ArenaConfig& config) {
        if (bytes == 0) return;
#if defined(_MSC_VER)
        arena_base_ = _aligned_malloc(bytes, config.alignment);
#else
        if (posix_memalign(&arena_base_, config.alignment, bytes) != 0) {
            arena_base_ = nullptr;
        }
#endif
        if (arena_base_ == nullptr) throw std::bad_alloc();
        total_size_     = bytes;
        current_offset_ = 0;
    }

    void deallocate_arena() {
        if (arena_base_ == nullptr) return;
#if defined(_MSC_VER)
        _aligned_free(arena_base_);
#else
        std::free(arena_base_);
#endif
        arena_base_ = nullptr;
        total_size_ = 0;
        current_offset_ = 0;
    }

    void move_from(MemoryArena&& other) noexcept {
        config_          = other.config_;
        fill_policy_     = other.fill_policy_;
        arena_base_      = other.arena_base_;
        total_size_      = other.total_size_;
        current_offset_  = other.current_offset_;
        peak_usage_      = other.peak_usage_;
        other.arena_base_     = nullptr;
        other.total_size_     = 0;
        other.current_offset_ = 0;
        other.peak_usage_     = 0;
    }

    ArenaConfig     config_{};
    ArenaFillPolicy fill_policy_   = ArenaFillPolicy::None;
    void*           arena_base_    = nullptr;
    size_t          total_size_    = 0;
    size_t          current_offset_= 0;
    size_t          peak_usage_    = 0;
};

// =============================================================================
// Scoped checkpoint: restores arena on scope exit
// =============================================================================

class ArenaCheckpoint {
public:
    explicit ArenaCheckpoint(MemoryArena& a) noexcept
        : arena_(a), checkpoint_(a.save_point()) {}
    ~ArenaCheckpoint() { arena_.restore(checkpoint_); }
    ArenaCheckpoint(const ArenaCheckpoint&) = delete;
    ArenaCheckpoint& operator=(const ArenaCheckpoint&) = delete;
private:
    MemoryArena& arena_;
    size_t       checkpoint_;
};

} // namespace opg::common

#endif // OPG_COMMON_MEMORY_ARENA_HPP
