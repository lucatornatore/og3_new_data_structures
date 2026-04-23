/**
 * @file memory_arena.hpp
 * @brief Single-buffer arena with persistent allocations, recyclable regions,
 *        and a top-of-buffer temporary stack.
 *
 * Design:
 *   - One aligned allocation at arena construction.
 *   - Persistent allocations grow upward from the bottom of the buffer.
 *   - Temporary allocations grow downward from the top of the buffer.
 *   - Persistent regions may be returned to the arena and reused later.
 *   - Temporary regions are reclaimed wholesale with stack marks.
 *
 * This gives the resident layer a stable single-buffer ownership model while
 * still supporting two common HPC needs:
 *   1. Long-lived recyclable regions (auxiliary buffers, staging ranges,
 *      scratch pools that outlive one scope).
 *   2. Cheap throw-away temporaries for per-rebuild / per-kernel work.
 */

#ifndef OPG_COMMON_MEMORY_ARENA_HPP
#define OPG_COMMON_MEMORY_ARENA_HPP

#include "../types/scalar_types.hpp"
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace opg::common {

// =============================================================================
// ArenaConfig
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
// Fill policy
// =============================================================================

enum class ArenaFillPolicy {
    None,
    Zero,
    Poison
};

// =============================================================================
// Region descriptors
// =============================================================================

struct ArenaRegion {
    void*  ptr       = nullptr;
    size_t offset    = 0;
    size_t bytes     = 0;
    size_t alignment = 1;

    constexpr explicit operator bool() const noexcept { return ptr != nullptr; }
};

template<typename T>
struct ArenaArray {
    T*     ptr       = nullptr;
    size_t offset    = 0;
    count_t count    = 0;
    size_t bytes     = 0;
    size_t alignment = alignof(T);

    constexpr explicit operator bool() const noexcept { return ptr != nullptr; }
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

    // -------------------------------------------------------------------------
    // Persistent allocations (bottom-up)
    // -------------------------------------------------------------------------

    template<ParticleComponent T>
    T* allocate(count_t n, const char* name = "") {
        return allocate_region<T>(n, name).ptr;
    }

    template<ParticleComponent T>
    ArenaArray<T> allocate_region(count_t n, const char* name = "") {
        const size_t bytes = static_cast<size_t>(n) * sizeof(T);
        const size_t align = max_alignment(alignof(T), config_.alignment);
        const ArenaRegion region = allocate_bytes(bytes, align, name);
        return ArenaArray<T>{
            static_cast<T*>(region.ptr),
            region.offset,
            n,
            region.bytes,
            region.alignment
        };
    }

    ArenaRegion allocate_bytes(size_t bytes,
                               size_t alignment,
                               const char* /*name*/ = "")
    {
        if (bytes == 0) return {};

        const size_t align = max_alignment(alignment, config_.alignment);
        assert(is_power_of_two(align));

        // 1. Reuse a previously freed persistent block if possible.
        for (size_t i = 0; i < free_count_; ++i) {
            const size_t block_begin = free_blocks_[i].offset;
            const size_t block_end   = free_blocks_[i].offset + free_blocks_[i].bytes;
            const size_t aligned     = align_up(block_begin, align);

            if (aligned > block_end) continue;
            if (aligned + bytes > block_end) continue;
            if (aligned + bytes > temp_begin_) continue;

            const size_t prefix = aligned - block_begin;
            const size_t suffix = block_end - (aligned + bytes);

            if (prefix > 0 && suffix > 0) {
                assert(free_count_ < MaxFreeBlocks);
                if (free_count_ >= MaxFreeBlocks) throw std::bad_alloc();
                for (size_t j = free_count_; j > i + 1; --j) {
                    free_blocks_[j] = free_blocks_[j - 1];
                }
                free_blocks_[i]     = FreeBlock{block_begin, prefix};
                free_blocks_[i + 1] = FreeBlock{aligned + bytes, suffix};
                ++free_count_;
            } else if (prefix > 0) {
                free_blocks_[i] = FreeBlock{block_begin, prefix};
            } else if (suffix > 0) {
                free_blocks_[i] = FreeBlock{aligned + bytes, suffix};
            } else {
                erase_free_block(i);
            }

            persistent_live_bytes_ += bytes;
            apply_fill(static_cast<char*>(arena_base_) + aligned, bytes);
            update_peak_usage();
            return ArenaRegion{
                static_cast<char*>(arena_base_) + aligned,
                aligned,
                bytes,
                align
            };
        }

        // 2. Fall back to the bump pointer.
        const size_t aligned = align_up(persistent_end_, align);
        if (aligned + bytes > temp_begin_) throw std::bad_alloc();

        persistent_end_ = aligned + bytes;
        persistent_live_bytes_ += bytes;
        apply_fill(static_cast<char*>(arena_base_) + aligned, bytes);
        update_peak_usage();

        return ArenaRegion{
            static_cast<char*>(arena_base_) + aligned,
            aligned,
            bytes,
            align
        };
    }

    void deallocate(const ArenaRegion& region) noexcept {
        if (!region || region.bytes == 0) return;
        assert(region.offset + region.bytes <= persistent_end_);
        if (region.offset + region.bytes > persistent_end_) return;

        persistent_live_bytes_ -= region.bytes;
        insert_free_block(region.offset, region.bytes);
        release_trailing_free_blocks();
    }

    template<typename T>
    void deallocate(const ArenaArray<T>& region) noexcept {
        deallocate(ArenaRegion{region.ptr, region.offset, region.bytes, region.alignment});
    }

    // -------------------------------------------------------------------------
    // Temporary stack allocations (top-down)
    // -------------------------------------------------------------------------

    template<ParticleComponent T>
    T* allocate_temp(count_t n, const char* name = "") {
        return allocate_temp_region<T>(n, name).ptr;
    }

    template<ParticleComponent T>
    ArenaArray<T> allocate_temp_region(count_t n, const char* name = "") {
        const size_t bytes = static_cast<size_t>(n) * sizeof(T);
        const size_t align = max_alignment(alignof(T), config_.alignment);
        const ArenaRegion region = allocate_temp_bytes(bytes, align, name);
        return ArenaArray<T>{
            static_cast<T*>(region.ptr),
            region.offset,
            n,
            region.bytes,
            region.alignment
        };
    }

    ArenaRegion allocate_temp_bytes(size_t bytes,
                                    size_t alignment,
                                    const char* /*name*/ = "")
    {
        if (bytes == 0) return {};

        const size_t align = max_alignment(alignment, config_.alignment);
        assert(is_power_of_two(align));
        const size_t tentative_begin = temp_begin_ - bytes;
        const size_t aligned_begin   = align_down(tentative_begin, align);

        if (aligned_begin < persistent_end_) throw std::bad_alloc();

        temp_begin_ = aligned_begin;
        apply_fill(static_cast<char*>(arena_base_) + temp_begin_, bytes);
        update_peak_usage();

        return ArenaRegion{
            static_cast<char*>(arena_base_) + temp_begin_,
            temp_begin_,
            bytes,
            align
        };
    }

    size_t stack_save_point() const noexcept { return temp_begin_; }
    void   stack_restore(size_t p) noexcept {
        assert(p >= persistent_end_ && p <= total_size_);
        if (p < persistent_end_ || p > total_size_) return;
        temp_begin_ = p;
    }

    // Backward-compatible aliases: checkpoints now refer to the temp stack.
    size_t save_point() const noexcept { return stack_save_point(); }
    void   restore(size_t p) noexcept  { stack_restore(p); }

    // -------------------------------------------------------------------------
    // Introspection
    // -------------------------------------------------------------------------

    size_t total_size() const noexcept { return total_size_; }
    size_t used() const noexcept { return persistent_live_bytes_ + temp_bytes_in_use(); }
    size_t peak_usage() const noexcept { return peak_usage_; }
    size_t remaining() const noexcept { return total_size_ - used(); }
    size_t persistent_reserved_bytes() const noexcept { return persistent_end_; }
    size_t persistent_live_bytes() const noexcept { return persistent_live_bytes_; }
    size_t temp_bytes_in_use() const noexcept { return total_size_ - temp_begin_; }
    const ArenaConfig& config() const noexcept { return config_; }

private:
    static constexpr size_t MaxFreeBlocks = 256;

    struct FreeBlock {
        size_t offset = 0;
        size_t bytes  = 0;
    };

    static constexpr size_t max_alignment(size_t a, size_t b) noexcept {
        return a > b ? a : b;
    }

    static constexpr bool is_power_of_two(size_t x) noexcept {
        return x != 0 && (x & (x - 1u)) == 0;
    }

    static constexpr size_t align_up(size_t value, size_t alignment) noexcept {
        return (value + alignment - 1u) & ~(alignment - 1u);
    }

    static constexpr size_t align_down(size_t value, size_t alignment) noexcept {
        return value & ~(alignment - 1u);
    }

    void apply_fill(void* ptr, size_t bytes) noexcept {
        switch (fill_policy_) {
            case ArenaFillPolicy::Zero:
                std::memset(ptr, 0, bytes);
                break;
            case ArenaFillPolicy::Poison:
                std::memset(ptr, 0xCD, bytes);
                break;
            default:
                break;
        }
    }

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
        total_size_            = bytes;
        persistent_end_        = 0;
        persistent_live_bytes_ = 0;
        temp_begin_            = bytes;
        free_count_            = 0;
    }

    void deallocate_arena() {
        if (arena_base_ == nullptr) return;
#if defined(_MSC_VER)
        _aligned_free(arena_base_);
#else
        std::free(arena_base_);
#endif
        arena_base_            = nullptr;
        total_size_            = 0;
        persistent_end_        = 0;
        persistent_live_bytes_ = 0;
        temp_begin_            = 0;
        peak_usage_            = 0;
        free_count_            = 0;
    }

    void move_from(MemoryArena&& other) noexcept {
        config_                = other.config_;
        fill_policy_           = other.fill_policy_;
        arena_base_            = other.arena_base_;
        total_size_            = other.total_size_;
        persistent_end_        = other.persistent_end_;
        persistent_live_bytes_ = other.persistent_live_bytes_;
        temp_begin_            = other.temp_begin_;
        peak_usage_            = other.peak_usage_;
        free_count_            = other.free_count_;
        for (size_t i = 0; i < other.free_count_; ++i) {
            free_blocks_[i] = other.free_blocks_[i];
        }

        other.arena_base_            = nullptr;
        other.total_size_            = 0;
        other.persistent_end_        = 0;
        other.persistent_live_bytes_ = 0;
        other.temp_begin_            = 0;
        other.peak_usage_            = 0;
        other.free_count_            = 0;
    }

    void erase_free_block(size_t idx) noexcept {
        assert(idx < free_count_);
        for (size_t i = idx + 1; i < free_count_; ++i) {
            free_blocks_[i - 1] = free_blocks_[i];
        }
        --free_count_;
    }

    void insert_free_block(size_t offset, size_t bytes) noexcept {
        if (bytes == 0) return;
        assert(free_count_ < MaxFreeBlocks);
        if (free_count_ >= MaxFreeBlocks) return;

        size_t pos = 0;
        while (pos < free_count_ && free_blocks_[pos].offset < offset) ++pos;
        for (size_t i = free_count_; i > pos; --i) {
            free_blocks_[i] = free_blocks_[i - 1];
        }
        free_blocks_[pos] = FreeBlock{offset, bytes};
        ++free_count_;
        coalesce_free_blocks();
    }

    void coalesce_free_blocks() noexcept {
        if (free_count_ < 2) return;
        size_t out = 0;
        for (size_t i = 1; i < free_count_; ++i) {
            const size_t current_end = free_blocks_[out].offset + free_blocks_[out].bytes;
            if (current_end >= free_blocks_[i].offset) {
                const size_t next_end = free_blocks_[i].offset + free_blocks_[i].bytes;
                if (next_end > current_end) {
                    free_blocks_[out].bytes = next_end - free_blocks_[out].offset;
                }
            } else {
                ++out;
                free_blocks_[out] = free_blocks_[i];
            }
        }
        free_count_ = out + 1;
    }

    void release_trailing_free_blocks() noexcept {
        bool changed = true;
        while (changed) {
            changed = false;
            for (size_t i = 0; i < free_count_; ++i) {
                const size_t block_end = free_blocks_[i].offset + free_blocks_[i].bytes;
                if (block_end == persistent_end_) {
                    persistent_end_ = free_blocks_[i].offset;
                    erase_free_block(i);
                    changed = true;
                    break;
                }
            }
        }
    }

    void update_peak_usage() noexcept {
        const size_t current_used = used();
        if (current_used > peak_usage_) peak_usage_ = current_used;
    }

    ArenaConfig     config_{};
    ArenaFillPolicy fill_policy_ = ArenaFillPolicy::None;
    void*           arena_base_ = nullptr;
    size_t          total_size_ = 0;

    size_t persistent_end_ = 0;
    size_t persistent_live_bytes_ = 0;
    size_t temp_begin_ = 0;
    size_t peak_usage_ = 0;

    FreeBlock free_blocks_[MaxFreeBlocks]{};
    size_t    free_count_ = 0;
};

// =============================================================================
// Scoped temp-stack checkpoint
// =============================================================================

class ArenaCheckpoint {
public:
    explicit ArenaCheckpoint(MemoryArena& a) noexcept
        : arena_(a), checkpoint_(a.stack_save_point()) {}
    ~ArenaCheckpoint() { arena_.stack_restore(checkpoint_); }

    ArenaCheckpoint(const ArenaCheckpoint&) = delete;
    ArenaCheckpoint& operator=(const ArenaCheckpoint&) = delete;

private:
    MemoryArena& arena_;
    size_t       checkpoint_;
};

} // namespace opg::common

#endif // OPG_COMMON_MEMORY_ARENA_HPP
