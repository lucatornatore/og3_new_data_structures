#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"

cmake -S . -B build >/dev/null
cmake --build build -j4 >/dev/null

./build/example_layout_b >/dev/null
./build/example_layout_c >/dev/null
./build/example_layout_b_coarse >/dev/null

INC='-I. -Iopg_common/include -Iopg_layout_b/include -Iopg_layout_c/include -Iopg_layout_b_coarse/include'

cat > /tmp/opg_compile_fail_gas_core.cpp <<'CPP'
#include <opg_common/types/physics_config.hpp>
#include <opg_common/memory/memory_arena.hpp>
#include <opg_layout_b/particle_container.hpp>
using namespace opg::common;
using namespace opg::layout_b;
int main() {
    MemoryArena arena(1024*1024, ArenaConfig{});
    Capacity cap{16,0,0,0};
    ParticleContainer<physics_configs::GravityOnly> c(arena, cap);
    (void)c.gas_core();
    return 0;
}
CPP
if g++ -std=c++20 -Wall -Wextra -Wpedantic $INC /tmp/opg_compile_fail_gas_core.cpp -c -o /tmp/opg_compile_fail_gas_core.o >/dev/null 2>/tmp/opg_compile_fail_gas_core.err; then
    echo '[FAIL] gas_core compiled for GravityOnly' >&2
    exit 1
fi

echo "[OK] gas_core absent for GravityOnly (B')"

cat > /tmp/opg_compile_fail_gas_all.cpp <<'CPP'
#include <opg_common/types/physics_config.hpp>
#include <opg_common/memory/memory_arena.hpp>
#include <opg_layout_b_coarse/particle_container.hpp>
using namespace opg::common;
using namespace opg::layout_b_coarse;
int main() {
    MemoryArena arena(1024*1024, ArenaConfig{});
    Capacity cap{16,0,0,0};
    ParticleContainer<physics_configs::GravityOnly> c(arena, cap);
    (void)c.gas_all();
    return 0;
}
CPP
if g++ -std=c++20 -Wall -Wextra -Wpedantic $INC /tmp/opg_compile_fail_gas_all.cpp -c -o /tmp/opg_compile_fail_gas_all.o >/dev/null 2>/tmp/opg_compile_fail_gas_all.err; then
    echo '[FAIL] gas_all compiled for GravityOnly' >&2
    exit 1
fi

echo '[OK] gas_all absent for GravityOnly (B coarse)'

cat > /tmp/opg_check_nmet0.cpp <<'CPP'
#include <iostream>
#include <opg_layout_b/particles.hpp>
#include <opg_layout_c/particles.hpp>
#include <opg_layout_b_coarse/particles.hpp>
int main() {
    std::cout << sizeof(opg::layout_b::StarCore<0>) << ' ' << sizeof(opg::layout_b::StarCore<11>) << '\n';
    std::cout << sizeof(opg::layout_c::StarCore<0>) << ' ' << sizeof(opg::layout_c::StarCore<11>) << '\n';
    std::cout << sizeof(opg::layout_b_coarse::StarCore<0>) << ' ' << sizeof(opg::layout_b_coarse::StarCore<11>) << '\n';
    std::cout << sizeof(opg::layout_b::GasMetal<0>) << ' ' << sizeof(opg::layout_b::GasMetal<11>) << '\n';
    return 0;
}
CPP
g++ -std=c++20 -Wall -Wextra -Wpedantic $INC /tmp/opg_check_nmet0.cpp -o /tmp/opg_check_nmet0
NMET_OUT=$(/tmp/opg_check_nmet0)
echo "$NMET_OUT" | awk '{ if ($1 >= $2) exit 1; }'
echo '[OK] NMET=0 specialisations are smaller than NMET=11'

cat > /tmp/opg_check_arena.cpp <<'CPP'
#include <cassert>
#include <opg_common/memory/memory_arena.hpp>
using namespace opg::common;
struct alignas(64) X { int v; };
int main() {
    MemoryArena arena(4096, ArenaConfig{});
    auto a = arena.allocate_region<X>(4, "a");
    auto b = arena.allocate_region<X>(2, "b");
    auto before_free = arena.persistent_live_bytes();
    arena.deallocate(b);
    auto c = arena.allocate_region<X>(2, "c");
    assert(c.offset == b.offset);
    {
        ArenaCheckpoint ck(arena);
        auto t1 = arena.allocate_temp_region<X>(8, "t1");
        auto t2 = arena.allocate_temp_region<X>(4, "t2");
        (void)t1;
        (void)t2;
        assert(arena.temp_bytes_in_use() > 0);
    }
    assert(arena.temp_bytes_in_use() == 0);
    assert(arena.persistent_live_bytes() == before_free);
    (void)a;
    (void)c;
    return 0;
}
CPP
g++ -std=c++20 -Wall -Wextra -Wpedantic $INC /tmp/opg_check_arena.cpp -o /tmp/opg_check_arena
/tmp/opg_check_arena

echo '[OK] arena free-list reuse and temp stack restore'
echo '[OK] all regression checks passed'
