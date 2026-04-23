cmake_minimum_required(VERSION 3.20)
project(opg_v2 LANGUAGES CXX)

# C++20 required: class-type NTTPs, consteval, concepts, [[no_unique_address]].
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra -Wpedantic)

# -----------------------------------------------------------------------------
# opg_common — shared infrastructure (header-only)
#   - scalar types, PhysicsConfig + presets + concepts, optional_ptr/field,
#     memory arena, permutation primitives (ArrayRegistry, PHKeySorter),
#     generic box-leaf descriptor base.
#   - NO particle struct definitions live here. Each layout owns its own data
#     model in its own namespace (opg::layout_b, opg::layout_c,
#     opg::layout_b_coarse).
# -----------------------------------------------------------------------------
add_library(opg_common INTERFACE)
target_include_directories(opg_common INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/opg_common/include)
target_compile_features(opg_common INTERFACE cxx_std_20)

# -----------------------------------------------------------------------------
# opg_layout_b — B' (fine-grained cross-index) prototype (header-only)
#   Each access-pattern cluster is a separate small struct in its own array,
#   linked to type-specific arrays via explicit PLinkage. ~17 arrays.
# -----------------------------------------------------------------------------
add_library(opg_layout_b INTERFACE)
target_include_directories(opg_layout_b INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/opg_layout_b/include)
target_link_libraries(opg_layout_b INTERFACE opg_common)

# -----------------------------------------------------------------------------
# opg_layout_c — C (fine-grained positional mapping) prototype (header-only)
#   Same granularity as B'; connection via positional mapping through the
#   box-leaf descriptor (no stored cross-index).
# -----------------------------------------------------------------------------
add_library(opg_layout_c INTERFACE)
target_include_directories(opg_layout_c INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/opg_layout_c/include)
target_link_libraries(opg_layout_c INTERFACE opg_common)

# -----------------------------------------------------------------------------
# opg_layout_b_coarse — B (coarse cross-index) prototype (header-only)
#   3-way common split + one type-specific struct per type. Same connection
#   scheme as B' (PLinkageB cross-index). 7 arrays for FullMHD.
#   See opg_layout_b_coarse/README.md.
# -----------------------------------------------------------------------------
add_library(opg_layout_b_coarse INTERFACE)
target_include_directories(opg_layout_b_coarse INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/opg_layout_b_coarse/include)
target_link_libraries(opg_layout_b_coarse INTERFACE opg_common)

# -----------------------------------------------------------------------------
# Examples
# -----------------------------------------------------------------------------
add_executable(example_layout_b        examples/example_layout_b.cpp)
target_link_libraries(example_layout_b PRIVATE opg_layout_b)

add_executable(example_layout_c        examples/example_layout_c.cpp)
target_link_libraries(example_layout_c PRIVATE opg_layout_c)

add_executable(example_layout_b_coarse examples/example_layout_b_coarse.cpp)
target_link_libraries(example_layout_b_coarse PRIVATE opg_layout_b_coarse)

# Utility: prints sizeof(ParticleRecord<Cfg>) so distribute.c can be built
# with -DDATA_SIZE matching the layout-B exchange unit. See
# opg_layout_b_coarse/README.md for the wiring.
add_executable(print_record_size       examples/print_record_size.cpp)
target_link_libraries(print_record_size PRIVATE opg_layout_b_coarse)

message(STATUS "opg_v2 configured. Build with:  cmake --build .")
message(STATUS "Run examples:  ./example_layout_b  ./example_layout_c  ./example_layout_b_coarse")
