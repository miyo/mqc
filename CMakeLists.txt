# CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(qc_sim LANGUAGES CXX)

# Use C++20 (for <numbers> etc.)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Warnings
if (MSVC)
  add_compile_options(/W4 /permissive-)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# ---- Option A: simple one-target build (quickest) ----
add_executable(qc_sim
  sources/main.cc
  sources/gate.cc
  sources/qc.cc
  sources/utils.cc
  sources/surface_code.cc
)

add_executable(qc_surface
  sources/main_surface.cc
  sources/gate.cc
  sources/qc.cc
  sources/utils.cc
  sources/surface_code.cc
)
target_include_directories(qc_surface PRIVATE sources)

target_include_directories(qc_sim PRIVATE sources)

# For Google Test
include(FetchContent)
include(CTest)  # enables BUILD_TESTING

if (BUILD_TESTING)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
  )
  # For older CMake: avoid overriding parent compiler flags
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  add_executable(qc_tests
    tests/core_ops_test.cc
    tests/measure_and_ctrl_test.cc
    tests/surface_test.cc
    sources/gate.cc
    sources/qc.cc
    sources/utils.cc
    sources/surface_code.cc
  )
  target_include_directories(qc_tests PRIVATE sources)
  target_link_libraries(qc_tests
    GTest::gtest_main
  )

  include(GoogleTest)
  gtest_discover_tests(qc_tests)
endif()
