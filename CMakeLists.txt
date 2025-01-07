cmake_minimum_required(VERSION 3.14)
project(FloatingExp2Integer)

# GoogleTest requires at least C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/f3c355f9dd382bc2c323be2713e351a578b68c61.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


enable_testing()

add_executable(
  DoubleExp2Int.test
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/DoubleExp2Int.test.cpp
)
target_link_libraries(
  DoubleExp2Int.test
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(DoubleExp2Int.test)
