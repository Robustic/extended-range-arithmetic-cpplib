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

# Float64Exp2Int64 & Float64PosExp2Int64
ADD_LIBRARY(floatingExp2IntegerModule 
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Float64Exp2Int64/Float64Exp2Int64.cpp
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Float64PosExp2Int64/Float64PosExp2Int64.cpp
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Float32PosExp2Int32/Float32PosExp2Int32.cpp
)

add_executable(
  Float64Exp2Int64.try
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Float64Exp2Int64/Float64Exp2Int64.try.cpp
)
target_link_libraries(
  Float64Exp2Int64.try
  floatingExp2IntegerModule
)

add_executable(
  Float64Exp2Int64.test
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Float64Exp2Int64/Float64Exp2Int64.test.cpp
)
target_link_libraries(
  Float64Exp2Int64.test
  floatingExp2IntegerModule
)
target_link_libraries(
  Float64Exp2Int64.test
  GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(Float64Exp2Int64.test)


add_executable(
  Float64PosExp2Int64.try
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Float64PosExp2Int64/Float64PosExp2Int64.try.cpp
)
target_link_libraries(
  Float64PosExp2Int64.try
  floatingExp2IntegerModule
)

add_executable(
  Float64PosExp2Int64.test
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Float64PosExp2Int64/Float64PosExp2Int64.test.cpp
)
target_link_libraries(
  Float64PosExp2Int64.test
  floatingExp2IntegerModule
)
target_link_libraries(
  Float64PosExp2Int64.test
  GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(Float64PosExp2Int64.test)


add_executable(
  Float32PosExp2Int32.try
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Float32PosExp2Int32/Float32PosExp2Int32.try.cpp
)
target_link_libraries(
  Float32PosExp2Int32.try
  floatingExp2IntegerModule
)

add_executable(
  Float32PosExp2Int32.test
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Float32PosExp2Int32/Float32PosExp2Int32.test.cpp
)
target_link_libraries(
  Float32PosExp2Int32.test
  floatingExp2IntegerModule
)
target_link_libraries(
  Float32PosExp2Int32.test
  GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(Float32PosExp2Int32.test)
