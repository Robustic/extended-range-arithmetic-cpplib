cmake_minimum_required(VERSION 3.14)
project(extended_range_arithmetic)

# GoogleTest requires at least C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  enable_testing()
endif()


ADD_LIBRARY(extended_range_arithmetic_module 
  ${PROJECT_SOURCE_DIR}/src/extended_range_arithmetic/FloatExp2Int64/FloatExp2Int64.cpp
  ${PROJECT_SOURCE_DIR}/src/extended_range_arithmetic/WideRangeNumber64/WideRangeNumber64.cpp
  ${PROJECT_SOURCE_DIR}/src/extended_range_arithmetic/Xnumber/Xnumber.cpp
  ${PROJECT_SOURCE_DIR}/src/extended_range_arithmetic/IntExp2Int64/IntExp2Int64.cpp
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_executable(
    FloatExp2Int64.test
    ${PROJECT_SOURCE_DIR}/src/extended_range_arithmetic/FloatExp2Int64/FloatExp2Int64.test.cpp
  )
  target_link_libraries(
    FloatExp2Int64.test
    extended_range_arithmetic_module
  )
  target_link_libraries(
    FloatExp2Int64.test
    GTest::gtest_main
  )
  include(GoogleTest)
  gtest_discover_tests(FloatExp2Int64.test)
endif()


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_executable(
    IntExp2Int64.test
    ${PROJECT_SOURCE_DIR}/src/extended_range_arithmetic/IntExp2Int64/IntExp2Int64.test.cpp
  )
  target_link_libraries(
    IntExp2Int64.test
    extended_range_arithmetic_module
  )
  target_link_libraries(
    IntExp2Int64.test
    GTest::gtest_main
  )
  include(GoogleTest)
  gtest_discover_tests(IntExp2Int64.test)
endif()


add_executable(
  Float64.try
  ${PROJECT_SOURCE_DIR}/src/Float64.try.cpp
)
target_link_libraries(
  Float64.try
  extended_range_arithmetic_module
)
