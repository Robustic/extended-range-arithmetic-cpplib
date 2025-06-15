cmake_minimum_required(VERSION 3.14)
project(FloatingExp2Integer)

# GoogleTest requires at least C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    URL C:/msys64/home/FloatingExp2Integer/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  enable_testing()
endif()


ADD_LIBRARY(floatingExp2IntegerModule 
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Float64PosExp2Int64/Float64PosExp2Int64.cpp
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/WideRangeNumber64/WideRangeNumber64.cpp
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Fukushima/Fukushima.cpp
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Int64PosExp2Int64/Int64PosExp2Int64.cpp
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
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
endif()


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_executable(
    Int64PosExp2Int64.test
    ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Int64PosExp2Int64/Int64PosExp2Int64.test.cpp
  )
  target_link_libraries(
    Int64PosExp2Int64.test
    floatingExp2IntegerModule
  )
  target_link_libraries(
    Int64PosExp2Int64.test
    GTest::gtest_main
  )
  include(GoogleTest)
  gtest_discover_tests(Int64PosExp2Int64.test)
endif()


add_executable(
  Float64.try
  ${PROJECT_SOURCE_DIR}/FloatingExp2Integer/Float64.try.cpp
)
target_link_libraries(
  Float64.try
  floatingExp2IntegerModule
)
