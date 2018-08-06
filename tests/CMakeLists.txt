cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

set(GTEST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest/googletest/)
include_directories(${GTEST_ROOT}/include ${CMAKE_CURRENT_SOURCE_DIR}/../include)
link_directories(${GTEST_ROOT})
add_subdirectory(${GTEST_ROOT})

add_executable(test_naive_bit_vector
  ${CMAKE_CURRENT_SOURCE_DIR}/test_naive_bit_vector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/naive_bit_vector.cc)
target_link_libraries(test_naive_bit_vector gtest gtest_main pthread)

add_executable(test_bit_vector
  ${CMAKE_CURRENT_SOURCE_DIR}/test_bit_vector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/bit_vector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/naive_bit_vector.cc)
target_link_libraries(test_bit_vector gtest gtest_main pthread)
