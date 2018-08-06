cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(libsuccient)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

set (CMAKE_CXX_FLAGS "-Wall -O3 -DNDEBUG -march=native -msse4.2 -mpopcnt")

add_subdirectory(src)
add_subdirectory(tests)
enable_testing()

add_test(NAME NaiveBitVectorTest COMMAND test_naive_bit_vector)
add_test(NAME BitVectorTest COMMAND test_bit_vector)

