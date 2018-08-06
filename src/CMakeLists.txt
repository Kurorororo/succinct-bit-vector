cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)
add_library(succient_bv STATIC bit_vector.cc naive_bit_vector.cc)
