cmake_minimum_required(VERSION 3.16)
project(radae_nopy C)

set(CMAKE_C_STANDARD 11)

# Build opus with FARGAN/LPCNet support
option(AVX "Enable AVX CPU optimizations." ON)
include(cmake/BuildOpus.cmake)

add_subdirectory(src)

enable_testing()
add_subdirectory(test)
