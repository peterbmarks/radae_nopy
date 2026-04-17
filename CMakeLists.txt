cmake_minimum_required(VERSION 3.16)
cmake_policy(SET CMP0114 NEW)
cmake_policy(SET CMP0135 NEW)
project(radae_nopy C)

set(CMAKE_C_STANDARD 11)

# Build opus with FARGAN/LPCNet support
include(cmake/BuildOpus.cmake)

add_subdirectory(src)
