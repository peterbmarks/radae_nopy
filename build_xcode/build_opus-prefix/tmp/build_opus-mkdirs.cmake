# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/Volumes/LaCie2TB/Developer/radae_nopy/build_xcode/build_opus-prefix/src/build_opus")
  file(MAKE_DIRECTORY "/Volumes/LaCie2TB/Developer/radae_nopy/build_xcode/build_opus-prefix/src/build_opus")
endif()
file(MAKE_DIRECTORY
  "/Volumes/LaCie2TB/Developer/radae_nopy/build_xcode/build_opus-prefix/src/build_opus-build"
  "/Volumes/LaCie2TB/Developer/radae_nopy/build_xcode/build_opus-prefix"
  "/Volumes/LaCie2TB/Developer/radae_nopy/build_xcode/build_opus-prefix/tmp"
  "/Volumes/LaCie2TB/Developer/radae_nopy/build_xcode/build_opus-prefix/src/build_opus-stamp"
  "/Volumes/LaCie2TB/Developer/radae_nopy/build_xcode/build_opus-prefix/src"
  "/Volumes/LaCie2TB/Developer/radae_nopy/build_xcode/build_opus-prefix/src/build_opus-stamp"
)

set(configSubDirs Debug;Release;MinSizeRel;RelWithDebInfo)
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Volumes/LaCie2TB/Developer/radae_nopy/build_xcode/build_opus-prefix/src/build_opus-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Volumes/LaCie2TB/Developer/radae_nopy/build_xcode/build_opus-prefix/src/build_opus-stamp${cfgdir}") # cfgdir has leading slash
endif()
