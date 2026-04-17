message(STATUS "Will build opus with FARGAN")

# Locate autoreconf at CMake-configure time (when PATH is the user's shell PATH).
# The path is baked into the ExternalProject configure command so it works inside
# Xcode's sandboxed build environment, which does not inherit /opt/homebrew/bin.
find_program(AUTORECONF autoreconf HINTS /opt/homebrew/bin /usr/local/bin REQUIRED)
get_filename_component(AUTOTOOLS_BIN_DIR "${AUTORECONF}" DIRECTORY)

set(OPUS_CONFIGURE_FLAGS --with-pic --enable-osce --enable-dred --disable-shared --disable-doc --disable-extra-programs)

if(CMAKE_CROSSCOMPILING)
    list(APPEND OPUS_CONFIGURE_FLAGS --host=${CMAKE_C_COMPILER_TARGET} --target=${CMAKE_C_COMPILER_TARGET})
endif()

# Join flags into a single string for embedding in a sh -c command.
list(JOIN OPUS_CONFIGURE_FLAGS " " OPUS_CONFIGURE_FLAGS_STR)

# CONFIGURE_COMMAND wraps autogen+configure in a single sh -c call so that &&
# works as a shell operator and PATH is set for both commands.
set(CONFIGURE_COMMAND
    sh -c "PATH=${AUTOTOOLS_BIN_DIR}:/usr/bin:/bin ./autogen.sh && ./configure ${OPUS_CONFIGURE_FLAGS_STR}"
)

if(NOT DEFINED OPUS_URL)
    set(OPUS_URL https://github.com/xiph/opus/archive/940d4e5af64351ca8ba8390df3f555484c567fbb.zip)
endif()

include(ExternalProject)
if(APPLE AND BUILD_OSX_UNIVERSAL)
# Opus ./configure doesn't behave properly when built as a universal binary;
# build it twice and use lipo to create a universal libopus.a instead.
ExternalProject_Add(build_opus_x86
    DOWNLOAD_EXTRACT_TIMESTAMP NO
    BUILD_IN_SOURCE 1
    PATCH_COMMAND sh -c "patch dnn/nnet.h < ${CMAKE_SOURCE_DIR}/src/opus-nnet.h.diff"
    CONFIGURE_COMMAND sh -c "PATH=${AUTOTOOLS_BIN_DIR}:/usr/bin:/bin ./autogen.sh && ./configure --with-pic --enable-osce --enable-dred --disable-shared --disable-doc --disable-extra-programs --host=x86_64-apple-darwin --target=x86_64-apple-darwin CFLAGS='-arch x86_64 -O2 -mmacosx-version-min=10.11'"
    BUILD_COMMAND /usr/bin/make
    INSTALL_COMMAND ""
    URL ${OPUS_URL}
)
ExternalProject_Add(build_opus_arm
    DOWNLOAD_EXTRACT_TIMESTAMP NO
    BUILD_IN_SOURCE 1
    PATCH_COMMAND sh -c "patch dnn/nnet.h < ${CMAKE_SOURCE_DIR}/src/opus-nnet.h.diff"
    CONFIGURE_COMMAND sh -c "PATH=${AUTOTOOLS_BIN_DIR}:/usr/bin:/bin ./autogen.sh && ./configure --with-pic --enable-osce --enable-dred --disable-shared --disable-doc --disable-extra-programs --host=aarch64-apple-darwin --target=aarch64-apple-darwin CFLAGS='-arch arm64 -O2 -mmacosx-version-min=10.11'"
    BUILD_COMMAND /usr/bin/make
    INSTALL_COMMAND ""
    URL ${OPUS_URL}
)

ExternalProject_Get_Property(build_opus_arm BINARY_DIR)
ExternalProject_Get_Property(build_opus_arm SOURCE_DIR)
set(OPUS_ARM_BINARY_DIR ${BINARY_DIR})
ExternalProject_Get_Property(build_opus_x86 BINARY_DIR)
set(OPUS_X86_BINARY_DIR ${BINARY_DIR})

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libopus${CMAKE_STATIC_LIBRARY_SUFFIX}
    COMMAND lipo ${OPUS_ARM_BINARY_DIR}/.libs/libopus${CMAKE_STATIC_LIBRARY_SUFFIX} ${OPUS_X86_BINARY_DIR}/.libs/libopus${CMAKE_STATIC_LIBRARY_SUFFIX} -output ${CMAKE_CURRENT_BINARY_DIR}/libopus${CMAKE_STATIC_LIBRARY_SUFFIX} -create
    DEPENDS build_opus_arm build_opus_x86)

add_custom_target(
    libopus.a
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/libopus${CMAKE_STATIC_LIBRARY_SUFFIX})

include_directories(${SOURCE_DIR}/dnn ${SOURCE_DIR}/celt ${SOURCE_DIR}/include ${SOURCE_DIR})

add_library(opus STATIC IMPORTED)
add_dependencies(opus libopus.a)
set_target_properties(opus PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libopus${CMAKE_STATIC_LIBRARY_SUFFIX}"
)

else(APPLE AND BUILD_OSX_UNIVERSAL)
ExternalProject_Add(build_opus
    BUILD_IN_SOURCE 1
    PATCH_COMMAND sh -c "patch dnn/nnet.h < ${CMAKE_SOURCE_DIR}/src/opus-nnet.h.diff"
    CONFIGURE_COMMAND ${CONFIGURE_COMMAND}
    BUILD_COMMAND /usr/bin/make
    INSTALL_COMMAND ""
    URL ${OPUS_URL}
)

ExternalProject_Get_Property(build_opus BINARY_DIR)
ExternalProject_Get_Property(build_opus SOURCE_DIR)
add_library(opus STATIC IMPORTED)
add_dependencies(opus build_opus)

set_target_properties(opus PROPERTIES
    IMPORTED_LOCATION "${BINARY_DIR}/.libs/libopus${CMAKE_STATIC_LIBRARY_SUFFIX}"
    IMPORTED_IMPLIB   "${BINARY_DIR}/.libs/libopus${CMAKE_STATIC_LIBRARY_SUFFIX}"
)

include_directories(${SOURCE_DIR}/dnn ${SOURCE_DIR}/celt ${SOURCE_DIR}/include ${SOURCE_DIR})
endif(APPLE AND BUILD_OSX_UNIVERSAL)
