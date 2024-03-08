cmake_minimum_required(VERSION 3.12)

project(tadah.LAMMPS)
set(TADAH "Ta-dah!")

# Here are the options... set them with -D,
# for example:
# cmake .. -DTADAH_BUILD_TESTS=ON
option(TADAH_MARCH_NATIVE "Build with -march=native" OFF)
option(TADAH_ENABLE_FAST_MATH "Enable fast math optimizations" OFF)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt and CMakeFiles dir. ")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(BUILD_SHARED_LIBS)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

if (BUILD_SHARED_LIBS)
    # Always full RPATH
    # use, i.e. don't skip the full RPATH for the build tree
    set(CMAKE_SKIP_BUILD_RPATH FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    # the RPATH to be used when installing, but only if it's not a system directory
    list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    if("${isSystemDir}" STREQUAL "-1")
        set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    endif("${isSystemDir}" STREQUAL "-1")
endif()

set(CMAKE_CXX_FLAGS_DEBUG
    "-g -O2 -fsanitize=address -Wall -Wpedantic -pedantic-errors -Wextra"
    CACHE STRING
    "Debug build."
    FORCE)

set(CMAKE_CXX_FLAGS_RELEASE
    "-g -O3 -Wall -Wpedantic -Wextra"
    CACHE STRING
    "Release build."
    FORCE)

set(CMAKE_CXX_FLAGS_PROFILE
    "-g -O3 -Wall -Wpedantic -Wextra -pg -fno-inline "
    CACHE STRING
    "Flags used by the C++ compiler during profile builds."
    FORCE)

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Check is build supported
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)
if (CMAKE_BUILD_TYPE STREQUAL "DEBUG" 
        OR
        CMAKE_BUILD_TYPE STREQUAL "RELEASE"
        OR
        CMAKE_BUILD_TYPE STREQUAL "PROFILE"
        )
    message(STATUS "${TADAH} Build type:  ${CMAKE_BUILD_TYPE}.")
else()
    message(FATAL_ERROR "
    Selected build type is not supported: ${CMAKE_BUILD_TYPE}.
    Supported build types: Release, Debug, Profile"
    )
endif()

# Update the documentation string of CMAKE_BUILD_TYPE for GUIs
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release Profile."
    FORCE)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)

FetchContent_Declare(
    Tadah.MD
    GIT_REPOSITORY git@git.ecdf.ed.ac.uk:tadah/md.git
    GIT_TAG origin/main
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/MD
    PATCH_COMMAND  mkdir -p build
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/MD/build
    )
FetchContent_MakeAvailable(Tadah.MD)

# Set CONFIG_FILE_DIR before adding CORE module
# otherwise CORE module will set it itself to src dir
# set(CONFIG_FILE_DIR "${CMAKE_INSTALL_PREFIX}/include/tadah/config/")

add_library(tadah.lammps
    $<TARGET_OBJECTS:tadah.core.obj>
    $<TARGET_OBJECTS:tadah.models.obj>
    $<TARGET_OBJECTS:tadah.md.obj>
    )
target_link_libraries(tadah.lammps PRIVATE config::rc)

#########################################################################
if(TADAH_MARCH_NATIVE)
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
    if(COMPILER_SUPPORTS_MARCH_NATIVE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    else()
        message(STATUS "${TADAH} ...but -march=native is not supported")
    endif()
endif()
message(STATUS "${TADAH} Build with -march=native in ${TADAH_MARCH_NATIVE}")
#########################################################################

#########################################################################
if(TADAH_ENABLE_FAST_MATH)
    if((CMAKE_CXX_COMPILER_ID MATCHES "Clang") OR (CMAKE_CXX_COMPILER_ID MATCHES "GNU"))
        add_compile_options(-ffast-math)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        add_compile_options(/fp:fast)
    endif()
else()
    if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        add_compile_options(/fp:precise)
    endif()
endif()
message(STATUS "${TADAH}: Build with fast math optimisation is ${TADAH_ENABLE_FAST_MATH}")
#########################################################################

#########################################################################

if(TADAH_BUILD_TESTS)
    include(CTest) 
    add_subdirectory(tests)
endif()

message(STATUS "===================================================================")
message(STATUS "                                                                   ")
message(STATUS "${TADAH}.MD in now configured and ready for compilation.           ")
message(STATUS "                                                                   ")
message(STATUS "--------------+----------------------------------------------------")
message(STATUS "Command       |   Description                                      ")
message(STATUS "--------------+----------------------------------------------------")
message(STATUS "              |                                                    ")
message(STATUS "make          |   Compile                                          ")
message(STATUS "              |                                                    ")
message(STATUS "make test     |   Run testing framework                            ")
message(STATUS "              |                                                    ")
message(STATUS "===================================================================")
