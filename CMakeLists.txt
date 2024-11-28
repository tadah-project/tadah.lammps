cmake_minimum_required(VERSION 3.12)

project(tadah.LAMMPS)
set(TADAH "Tadah!")

option(TADAH_MARCH_NATIVE "Build with -march=native" OFF)
option(TADAH_ENABLE_FAST_MATH "Enable fast math optimizations" OFF)

# Guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt and CMakeFiles dir. ")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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

# Required for GitLab CI/CD
if(DEFINED ENV{CI_COMMIT_REF_NAME})
  set(GIT_BRANCH "$ENV{CI_COMMIT_REF_NAME}")
else()
  execute_process(
        COMMAND git rev-parse --symbolic-full-name --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
message(STATUS "${TADAH}.LAMMPS Using GIT_BRANCH:  ${GIT_BRANCH}.")

include(FetchContent)
FetchContent_Declare(
    Tadah.MD
    GIT_REPOSITORY https://git.ecdf.ed.ac.uk/tadah/md.git
    GIT_TAG origin/${GIT_BRANCH}
)
FetchContent_MakeAvailable(Tadah.MD)

add_library(tadah.lammps STATIC lammps_tadah.cpp)
target_link_libraries(tadah.lammps PUBLIC tadah.md)
target_link_libraries(tadah.lammps PUBLIC tadah.core)

find_package(LAPACK)
target_link_libraries (tadah.lammps PUBLIC ${LAPACK_LIBRARIES})

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
