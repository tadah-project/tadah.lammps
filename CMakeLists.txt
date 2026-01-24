cmake_minimum_required(VERSION 3.12)

project(Tadah.LAMMPS LANGUAGES CXX)

option(TADAH_MARCH_NATIVE "Build with -march=native" OFF)
option(TADAH_ENABLE_FAST_MATH "Enable fast math optimizations" OFF)
option(TADAH_LAMMPS_QUIET "Reduce configure output when used as a dependency" OFF)
option(TADAH_LAMMPS_QUIET "Reduce configure output when used as a dependency" OFF)

function(tadah_lammps_status)
  if(NOT TADAH_LAMMPS_QUIET)
    message(STATUS "${ARGV}")
  endif()
endfunction()

# True when this project is the top-level build (standalone)
set(TADAH_LAMMPS_STANDALONE OFF)
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(TADAH_LAMMPS_STANDALONE ON)
endif()

# 0. Guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR
    "In-source builds are not allowed.\n"
    "Create a separate build directory and run CMake from there.")
endif()

# 1. Global settings
if(TADAH_LAMMPS_STANDALONE AND NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# default build type for single-config generators (Ninja/Make)
if(TADAH_LAMMPS_STANDALONE)
  if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE)
  endif()
endif()

if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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
endif()

# 2. Bring in the shared compile/link options
include(FetchContent)
if(NOT TARGET tadah::compile_options)
  FetchContent_Declare (
    tadah_cmake
    GIT_REPOSITORY  https://git.ecdf.ed.ac.uk/tadah/cmake-tools.git
    GIT_TAG         main
  )
  FetchContent_MakeAvailable(tadah_cmake)   # defines tadah::compile_options
endif()

# (option) Developer builds -DTADAH_DEVELOP=ON
tadah_use_local_module(Tadah.MD MD)  # FetchContent name + folder name

# 3. Git branch information
tadah_get_git_branch(TADAH_GIT_BRANCH)
tadah_lammps_status("${PROJECT_NAME} is using TADAH_GIT_BRANCH:  ${TADAH_GIT_BRANCH}")

#4. Fetch Tadah! modules
FetchContent_Declare(
    Tadah.MD
    GIT_REPOSITORY https://git.ecdf.ed.ac.uk/tadah/md.git
    GIT_TAG ${TADAH_GIT_BRANCH}
)
FetchContent_MakeAvailable(Tadah.MD)

add_library(tadah.lammps STATIC src/lammps_tadah.cpp)
target_include_directories(tadah.lammps
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(tadah.lammps PUBLIC tadah::compile_options tadah.core tadah.models tadah.md)

if (NOT TADAH_BUILD_DOCS_ONLY)
  find_package(LAPACK)
  target_link_libraries (tadah.lammps PUBLIC ${LAPACK_LIBRARIES})
endif()

#########################################################################
if(TADAH_MARCH_NATIVE)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  else()
    tadah_lammps_status("${PROJECT_NAME} ...but -march=native is not supported")
  endif()
endif()
tadah_lammps_status("${PROJECT_NAME} Build with -march=native in ${TADAH_MARCH_NATIVE}")
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
tadah_lammps_status("${PROJECT_NAME}: Build with fast math optimisation is ${TADAH_ENABLE_FAST_MATH}")
#########################################################################

set(BUILD_TESTING OFF CACHE BOOL "Build tests")
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

tadah_lammps_status("===================================================================")
tadah_lammps_status("                                                                   ")
tadah_lammps_status("${PROJECT_NAME} in now configured and ready for compilation.       ")
tadah_lammps_status("                                                                   ")
tadah_lammps_status("--------------+----------------------------------------------------")
tadah_lammps_status("Command       |   Description                                      ")
tadah_lammps_status("--------------+----------------------------------------------------")
tadah_lammps_status("              |                                                    ")
tadah_lammps_status("make          |   Compile                                          ")
tadah_lammps_status("              |                                                    ")
tadah_lammps_status("make test     |   Run testing framework                            ")
tadah_lammps_status("              |                                                    ")
tadah_lammps_status("===================================================================")
