cmake_minimum_required(VERSION 3.20)

option(PKG_ML-TADAH "Enable Tadah ML-TADAH integration (external)" OFF)

set(TADAH_LAMMPS_SOURCE_DIR
    ""
    CACHE PATH "Path to a local Tadah.LAMMPS checkout")

set(TADAH_LAMMPS_QUIET
    ON
    CACHE BOOL "Minimize Tadah.LAMMPS configure output when built as dependency")

set(TADAH_LAMMPS_PAIR_STYLES
    "tadah:PairTadah"
    CACHE STRING "Semicolon-separated list of pair styles style:Class[:header]")

# Capture the directory of THIS script now (before any DEFER happens).
# This must be stable and must NOT rely on CMAKE_CURRENT_LIST_DIR later.
set(_TADAH_INTEGRATION_DIR "${CMAKE_CURRENT_LIST_DIR}")

function(_tadah_require_dir varname)
  if(NOT ${varname})
    message(FATAL_ERROR "ML-TADAH: ${varname} is not set")
  endif()
  if(NOT IS_DIRECTORY "${${varname}}")
    message(FATAL_ERROR "ML-TADAH: ${varname} is not a directory: ${${varname}}")
  endif()
endfunction()

function(_tadah_patch_style_pair_h header_name)
  if(NOT DEFINED LAMMPS_STYLE_HEADERS_DIR)
    message(FATAL_ERROR "ML-TADAH: LAMMPS_STYLE_HEADERS_DIR is not defined - cannot patch style headers")
  endif()

  set(style_file "${LAMMPS_STYLE_HEADERS_DIR}/style_pair.h")
  if(NOT EXISTS "${style_file}")
    message(FATAL_ERROR "ML-TADAH: Cannot find generated file: ${style_file}")
  endif()

  set(include_line "#include \"${header_name}\"")

  file(READ "${style_file}" contents)
  string(FIND "${contents}" "${include_line}" inc_pos)
  if(inc_pos EQUAL -1)
    file(APPEND "${style_file}" "\n${include_line}\n")
    message(STATUS "ML-TADAH: Patched style_pair.h: added ${include_line}")
  endif()
endfunction()

function(_tadah_apply_after_lammps)
  if(NOT PKG_ML-TADAH)
    return()
  endif()

  if(NOT TARGET lammps)
    message(FATAL_ERROR "ML-TADAH: LAMMPS target 'lammps' not found (unexpected)")
  endif()

  if(NOT TARGET tadah.lammps)
    message(FATAL_ERROR "ML-TADAH: target 'tadah.lammps' not found (Tadah.LAMMPS did not configure?)")
  endif()

  # Use the saved integration directory, not CMAKE_CURRENT_LIST_DIR
  if(NOT EXISTS "${_TADAH_INTEGRATION_DIR}/ML-TADAH.cmake")
    message(FATAL_ERROR "ML-TADAH: Cannot find ${_TADAH_INTEGRATION_DIR}/ML-TADAH.cmake")
  endif()
  include("${_TADAH_INTEGRATION_DIR}/ML-TADAH.cmake")

  foreach(spec IN LISTS TADAH_LAMMPS_PAIR_STYLES)
    string(REPLACE ":" ";" parts "${spec}")
    list(LENGTH parts n)
    if(n LESS 2 OR n GREATER 3)
      message(FATAL_ERROR "ML-TADAH: Bad entry '${spec}' (expected style:Class[:header])")
    endif()

    list(GET parts 0 style_name)
    list(GET parts 1 class_name)

    if(n EQUAL 3)
      list(GET parts 2 header_name)
    else()
      set(header_name "pair_${style_name}.h")
    endif()

    _tadah_patch_style_pair_h("${header_name}")
  endforeach()
endfunction()

get_property(_tadah_done GLOBAL PROPERTY TADAH_ML_TADAH_INTEGRATION_DONE)
if(_tadah_done)
  return()
endif()
set_property(GLOBAL PROPERTY TADAH_ML_TADAH_INTEGRATION_DONE TRUE)

if(PKG_ML-TADAH)
  _tadah_require_dir(TADAH_LAMMPS_SOURCE_DIR)

  set(TADAH_LAMMPS_QUIET "${TADAH_LAMMPS_QUIET}" CACHE BOOL "" FORCE)

  add_subdirectory("${TADAH_LAMMPS_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/tadah_lammps_build")

  set(Tadah_LAMMPS_SOURCE_DIR "${TADAH_LAMMPS_SOURCE_DIR}")

  cmake_language(DEFER CALL _tadah_apply_after_lammps)
endif()

