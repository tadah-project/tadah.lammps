# Tadah ML-TADAH integration for LAMMPS (no LAMMPS file edits)

This integration enables the ML-TADAH pair style in LAMMPS using CMake injection
via CMAKE_PROJECT_INCLUDE.

It requires no edits to the LAMMPS source tree and does not require copying
ML-TADAH into lammps/src. Instead, ML-TADAH sources are compiled directly from
the Tadah.LAMMPS repository into the LAMMPS build.

## Requirements

- CMake >= 3.20
- A working C++ toolchain for LAMMPS
- A local checkout of Tadah.LAMMPS (recommended workflow)
- Tadah.LAMMPS must NOT force CMAKE_CXX_STANDARD=11 globally when used as a dependency
  (LAMMPS requires C++17)

## Recommended workflow (local Tadah.LAMMPS checkout)

1. Have a local Tadah.LAMMPS checkout:

   /path/to/Tadah.LAMMPS

   The integration scripts are inside:

   /path/to/Tadah.LAMMPS/integration/

2. Configure LAMMPS (from the LAMMPS root directory):

```sh
cmake -S cmake -B build \
  -DPKG_ML-TADAH=ON \
  -DCMAKE_PROJECT_INCLUDE=/path/to/Tadah.LAMMPS/integration/TadahLAMMPSProjectInclude.cmake \
  -DTADAH_LAMMPS_SOURCE_DIR=/path/to/Tadah.LAMMPS
```

3. Build:

```sh
cmake --build build
```

## How it works (latest LAMMPS behavior)

- LAMMPS generates build/styles/style_pair.h as an include-only file that lists:
  #include "pair_*.h"
- This integration appends one or more lines like:
  #include "pair_tadah.h"
  to that generated file.
- Your pair_tadah.h contains:
  #ifdef PAIR_CLASS
  PairStyle(tadah,PairTadah)
  #else
  ...
  #endif
  so the style is registered automatically when the header is included.
- ML-TADAH/*.cpp from Tadah.LAMMPS is compiled directly into the lammps target.

## Script options (all available options)

All options below are passed to the LAMMPS CMake configure command (cmake -S cmake -B build ...).

### Option 1: PKG_ML-TADAH (ON/OFF)

Enable or disable the integration.

Example:

```sh
-DPKG_ML-TADAH=ON
```

Disable:

```sh
-DPKG_ML-TADAH=OFF
```

### Option 2: CMAKE_PROJECT_INCLUDE (path)

CMake hook that injects the integration into the LAMMPS configure run.

Required when using this integration.

Example:

```sh
-DCMAKE_PROJECT_INCLUDE=/path/to/Tadah.LAMMPS/integration/TadahLAMMPSProjectInclude.cmake
```

### Option 3: TADAH_LAMMPS_SOURCE_DIR (path)

Path to a local checkout of Tadah.LAMMPS.

This is required for the recommended workflow (latest-only integration). If this is
empty, the current scripts will error out (by design, to avoid unexpected behavior).

Example:

```sh
-DTADAH_LAMMPS_SOURCE_DIR=/path/to/Tadah.LAMMPS
```

### Option 4: TADAH_LAMMPS_QUIET (ON/OFF)

Request reduced configure output from Tadah.LAMMPS when it is built as a dependency.

Note: Tadah.LAMMPS must define and honor this option in its top-level CMakeLists.txt
for it to have any effect.

Default in the integration scripts: ON

Examples:

Quiet (recommended):

```sh
-DTADAH_LAMMPS_QUIET=ON
```

Verbose:

```sh
-DTADAH_LAMMPS_QUIET=OFF
```

### Option 5: TADAH_LAMMPS_PAIR_STYLES (list)

Semicolon-separated list of pair styles to register.

Format:

- style:Class[:header]

Latest-only behavior:
- the script uses this list to determine which header(s) to include into style_pair.h
- the Class field is currently not used by the script, but is kept for future expansion
- if header is omitted it defaults to: pair_<style>.h

Default:

- tadah:PairTadah

Examples:

A) Default (single style):

```sh
-DTADAH_LAMMPS_PAIR_STYLES="tadah:PairTadah"
```

B) Explicit header name:

```sh
-DTADAH_LAMMPS_PAIR_STYLES="tadah:PairTadah:pair_tadah.h"
```

C) Multiple styles (future use):

```sh
-DTADAH_LAMMPS_PAIR_STYLES="tadah:PairTadah;other:PairOther:pair_other.h"
```

## Full example command lines

### Minimal recommended command

```sh
cmake -S cmake -B build \
  -DPKG_ML-TADAH=ON \
  -DCMAKE_PROJECT_INCLUDE=/path/to/Tadah.LAMMPS/integration/TadahLAMMPSProjectInclude.cmake \
  -DTADAH_LAMMPS_SOURCE_DIR=/path/to/Tadah.LAMMPS
```

### Verbose Tadah.LAMMPS output (debugging)

```sh
cmake -S cmake -B build \
  -DPKG_ML-TADAH=ON \
  -DCMAKE_PROJECT_INCLUDE=/path/to/Tadah.LAMMPS/integration/TadahLAMMPSProjectInclude.cmake \
  -DTADAH_LAMMPS_SOURCE_DIR=/path/to/Tadah.LAMMPS \
  -DTADAH_LAMMPS_QUIET=OFF
```

### Custom style list example

```sh
cmake -S cmake -B build \
  -DPKG_ML-TADAH=ON \
  -DCMAKE_PROJECT_INCLUDE=/path/to/Tadah.LAMMPS/integration/TadahLAMMPSProjectInclude.cmake \
  -DTADAH_LAMMPS_SOURCE_DIR=/path/to/Tadah.LAMMPS \
  -DTADAH_LAMMPS_PAIR_STYLES="tadah:PairTadah:pair_tadah.h"
```
