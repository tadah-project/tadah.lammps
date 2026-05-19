[![pipeline](https://git.ecdf.ed.ac.uk/tadah/tadah.lammps/badges/develop/pipeline.svg)](https://git.ecdf.ed.ac.uk/tadah/tadah.lammps/-/commits/develop)
[![coverage](https://git.ecdf.ed.ac.uk/tadah/tadah.lammps/badges/develop/coverage.svg)](https://git.ecdf.ed.ac.uk/tadah/tadah.lammps/-/commits/develop)
[![license](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![version](https://img.shields.io/badge/version-1.3.0--beta.1-blue.svg)](https://git.ecdf.ed.ac.uk/tadah/tadah.lammps/-/tags)

# ML‑TADAH for LAMMPS (Tadah.LAMMPS)

This repository adds **ML‑TADAH** to **LAMMPS** as a **pair style** (`pair_style tadah`) for running interatomic potentials developed with **Tadah!MLIP**.

> **Note:** ML‑TADAH / Tadah.LAMMPS is **not part of the official LAMMPS distribution**. You must clone this repository yourself and enable it in your LAMMPS build.

Documentation (recommended starting point, incl. extended instructions):
https://tadah.readthedocs.io/

---

## Installation overview

There are two supported ways to build LAMMPS with ML‑TADAH:

1. **CMake build with LAMMPS (recommended)**
2. **Traditional make-based build with LAMMPS (deprecated)**

In both cases, you first **clone Tadah.LAMMPS** (this repository), since it is not shipped with LAMMPS.

For general LAMMPS build options and prerequisites, see:
https://docs.lammps.org/Build_cmake.html

---

## 1) Build with CMake (recommended)

### Step 1 — Clone Tadah.LAMMPS

Clone this repository anywhere you like:

```bash
git clone https://git.ecdf.ed.ac.uk/tadah/tadah.lammps.git /path/to/Tadah.LAMMPS
```

Optional: choose the Tadah.LAMMPS branch to use when configuring (default is `main`):

- `-DTADAH_GIT_BRANCH=main` (default)
- `-DTADAH_GIT_BRANCH=develop`

You can also `git checkout develop` in `/path/to/Tadah.LAMMPS` if you prefer managing the branch directly.

---

### Step 2 — Configure LAMMPS and enable ML‑TADAH

From the **LAMMPS root directory**, run your usual LAMMPS CMake configure command **and add** the ML‑TADAH integration options below **to that command**:

```bash
  -DPKG_ML-TADAH=ON \
  -DCMAKE_PROJECT_INCLUDE=/path/to/Tadah.LAMMPS/integration/TadahLAMMPSProjectInclude.cmake \
  -DTADAH_LAMMPS_SOURCE_DIR=/path/to/Tadah.LAMMPS
```

Those three options are the "integration hook": they enable the `ML-TADAH` package and tell LAMMPS where to find the Tadah.LAMMPS sources.

#### Example A: minimal LAMMPS configure + ML‑TADAH
Run from the LAMMPS root directory:
```bash
cmake -S cmake -B build \
  -DPKG_ML-TADAH=ON \
  -DCMAKE_PROJECT_INCLUDE=/path/to/Tadah.LAMMPS/integration/TadahLAMMPSProjectInclude.cmake \
  -DTADAH_LAMMPS_SOURCE_DIR=/path/to/Tadah.LAMMPS
```

#### Example B: use the LAMMPS GCC preset + build a serial binary + ML‑TADAH

This example uses the LAMMPS preset file `lammps/cmake/presets/gcc.cmake`, enables ML‑TADAH, selects the Tadah branch `develop`, and sets the LAMMPS machine target to `serial`:

```bash
cmake -C cmake/presets/gcc.cmake -S cmake -B build \
  -DPKG_ML-TADAH=ON \
  -DCMAKE_PROJECT_INCLUDE=/home/s1351949/PROJECTS/TADAH/Tadah.LAMMPS/integration/TadahLAMMPSProjectInclude.cmake \
  -DTADAH_LAMMPS_SOURCE_DIR=/home/s1351949/PROJECTS/TADAH/Tadah.LAMMPS \
  -DTADAH_GIT_BRANCH=develop \
  -DLAMMPS_MACHINE=serial
```

---

### Step 3 — Build

```bash
cmake --build build
```

Parallel builds are supported; on Linux a common choice is:

```bash
cmake --build build --parallel $(nproc)
```

Result: the resulting LAMMPS binary includes **ML‑TADAH** and provides `pair_style tadah`.

---

## 2) Traditional make-based build (deprecated)

This follows the classic LAMMPS `make` workflow. LAMMPS documentation for building “extras”:
https://docs.lammps.org/Build_extras.html

### Step 1 — Clone Tadah.LAMMPS into `lammps/lib`

From your LAMMPS source tree:

```bash
cd /path/to/lammps/lib
git clone https://git.ecdf.ed.ac.uk/tadah/tadah.lammps.git tadah.lammps
```

(Optionally select a branch, e.g. `develop`, if desired.)

### Step 2 — Build the library and enable the ML‑TADAH package

```bash
cd /path/to/lammps/src
make lib-tadah.lammps
make yes-ml-tadah
make serial
```

Choose the appropriate LAMMPS make target for your system (MPI/OpenMP, etc.) instead of `serial` if needed.

---

## What’s in this repository?

- **LAMMPS pair style sources**: **ML‑TADAH** (`pair_style tadah`)
- Supporting **C++ library** used by the integration: `tadah.lammps`
- CMake integration hook for LAMMPS:
  `integration/TadahLAMMPSProjectInclude.cmake`

For usage, examples, and troubleshooting, see:
https://tadah.readthedocs.io/
