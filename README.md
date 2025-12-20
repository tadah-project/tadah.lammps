# Tadah!LAMMPS

This repository provides:

- a C++ library for LAMMPS: `tadah.lammps`
- a LAMMPS pair style (plugin sources): `ML-TADAH` (`pair_style tadah`)

It is used to run interatomic potentials developed with Tadah!MLIP.

Documentation:
[https://tadah.readthedocs.io/](https://tadah.readthedocs.io/)

## Build with CMake (preferred)

This is the recommended way. It does not require editing any LAMMPS files and does not require copying `ML-TADAH` into `lammps/src`.

1. Clone Tadah!LAMMPS anywhere (example):
   ```bash
   git clone https://git.ecdf.ed.ac.uk/tadah/tadah.lammps.git /path/to/Tadah.LAMMPS
   ```

2. Configure LAMMPS with the integration hook (from the LAMMPS root directory):
   ```bash
   cmake -S cmake -B build \
     -DPKG_ML-TADAH=ON \
     -DCMAKE_PROJECT_INCLUDE=/path/to/Tadah.LAMMPS/integration/TadahLAMMPSProjectInclude.cmake \
     -DTADAH_LAMMPS_SOURCE_DIR=/path/to/Tadah.LAMMPS
   ```

3. Build:
   ```bash
   cmake --build build --parallel 10
   ```

Result: the LAMMPS binary will include `pair_style tadah`.

## Traditional make-based build (supported, but CMake is preferred)

If you must use the classic make-based workflow, follow the LAMMPS "Build extras" docs:
https://docs.lammps.org/Build_extras.html

Typical steps:

1. Clone this repository into `lammps/lib` (directory name should be `tadah.lammps`):
   ```bash
   cd /path/to/lammps/lib
   git clone https://git.ecdf.ed.ac.uk/tadah/tadah.lammps.git
   ```

2. Build and enable the package from `lammps/src`:
   ```bash
   cd /path/to/lammps/src
   make lib-tadah.lammps
   make yes-ml-tadah
   make serial
   ```

Notes:
- `make serial` builds a non-MPI binary; use the appropriate make target for MPI/OpenMP on your system.
- The CMake workflow above is preferred because it avoids manual file copying and is easier to reproduce.
