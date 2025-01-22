# Tadah!LAMMPS

This repository provides the LAMMPS library (`tadah.lammps`) and plugin (`ML-TADAH`) required to run interatomic potentials developed with **Tadah!MLIP**.

**Note:** LAMMPS compiled with **Tadah!LAMMPS** is also needed for building **Tadah!MLIP** with hyperparameter optimization toolkit (HPO module). Please refer to the online documentation for more information.

[https://tadah.readthedocs.io/](https://tadah.readthedocs.io/)

In the near future, we intend to distribute this plugin with the LAMMPS source code. Until then, it must be cloned from this repository. Full installation instructions are available in the online documentation. However, if you intend to only run MD simulations using **Tadah!** potentials, the essential steps are listed here. The installation process follows the standard steps for building LAMMPS with external plugins: [https://docs.lammps.org/Build_extras.html](https://docs.lammps.org/Build_extras.html).

## Installation Instructions

1. Clone this repository to `lammps/lib`.  
   The directory name must match the repository name: `tadah.lammps`.
   ```bash
   # cd to /path/to/lammps/lib
   git clone https://git.ecdf.ed.ac.uk/tadah/tadah.lammps.git

2. In `lammps/src`, issue the following commands:

   ```bash
   make lib-tadah.lammps
   make yes-ml-tadah
   make serial

The first command compiles the **Tadah!LAMMPS** library and copies the `ML-TADAH` plugin to the `lammps/src` directory. The second command adds `ML-TADAH` to the list of packages compiled by LAMMPS. Finally, `make serial` recompiles the LAMMPS software. Please note that this command compiles a standard, non-parallel version of LAMMPS. Users are expected to adjust this as needed.

