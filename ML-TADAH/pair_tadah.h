/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government rgrid_2bins
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*----------------------------------------------------------------------*/
/* Contributing author: Marcin Kirsz | marcin.kirsz@ed.ac.uk            */
/* The University of Edinburgh       | https://git.ecdf.ed.ac.uk/tadah  */
/*----------------------------------------------------------------------*/

#ifdef PAIR_CLASS
// clang-format off
PairStyle(tadah,PairTadah)
// clang-format on
#else

#ifndef LMP_PAIR_TADAH_H
#define LMP_PAIR_TADAH_H

#include "pair.h"
#include <memory>
#include <tadah/tadah.lammps/lammps_tadah.h>

namespace LAMMPS_NS {

  class PairTadah : public Pair {
    private:
      /** Linear kernel + two-body descriptor */
      void compute_2b(int eflag, int vflag);

      /** (Non)Linear +  2b+mb or mb  */
      void compute_2b_mb_half(int eflag, int vflag);
      void compute_2b_mb_full(int eflag, int vflag);

      /** dimer-dimer: */
      void compute_dimers(int eflag, int vflag);
      std::vector<std::pair<int,int>> midx;

      std::unique_ptr<tadah::LammpsTadah> lt;

      int nmax=0;
      void allocate();

      double single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                double /*factor_coul*/, double /*factor_lj*/, double &fforce);

    public:
      void say_hello() const {
        std::cout << "Hello from PairTadah!" << std::endl;
      }

      bool load(int argc, char **argv);                          // CLI
      bool load(const tadah::core::Context& ctx);                // copy
      bool load(tadah::core::Context ctx);                       // move

      PairTadah(class LAMMPS *);
      ~PairTadah() override;

      // REQUIRED
      void compute(int, int) override;
      void settings(int, char **) override;
      void coeff(int, char **) override;

      // OPTIONAL
      void init_style() override;
      double init_one(int i, int j) override;
      int pack_forward_comm(int, int *, double *, int, int *) override;
      void unpack_forward_comm(int, int, double *) override;
      int pack_reverse_comm(int, int, double *) override;
      void unpack_reverse_comm(int, int *, double *) override;
  };
}
#endif
#endif
