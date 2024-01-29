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

//#include "tadah/libs/Eigen/src/Core/util/Constants.h"
#ifdef PAIR_CLASS

PairStyle(tadah/tadah,PairTadah)

#else

#ifndef LMP_PAIR_TADAH_TADAH_H
#define LMP_PAIR_TADAH_TADAH_H

#include "pair.h"
#include "../../lib/tadah.md/lammps_tadah.h"

namespace LAMMPS_NS {


class PairTadah : public Pair {
    private:
        /** Linear kernel + two-body descriptor */
        void compute_2b(int eflag, int vflag);

        /** (Non)Linear +  2b+mb or mb  */
        void compute_2b_mb_half(int eflag, int vflag);
        void compute_2b_mb_full(int eflag, int vflag);

        /** dimer-dimer: quadrupole interactions */
        void compute_dimers(int eflag, int vflag);
        std::vector<std::pair<int,int>> midx;

        TADAH::LammpsTadah lt;
        int nmax=0;
    public:
        PairTadah(class LAMMPS *);
        virtual ~PairTadah();
        virtual void compute(int, int);
        void settings(int, char **);
        void coeff(int, char **);
        double init_one(int, int);
        void init_style();
        virtual int pack_forward_comm(int, int *, double *, int, int *);
        virtual void unpack_forward_comm(int, int, double *);
        int pack_reverse_comm(int, int, double *);
        void unpack_reverse_comm(int, int *, double *);

        virtual void allocate();
        void *extract(const char *, int &) override;

    protected:
        double **scale;
};
}
#endif
#endif
