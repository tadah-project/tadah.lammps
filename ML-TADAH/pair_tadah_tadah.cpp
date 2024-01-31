/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government rgrid_2bins
certain rights in this software.  This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Marcin Kirsz (The University of Edinburgh)
------------------------------------------------------------------------- */

#include "pair_tadah_tadah.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "molecule.h"
#include "domain.h"

using namespace LAMMPS_NS;
/* ---------------------------------------------------------------------- */

PairTadah::PairTadah(LAMMPS *lmp) : Pair(lmp)
{
    writedata = 0;
    restartinfo = 0;
    //single_enable = 0;
    manybody_flag = 0;  // TODO req if else
    scale = nullptr;
}
/* ---------------------------------------------------------------------- */

PairTadah::~PairTadah()
{
    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(scale);
    }
}
/* ---------------------------------------------------------------------- */
void PairTadah::compute_2b(int eflag, int vflag)
{
    // Works only for linear regression with two-body descriptor.
    // Assumes weights are "normalised" so we
    // do not need to normalise descriptors.
    //
    // Works for both newton on or off.
    int i,j,ii,jj,jnum;
    double xtmp,ytmp,ztmp,fpair,ei;
    double rij_sq;
    int *ilist,*jlist,*numneigh,**firstneigh;

    ev_init(eflag,vflag);

    double **x = atom->x;
    double **f = atom->f;

    int *type = atom->type;
    int nlocal = atom->nlocal;
    Vec3d delij;
    //double delij[3];

    aed_type2 aed(lt.c.get<size_t>("DSIZE"));
    fd_type fd(lt.c.get<size_t>("DSIZE"));

    int inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    int newton_pair = force->newton_pair;

    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];
        //double iweight = lt.weights[type[i]-1];

        for (jj = 0; jj < jnum; jj++) {
            aed.setZero();
            j = jlist[jj];
            j &= NEIGHMASK;
            delij[0] = xtmp - x[j][0];
            delij[1] = ytmp - x[j][1];
            delij[2] = ztmp - x[j][2];
            rij_sq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

            if (rij_sq > lt.cutoff_max_sq) continue;

            double jweight = lt.weights[type[j]-1];
            //double ijw = iweight+jweight;   // TODO is this correct?
            double rij = sqrt(rij_sq);
            //double fc_ij = ijw*lt.S.c2b->calc(rij);
            //double fcp_ij = ijw*lt.S.c2b->calc_prime(rij);
            double fc_ij = jweight*lt.S.c2b->calc(rij);
            double fcp_ij = jweight*lt.S.c2b->calc_prime(rij);
            lt.S.d2b->calc_all(rij,rij_sq,fc_ij,fcp_ij,aed, fd);

            fpair = scale[type[i]][type[j]] * lt.model->fpredict(fd,aed,0)/rij;
            ei = scale[type[i]][type[j]] * lt.model->epredict(aed);

            f[i][0] += fpair*delij[0];
            f[i][1] += fpair*delij[1];
            f[i][2] += fpair*delij[2];

            if (newton_pair || j < nlocal) {
                f[j][0] -= fpair*delij[0];
                f[j][1] -= fpair*delij[1];
                f[j][2] -= fpair*delij[2];
            }

            if (evflag) ev_tally(i,j,nlocal,newton_pair,ei,0.0,fpair,delij[0],delij[1],delij[2]);
        }

        aed.setZero();
        if (lt.bias)
            aed(0) = 1.0;
        ei = lt.model->epredict(aed);   // This is to get constant shift in E
        if (eflag_atom) eatom[i] += ei;
        if (eflag_global) eng_vdwl += ei;
    }

    if (vflag_fdotr)  virial_fdotr_compute();
}
void PairTadah::compute_2b_mb_half(int eflag, int vflag)
{
    // Works for both D2 and DM descriptor with Linear models
    int i,j,ii,jj,jnum;
    double xtmp,ytmp,ztmp;
    double rij_sq;
    int *ilist,*jlist,*numneigh,**firstneigh;

    ev_init(eflag,vflag);

    double **x = atom->x;
    double **f = atom->f;

    int *type = atom->type;
    int nlocal = atom->nlocal;
    //int nall = nlocal + atom->nghost;
    Vec3d delij;
    //double delij[3];

    int inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    int newton_pair = force->newton_pair;

    // fpsize x 3-directions
    fd_type fd(lt.c.get<size_t>("DSIZE"));

    // grow arrays if necessary
    if (atom->nmax > nmax) {
        nmax = atom->nmax;
        size_t s = lt.S.dmb->rhoi_size()+ lt.S.dmb->rhoip_size();
		lt.rhos.resize(s,nmax);
		lt.aeds.resize(lt.c.get<size_t>("DSIZE"),nmax);
    }
    lt.rhos.setZero();

    lt.aeds.setZero();
    if (lt.bias)
        for (size_t n=0; n<lt.aeds.cols(); ++n )lt.aeds(0,n) = 1.0;

    //if (newton_pair || lt.initmb) {
    //    for (i = 0; i < nall; i++) {
    //        lt.aeds[i].setZero();
    //        if (lt.bias)
    //            lt.aeds[i](0) = 1.0;
    //    }
    //}
    //else for (i = 0; i < nlocal; i++) {
    //    lt.aeds[i].setZero();
    //    if (lt.bias)
    //        lt.aeds[i](0) = 1.0;
    //}

    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        double iweight = lt.weights[type[i]-1];
        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;

            delij[0] = xtmp - x[j][0];
            delij[1] = ytmp - x[j][1];
            delij[2] = ztmp - x[j][2];
            rij_sq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

            if (rij_sq > lt.cutoff_max_sq) continue;
            double rij = sqrt(rij_sq);
            double jweight = lt.weights[type[j]-1];

            if (lt.S.c2b->get_rcut() > rij) {
                double fc2b = lt.S.c2b->calc(rij);
                lt.S.d2b->calc_aed(rij,rij_sq,jweight*fc2b,lt.aeds.col(i));

                if (newton_pair || j < nlocal) {
                    lt.S.d2b->calc_aed(rij,rij_sq,iweight*fc2b,lt.aeds.col(j));
                }
            }

            if (lt.S.cmb->get_rcut() > rij) {
                double fcmb = lt.S.cmb->calc(rij);
                lt.S.dmb->calc_rho(rij,rij_sq,jweight*fcmb,delij,lt.rhos.col(i));

                if (newton_pair || j < nlocal) {
                    lt.S.dmb->calc_rho(rij,rij_sq,iweight*fcmb,-delij,lt.rhos.col(j));
                }
            }
        }
    }

    // Communicate and sum densities for EAM.
    // Also communicate 2b aeds if newton is on.
    if (comm_reverse && newton_pair) {
        comm->reverse_comm(this);
        //comm->reverse_comm_pair(this);
    }

    // compute AEFP for all eam-types
    // for EAM the derivative of embedding energy at each atom
    // is also computed here
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        // we have to scale 2b energy by 0.5
        // we can't do it in the loop above because
        // if dummy cutoff is used or d2_eam with dm_eam
        // cutoff is ignored and that is the only way now to
        // introduce factor into 2b aed.
        // Hence we scale 2b part of descriptor by 0.5
        // before MB aed is calculated
        lt.aeds.col(i)*=0.5;

        lt.S.dmb->calc_aed(lt.rhos.col(i), lt.aeds.col(i));

        if (lt.norm)
            lt.model->norm.normalise(lt.aeds.col(i));

        // TODO The proper way is to scale manybody energy by scale[i][i]
        // and pairwise energy by scale[i][j]
        // In otherwords thermodynamic integration will
        // probably work for single species systems only
        // This applies to all tadah compute methods apart from compute_2b
        // which should be fine
        double ei = scale[type[i]][type[i]] * lt.model->epredict(lt.aeds.col(i));
        if (eflag_atom) eatom[i] = ei;
        if (eflag_global) eng_vdwl += ei;
    }

    // communicate the derivative of the embedding function
    if (comm_forward)
        comm->forward_comm(this);

    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];

        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];

        jlist = firstneigh[i];
        jnum = numneigh[i];
        double iweight = lt.weights[type[i]-1];
        for (jj = 0; jj < jnum; jj++) {
            fd.setZero();
            j = jlist[jj];
            j &= NEIGHMASK;
            delij[0] = xtmp - x[j][0];
            delij[1] = ytmp - x[j][1];
            delij[2] = ztmp - x[j][2];
            rij_sq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];
            if (rij_sq > lt.cutoff_max_sq) continue;
            double rij = sqrt(rij_sq);
            double jweight = lt.weights[type[j]-1];

            if (lt.S.c2b->get_rcut() > rij) {
                //double ijw = iweight+jweight;
                double fc2b = lt.S.c2b->calc(rij);
                double fcp2b = lt.S.c2b->calc_prime(rij);
                //lt.S.d2b->calc_dXijdri(rij,rij_sq,ijw*fc2b,ijw*fcp2b,fd);
                lt.S.d2b->calc_dXijdri(rij,rij_sq,jweight*fc2b,jweight*fcp2b,fd);
            }

            int mode=0;
            if (lt.S.cmb->get_rcut() > rij) {
                double fcmb = lt.S.cmb->calc(rij);
                double fcpmb = lt.S.cmb->calc_prime(rij);
                mode = lt.S.dmb->calc_dXijdri_dXjidri(rij,rij_sq,delij,fcmb,fcpmb,lt.rhos.col(i),lt.rhos.col(j), fd, iweight,jweight);
            }


            if (mode==0) {
                if (lt.norm)
                    lt.model->norm.normalise(fd,0);
                // mode==0 only x-component is calculated for both 2b and mb
                double fpair = scale[type[i]][type[j]] * lt.model->fpredict(fd,lt.aeds.col(i),0)/rij;

                f[i][0] += delij[0]*fpair;
                f[i][1] += delij[1]*fpair;
                f[i][2] += delij[2]*fpair;

                if (newton_pair || j < nlocal) {
                    f[j][0] -= delij[0]*fpair;
                    f[j][1] -= delij[1]*fpair;
                    f[j][2] -= delij[2]*fpair;
                }
                if (evflag) ev_tally(i,j,nlocal,newton_pair,0.0,0.0,fpair,delij[0],delij[1],delij[2]);
            }
            else {
                // mode!=0, S.d2b calculates x-comp only, S.dmb calcs all 3-directions

                // copy 2b descriptors first and
                // scale all 2b components by delij/rij
                // this is inefficient but will do for now... TODO
                for(size_t s=lt.bias; s<lt.S.d2b->size()+lt.bias; ++s) {
                    fd(s,0) /= rij;
                    fd(s,1) = fd(s,0);
                    fd(s,2) = fd(s,0);
                    fd(s,0) *= delij[0];
                    fd(s,1) *= delij[1];
                    fd(s,2) *= delij[2];
                }
                if (lt.norm)
                    lt.model->norm.normalise(fd);

                // for mode!=0 we don't need to scale mb part by delij/rij
                // as S.dmb calculator does it.
                Vec3d f_eam = scale[type[i]][type[j]] * lt.model->fpredict(fd,lt.aeds.col(i));

                f[i][0] += f_eam(0);
                f[i][1] += f_eam(1);
                f[i][2] += f_eam(2);

                if (newton_pair || j < nlocal) {
                    f[j][0] -= f_eam(0);
                    f[j][1] -= f_eam(1);
                    f[j][2] -= f_eam(2);
                }

                if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,f_eam(0),f_eam(1),f_eam(2),delij[0],delij[1],delij[2]);
            }
        }
    }

    if (vflag_fdotr)  virial_fdotr_compute();

}
void PairTadah::compute_2b_mb_full(int eflag, int vflag)
{
    int i,j,ii,jj,jnum;
    double xtmp,ytmp,ztmp;
    double rij_sq;
    int *ilist,*jlist,*numneigh,**firstneigh;

    ev_init(eflag,vflag);

    double **x = atom->x;
    double **f = atom->f;

    int *type = atom->type;
    int nlocal = atom->nlocal;
    //int nall = nlocal + atom->nghost;
    Vec3d delij;

    int inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    int newton_pair = force->newton_pair;

    // fpsize x 3-directions
    fd_type fd(lt.c.get<size_t>("DSIZE"));

    // grow arrays if necessary
    if (atom->nmax > nmax) {
        nmax = atom->nmax;
        size_t s = lt.S.dmb->rhoi_size()+ lt.S.dmb->rhoip_size();
		lt.rhos.resize(s,nmax);
		lt.aeds.resize(lt.c.get<size_t>("DSIZE"), nmax);
    }

    //for (i = 0; i < nlocal; i++) {
    //    lt.aeds[i].setZero();
    //    if (lt.bias)
    //        lt.aeds[i](0) = 1.0;
    //}
    lt.rhos.setZero();
    lt.aeds.setZero();
    if (lt.bias)
        for (size_t n=0; n<lt.aeds.cols(); ++n )lt.aeds(0,n) = 1.0;

    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;
            delij[0] = xtmp - x[j][0];
            delij[1] = ytmp - x[j][1];
            delij[2] = ztmp - x[j][2];
            rij_sq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

            if (rij_sq > lt.cutoff_max_sq) continue;
            double rij = sqrt(rij_sq);
            double jweight = lt.weights[type[j]-1];

            if (lt.S.c2b->get_rcut() > rij) {
                double fc2b = lt.S.c2b->calc(rij);
                lt.S.d2b->calc_aed(rij,rij_sq,jweight*fc2b,lt.aeds.col(i));
            }

            if (lt.S.cmb->get_rcut() > rij) {
                double fcmb = lt.S.cmb->calc(rij);
                lt.S.dmb->calc_rho(rij,rij_sq,jweight*fcmb,delij,lt.rhos.col(i));
            }
        }
    }

    // compute AEFP for all eam-types
    // for EAM the derivative of embedding energy at each atom
    // is also computed here
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        // we have to scale 2b energy by 0.5
        // we can't do it in the loop above because
        // if dummy cutoff is used or d2_eam with dm_eam
        // cutoff is ignored and that is the only way now to
        // introduce factor into 2b aed.
        // Hence we scale 2b part of descriptor by 0.5
        // before MB aed is calculated
        lt.aeds.col(i)*=0.5;

        lt.S.dmb->calc_aed(lt.rhos.col(i), lt.aeds.col(i));

        if (lt.norm)
            lt.model->norm.normalise(lt.aeds.col(i));

        double ei = scale[type[i]][type[i]] * lt.model->epredict(lt.aeds.col(i));
        if (eflag_atom) eatom[i] = ei;
        if (eflag_global) eng_vdwl += ei;
    }

    // communicate the derivative of the embedding function
    if (comm_forward)
        comm->forward_comm(this);

    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];

        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];

        jlist = firstneigh[i];
        jnum = numneigh[i];
        // double iweight = lt.weights[type[i]-1];
        for (jj = 0; jj < jnum; jj++) {
            fd.setZero();
            j = jlist[jj];
            j &= NEIGHMASK;
            delij[0] = xtmp - x[j][0];
            delij[1] = ytmp - x[j][1];
            delij[2] = ztmp - x[j][2];
            rij_sq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];
            if (rij_sq > lt.cutoff_max_sq) continue;
            double rij = sqrt(rij_sq);

            double jweight = lt.weights[type[j]-1];

            if (lt.S.c2b->get_rcut() > rij) {
                double fc2b = 0.5*jweight*lt.S.c2b->calc(rij);
                double fcp2b = 0.5*jweight*lt.S.c2b->calc_prime(rij);
                lt.S.d2b->calc_dXijdri(rij,rij_sq,fc2b,fcp2b,fd);
            }

            int mode=0;
            if (lt.S.cmb->get_rcut() > rij) {
                double fcmb = jweight*lt.S.cmb->calc(rij);
                double fcpmb = jweight*lt.S.cmb->calc_prime(rij);
                mode = lt.S.dmb->calc_dXijdri(rij,rij_sq,delij,fcmb,fcpmb,lt.rhos.col(i), fd);
                //double fcmb = lt.S.cmb->calc(rij);
                //double fcpmb = lt.S.cmb->calc_prime(rij);
                //mode = lt.S.dmb->calc_dXijdri_dXjidri(rij,rij_sq,delij,fcmb,fcpmb,rho[i],rho[j], fd, iweight,jweight);
            }


            if (mode==0) {
                if (lt.norm)
                    lt.model->norm.normalise(fd,0);
                // mode==0 only x-component is calculated for both 2b and mb
                double fpair = scale[type[i]][type[j]] * lt.model->fpredict(fd,lt.aeds.col(i),0)/rij;

                f[i][0] += delij[0]*fpair;
                f[i][1] += delij[1]*fpair;
                f[i][2] += delij[2]*fpair;

                f[j][0] -= delij[0]*fpair;
                f[j][1] -= delij[1]*fpair;
                f[j][2] -= delij[2]*fpair;

                if (evflag) ev_tally(i,j,nlocal,newton_pair,0.0,0.0,fpair,delij[0],delij[1],delij[2]);
            }
            else {
                // mode!=0, S.d2b calculates x-comp only, S.dmb calcs all 3-directions

                // copy 2b descriptors first and
                // scale all 2b components by delij/rij
                // this is inefficient but will do for now... TODO
                for(size_t s=lt.bias; s<lt.S.d2b->size()+lt.bias; ++s) {
                    fd(s,0) /= rij;
                    fd(s,1) = fd(s,0);
                    fd(s,2) = fd(s,0);
                    fd(s,0) *= delij[0];
                    fd(s,1) *= delij[1];
                    fd(s,2) *= delij[2];
                }
                if (lt.norm)
                    lt.model->norm.normalise(fd);

                // for mode!=0 we don't need to scale mb part by delij/rij
                // as S.dmb calculator does it.
                Vec3d f_eam = scale[type[i]][type[j]] * lt.model->fpredict(fd,lt.aeds.col(i));

                f[i][0] += f_eam(0);
                f[i][1] += f_eam(1);
                f[i][2] += f_eam(2);

                f[j][0] -= f_eam(0);
                f[j][1] -= f_eam(1);
                f[j][2] -= f_eam(2);

                if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,f_eam(0),f_eam(1),f_eam(2),delij[0],delij[1],delij[2]);
            }
        }
    }

    if (vflag_fdotr)  virial_fdotr_compute();

}
void PairTadah::compute_dimers(int eflag, int vflag)
{
    int ii,jj,jnum;
    int *ilist,*jlist,*numneigh,**firstneigh;

    ev_init(eflag,vflag);

    double **x = atom->x;
    double **f = atom->f;

    int nlocal = atom->nlocal;
    int *type = atom->type;

    int inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    int newton_pair = force->newton_pair;
    tagint *tag = atom->tag;

    tagint **special = atom->special;
    int **nspecial = atom->nspecial;

    fd_type fdIJ(lt.c.get<size_t>("DSIZE"));
    fd_type fdJI(lt.c.get<size_t>("DSIZE"));

    //n= 0, 1, 2, 3
    int i1,i2,j1,j2;
    Mat6R3C delM;
    double r_sq[6];
    double r[6];

    Vec3d delij_com;   // i-j CoM
    Vec3d i_com;   // i CoM
    Vec3d j_com;   // j CoM

    // if bond is not included begin summation from the third distance
    int Nstart= lt.dimer_bond_bool ? 0:2;
    double r_b = lt.dimer_r;
    double r_b_sq=r_b*r_b;

    double r_cut_com_sq =pow(lt.cutoff_max-r_b,2);

    // This is the tolerance for the bond distance
    double eps=1e-1;

    std::vector<int> markflag(atom->nlocal+atom->nghost,0);
    std::vector<int> jmarkflag(atom->nlocal+atom->nghost,0);
    //markflag.setZero();

    for (ii = 0; ii < inum; ii++) {
        i1 = ilist[ii];
        if (markflag[i1]) continue;

        i2 = atom->bond_atom[i1][0];
        i2 = atom->map(i2);
        i2 = domain->closest_image(i1,i2);

        //if (tag[i1]%2)
        //    i2 = atom->map(tag[i1] + 1);
        //else
        //    i2 = atom->map(tag[i1] - 1);
        //i2 = domain->closest_image(i1,i2);

        //int II = special[i1][0];
        //II = atom->map(II);
        //II = domain->closest_image(i1,II);
        //if (II!=i2)
        //    throw std::runtime_error("i1-II differ\n");

        //int III = atom->bond_atom[i1][0];
        //III = atom->map(III);
        //III = domain->closest_image(i1,III);
        //if (III!=i2)
        //    throw std::runtime_error("i1-III differ\n");

        if (markflag[i2]) continue;
        markflag[i2] = 1;
        markflag[i1] = 1;

        delM(0,0) = x[i1][0] - x[i2][0];
        delM(0,1) = x[i1][1] - x[i2][1];
        delM(0,2) = x[i1][2] - x[i2][2];
        r_sq[0] = delM.row(0)*delM.row(0);
        r[0] = sqrt(r_sq[0]);

        if (abs(r[0]-r_b)>eps) {
            throw std::runtime_error("i1-i2 bond is: " + std::to_string(r[0]));
        }

        jlist = firstneigh[i1];
        jnum = numneigh[i1];
        i_com[0] = 0.5*(x[i1][0]+x[i2][0]);
        i_com[1] = 0.5*(x[i1][1]+x[i2][1]);
        i_com[2] = 0.5*(x[i1][2]+x[i2][2]);

        jmarkflag = markflag;   // TODO just search markflag with j1 and j2
        for (jj = 0; jj < jnum; jj++) {
            j1 = jlist[jj];
            j1 &= NEIGHMASK;
            if (jmarkflag[j1]) continue;

            delM(2,0) = x[i1][0] - x[j1][0];
            delM(2,1) = x[i1][1] - x[j1][1];
            delM(2,2) = x[i1][2] - x[j1][2];
            r_sq[2] = delM.row(2)*delM.row(2);

            if (r_sq[2] > lt.cutoff_max_sq) continue;
            r[2] = sqrt(r_sq[2]);

            bool found=false;

            // but what if j1 is a ghost without a local copy?
            // we need a different approach in this case
            j2=-1;
            if (j1>=atom->nlocal) {
                // When full NN list is used than we should be able to search around i1
                // hence there is no need for ghost atoms lists
                for (int jj2 = 0; jj2 < jnum; jj2++) {
                    j2 = jlist[jj2];
                    j2 &= NEIGHMASK;
                    delM(1,0) = x[j1][0] - x[j2][0];
                    delM(1,1) = x[j1][1] - x[j2][1];
                    delM(1,2) = x[j1][2] - x[j2][2];
                    r_sq[1] = delM.row(1)*delM.row(1);
                    // the correct atom j2 is identified  when the distance
                    // between j1-j2 is approx equal to the bond length
                    if (abs(r_sq[1]-r_b_sq)<eps) {
                        found=true;
                        r[1] = sqrt(r_sq[1]);
                        break;
                    }
                }
            }
            else {
                int n = nspecial[j1][0]; // number of bonded atoms either 0 or 1

                if (n==0)
                    throw std::runtime_error("Ta-dah! No bonded atom found.\n");

                j2 = special[j1][0];    // j2 is a global index now, j1 is local
                j2 = atom->map(j2); // convert global j2 to local j2
                j2 &= NEIGHMASK;
                j2 = domain->closest_image(j1,j2);
                delM(1,0) = x[j1][0] - x[j2][0];
                delM(1,1) = x[j1][1] - x[j2][1];
                delM(1,2) = x[j1][2] - x[j2][2];
                r_sq[1] = delM.row(1)*delM.row(1);
                r[1] = sqrt(r_sq[1]);

            }

            if (j2==-1)
                throw std::runtime_error("Ta-dah! No bonded atom found v1.\n");

            // Finally mark but only once we have j1 and j2
            if (jmarkflag[j2]) continue;
            jmarkflag[j2] = 1;
            if (jmarkflag[j1]) continue;
            jmarkflag[j1] = 1;

            if (abs(r_sq[1]-r_b_sq)<eps) {
                found=true;
            }

            if (!found) {
                throw std::runtime_error("Ta-dah! No bonded atom found v2.\n");
            }


            // We should have all four atoms by now.
            // Each molecule-molecule interaction is computed just once
            // So here we apply forces the same way as we would do 
            // for the full NN list.

            // Check are ij CoM distance within the max CoM cutoff distance
            j_com[0] = 0.5*(x[j1][0]+x[j2][0]);
            j_com[1] = 0.5*(x[j1][1]+x[j2][1]);
            j_com[2] = 0.5*(x[j1][2]+x[j2][2]);
            delij_com = i_com - j_com;
            double r_com_sq = delij_com * delij_com;
            if (r_com_sq > r_cut_com_sq) {
                continue;
            }

            std::vector<int> IJ  = {i1,i2,j1,j2};

            if (lt.initmb) {
                lt.rhos.setZero();
            }
            lt.aeds.setZero();
            if (lt.bias)
                for (size_t  n=0; n<lt.aeds.cols(); ++n )lt.aeds(0,n) = 1.0;

            // Compute remaining 3 distances
            for (int n=3; n<6; ++n) {
                int I=IJ[midx[n].first];
                int J=IJ[midx[n].second];
                delM(n,0) = x[I][0] - x[J][0];
                delM(n,1) = x[I][1] - x[J][1];
                delM(n,2) = x[I][2] - x[J][2];
                r_sq[n] = delM.row(n)*delM.row(n);
                r[n] = sqrt(r_sq[n]);
            }

            // Compute densities for every atom.
            for (int n=Nstart; n<6; ++n) {

                if (lt.initmb && lt.S.cmb->get_rcut() >= r[n]) {
                    double fcmb = lt.S.cmb->calc(r[n]);
                    lt.S.dmb->calc_rho(r[n],r_sq[n],fcmb,delM.row(n),lt.rhos.col(midx[n].first));
                    lt.S.dmb->calc_rho(r[n],r_sq[n],fcmb,-delM.row(n),lt.rhos.col(midx[n].second));
                }
                if (lt.init2b && lt.S.c2b->get_rcut() >= r[n]) {
                    double fc2b = lt.S.c2b->calc(r[n]);
                    lt.S.d2b->calc_aed(r[n],r_sq[n],fc2b,lt.aeds.col(midx[n].first));
                    lt.S.d2b->calc_aed(r[n],r_sq[n],fc2b,lt.aeds.col(midx[n].second));
                }
            }

            // Compute two-body and many-body energy
            for (size_t n=0; n<4; ++n) {
                if (lt.initmb) {
                    lt.S.dmb->calc_aed(lt.rhos.col(n),lt.aeds.col(n));
                }
                if (lt.norm)
                    lt.model->norm.normalise(lt.aeds.col(n));

                if (IJ[n]<atom->nlocal) {
                    int I=IJ[midx[n].first];
                    int J=IJ[midx[n].second];
                    double ei = scale[type[I]][type[J]] *lt.model->epredict(lt.aeds.col(n));
                    if (eflag_atom) eatom[IJ[n]] = ei;
                    if (eflag_global) eng_vdwl += ei;
                }
            }

            // Compute forces
            for (int n=Nstart; n<6; ++n) {
                fdIJ.setZero();
                fdJI.setZero();
                int I=IJ[midx[n].first];
                int J=IJ[midx[n].second];
                if (lt.init2b && lt.S.c2b->get_rcut() > r[n]) {
                    double fc2b = lt.S.c2b->calc(r[n]);
                    double fcp2b = lt.S.c2b->calc_prime(r[n]);
                    lt.S.d2b->calc_dXijdri(r[n],r_sq[n],fc2b,fcp2b,fdIJ);
                }

                int mode=0;
                if (lt.initmb && lt.S.cmb->get_rcut() > r[n]) {
                    double fcmb = lt.S.cmb->calc(r[n]);
                    double fcpmb = lt.S.cmb->calc_prime(r[n]);
                    mode = lt.S.dmb->calc_dXijdri(r[n],r_sq[n],delM.row(n),
                            fcmb,fcpmb,lt.rhos.col(midx[n].first), fdIJ);
                    mode = lt.S.dmb->calc_dXijdri(r[n],r_sq[n],-delM.row(n),
                            fcmb,fcpmb,lt.rhos.col(midx[n].second), fdJI);
                }

                if (mode==0) {
                    // either compute 2b only or 2b+mb
                    if (lt.norm)
                        lt.model->norm.normalise(fdIJ,0);
                    double fpair = lt.model->fpredict(
                            fdIJ,lt.aeds.col(midx[n].first),0)/r[n];
                    fpair -= scale[type[I]][type[J]] *lt.model->fpredict(
                            -fdIJ,lt.aeds.col(midx[n].second),0)/r[n];

                    f[I][0] += fpair*delM(n,0);
                    f[I][1] += fpair*delM(n,1);
                    f[I][2] += fpair*delM(n,2);

                    f[J][0] -= fpair*delM(n,0);
                    f[J][1] -= fpair*delM(n,1);
                    f[J][2] -= fpair*delM(n,2);

                    if (evflag) ev_tally(I,J,nlocal,newton_pair,0,0,fpair,
                            delM(n,0),delM(n,1),delM(n,2));
                }
                else {
                    for(size_t s=lt.bias; s<lt.S.d2b->size()+lt.bias; ++s) {
                        fdIJ(s,0) /= r[n];
                        fdIJ(s,1) = fdIJ(s,0);
                        fdIJ(s,2) = fdIJ(s,0);
                        // fdIJ=fdJI for pairwise part
                        fdJI(s,0) = fdIJ(s,0);
                        fdJI(s,1) = fdIJ(s,0);
                        fdJI(s,2) = fdIJ(s,0);


                        fdIJ(s,0) *= delM(n,0);
                        fdIJ(s,1) *= delM(n,1);
                        fdIJ(s,2) *= delM(n,2);

                        fdJI(s,0) *= -delM(n,0);
                        fdJI(s,1) *= -delM(n,1);
                        fdJI(s,2) *= -delM(n,2);
                    }

                    if (lt.norm) {
                        lt.model->norm.normalise(fdIJ);
                        lt.model->norm.normalise(fdJI);
                    }

                    Vec3d f_eam = scale[type[I]][type[J]] *lt.model->fpredict(fdIJ,lt.aeds.col(midx[n].first));
                    f_eam -= lt.model->fpredict(fdJI,lt.aeds.col(midx[n].second));

                    f[I][0] += f_eam(0);
                    f[I][1] += f_eam(1);
                    f[I][2] += f_eam(2);

                    f[J][0] -= f_eam(0);
                    f[J][1] -= f_eam(1);
                    f[J][2] -= f_eam(2);

                    if (evflag) ev_tally_xyz(I,J,nlocal,newton_pair,0.0,0.0,
                            f_eam(0),f_eam(1),f_eam(2),
                            delM(n,0),delM(n,1),delM(n,2));
                }
            }
        }
    }

    if (vflag_fdotr)  virial_fdotr_compute();

}
void PairTadah::compute(int eflag, int vflag)
{
    if (lt.c.get<bool>("DIMER",0)) {
        compute_dimers(eflag,vflag);
    }
    else {
        if (lt.c.get<bool>("INITMB")) {
            if (lt.linear)
                compute_2b_mb_half(eflag,vflag);
            else
                compute_2b_mb_full(eflag,vflag);
        }
        else {
            if (lt.linear)
                compute_2b(eflag,vflag);
            else
                compute_2b_mb_full(eflag,vflag);
        }
    }

    // Check for NaN pressure.
    if (lt.c.get<bool>("CHECKPRESS")) {
        for (size_t i=0;i<6;++i) {
            if (std::isnan(virial[i])) {
                throw std::runtime_error("Ta-dah! ERROR: NaN pressure\n");
            }
        }
    }
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairTadah::allocate()
{
    if (!allocated) {
        int n = atom->ntypes;
        memory->create(setflag,n+1,n+1,"pair:setflag");
        memory->create(cutsq,n+1,n+1,"pair:cutsq");
        memory->create(scale,n+1,n+1,"pair:scale");
    }
    allocated = 1;
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairTadah::settings(int narg, char **arg)
{
    if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   SEt coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairTadah::coeff(int narg, char **arg)
{
    if (narg != 3 + atom->ntypes)
        error->all(FLERR,"Incorrect args for pair coefficients");

    if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
        error->all(FLERR,"Incorrect args for pair coefficients");

    if (!allocated) allocate();

    lt=TADAH::LammpsTadah(narg,arg);

    if (lt.c.get<bool>("INITMB"))
        manybody_flag++;
    if (lt.c.get<bool>("DIMER",0)) {
        manybody_flag=0;
        /*  [0] i1-i2    (i1)-[2]-(j1)
         *  [1] j1-j2     | \    / |
         *  [2] i1-j1     | [4] /  |
         *  [3] i2-j2    [0]  \/  [1]
         *  [4] i1-j2     |   /\   |
         *  [5] i2-j1     | [5] \  |   
         *                | /    \ |   
         *               (i2)-[3]-(j2)
         *
         *  This compute only works with rigid bonds.
         *  Style shake is your friend
         *  fix fixSHAKE all shake 1e-5 20 0 b ${bondlength}
         *  For now we assume that the distance between 
         *  atoms of two interacting molecules is always
         *  greater than the bond distance.
         * 
         *  i1---i2-----------------j1---j2
         *     |-------r_com-----------|   
         *   |---------r_max-------------|
         *
         */
        midx.resize(6);
        midx[0] = std::make_pair(0,1);
        midx[1] = std::make_pair(2,3);
        midx[2] = std::make_pair(0,2);
        midx[3] = std::make_pair(1,3);
        midx[4] = std::make_pair(0,3);
        midx[5] = std::make_pair(1,2);   
        lt.aeds.resize(lt.c.get<size_t>("DSIZE"),4);

        if (lt.initmb) {
            size_t s = lt.S.dmb->rhoi_size()+ lt.S.dmb->rhoip_size();
            lt.rhos.resize(s,4);
        }
    }
    if (lt.linear) {
        comm_reverse = lt.S.dmb->rhoi_size()+lt.S.d2b->size();
        comm_forward = lt.S.dmb->rhoip_size();
    }
    else {
        comm_reverse = lt.S.dmb->rhoi_size()+lt.S.d2b->size();
        comm_forward = lt.S.dmb->rhoip_size()+lt.S.d2b->size()+lt.S.dmb->size();
    }

    int ilo,ihi,jlo,jhi;
    utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
    utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);
    int count = 0;
    for (int i = ilo; i <= ihi; i++) {
        for (int j = MAX(jlo,i); j <= jhi; j++) {
            setflag[i][j] = 1;
            scale[i][j] = 1.0;
            count++;
        }
    }

    if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

double PairTadah::init_one(int i, int j)
{
    if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
    if (setflag[i][j] == 0) scale[i][j] = 1.0;
    scale[j][i] = scale[i][j];
    //return lt.cutoff(i,j);
    return lt.c.get<double>("RCUTMAX"); // TODO
}
void PairTadah::init_style()
{
    if (atom->tag_enable == 0)
        error->all(FLERR,"Pair style tadah/tadah requires atom IDs");

    neighbor->request(this,instance_me);

    if (lt.c.get<bool>("DIMER",0)) {
        neighbor->add_request(this, NeighConst::REQ_FULL);
        if (force->newton_pair)
            error->all(FLERR,"Pair style tadah/tadah requires newton pair off for this potential.");
    }
    else {

        if (lt.c.get<bool>("INITMB")) {
            if (lt.linear) {
                neighbor->add_request(this);
            }
            else {
                if (force->newton_pair == 0)
                    error->all(FLERR,"Pair style tadah/tadah requires newton pair on for this potential.");
                neighbor->add_request(this, NeighConst::REQ_FULL);
            }
        }
        else {
            if (lt.linear) {
                neighbor->add_request(this);
            }
            else {
                if (force->newton_pair == 0)
                    error->all(FLERR,"Pair style tadah/tadah requires newton pair on for this potential.");
                neighbor->add_request(this, NeighConst::REQ_FULL);
            }
        }
    }

}
/* ---------------------------------------------------------------------- */
int PairTadah::pack_reverse_comm(int n, int first, double *buf)
{
    if (lt.linear)
        return lt.pack_reverse_linear(n,first,buf);
    else
        return lt.pack_reverse_nonlinear(n,first,buf);
}
/* ---------------------------------------------------------------------- */
void PairTadah::unpack_reverse_comm(int n, int *list, double *buf)
{
    if (lt.linear)
        lt.unpack_reverse_linear(n,list,buf);
    else
        lt.unpack_reverse_nonlinear(n,list,buf);
}
/* ---------------------------------------------------------------------- */
int PairTadah::pack_forward_comm(int n, int *list, double *buf,
        int /*pbc_flag*/, int * /*pbc*/)
{
    if (lt.linear)
        return lt.pack_forward_linear(n,list,buf);
    else
        return lt.pack_forward_nonlinear(n,list,buf);
}
/* ---------------------------------------------------------------------- */
void PairTadah::unpack_forward_comm(int n, int first, double *buf)
{
    if (lt.linear)
        lt.unpack_forward_linear(n,first,buf);
    else
        lt.unpack_forward_nonlinear(n,first,buf);
}
void *PairTadah::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  return nullptr;
}
