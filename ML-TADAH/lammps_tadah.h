#ifndef LAMMPS_TADAH_H
#define LAMMPS_TADAH_H

#include <iostream>
#include <string>
#include "../libs/Eigen/Dense"
#include "../config/config.h"
#include "../utils/typedefs.h"
//#include "../utils/utils.h"
#include "../utils/registry.h"
#include "../utils/periodic_table.h"
#include "../descriptors/d_all.h"
#include "../cutoffs/cutoffs.h"
#include "../models/m_md_base.h"
#include "../functions/function_base.h"
#include "../dc_selector.h"

/**
 * LAMMPS INTERFACE
 */
namespace TADAH {
class LammpsTadah {

    public:
        Config c;
        bool norm;
        bool bias;
        double cutoff_max_sq;
        double cutoff_max;
        rhos_type rhos;
        aeds_type aeds;
        bool linear;
        bool dimer;
        bool dimer_bond_bool;
        double dimer_r;

        bool init2b;
        bool initmb;

        DC_Selector S;
        std::vector<double> weights;
        M_MD_Base *model=nullptr;
        Function_Base *fb=nullptr;

        // Const Assignment operator
        LammpsTadah& operator=(const LammpsTadah& lt)
        {
            c=lt.c;
            norm = lt.norm;
            bias=lt.bias;
            cutoff_max_sq=lt.cutoff_max_sq;
            cutoff_max=lt.cutoff_max;
            linear=lt.linear;
            dimer=lt.dimer;
            dimer_bond_bool=lt.dimer_bond_bool;
            dimer_r=lt.dimer_r;
            init2b=lt.init2b;
            initmb=lt.initmb;
            weights=lt.weights;
            S=lt.S;
            // We need deep copy for fb and model.
            // For now we just init them the same way
            // as constructor does
            fb = factory<Function_Base,Config&>(
                    c.get<std::string>("MODEL",1),c);
            model = factory<M_MD_Base,Function_Base&,Config&>(
                    c.get<std::string>("MODEL",0),*fb,c);
            aeds=lt.aeds;
            rhos=lt.rhos;
            return *this;
        }
        // Assignment operator
        LammpsTadah& operator=(LammpsTadah& lt)
        {
            std::swap(c,lt.c);
            std::swap(norm , lt.norm);
            std::swap(bias,lt.bias);
            std::swap(cutoff_max_sq,lt.cutoff_max_sq);
            std::swap(cutoff_max,lt.cutoff_max);
            std::swap(linear,lt.linear);
            std::swap(dimer,lt.dimer);
            std::swap(dimer_bond_bool,lt.dimer_bond_bool);
            std::swap(dimer_r,lt.dimer_r);
            std::swap(init2b,lt.init2b);
            std::swap(initmb,lt.initmb);
            std::swap(fb,lt.fb);
            std::swap(model,lt.model);
            std::swap(weights,lt.weights);
            S=lt.S;
            aeds=lt.aeds;
            rhos=lt.rhos;
            return *this;
        }
        // Constructor:
        LammpsTadah():
            model(nullptr),
            fb(nullptr)
    {};
        LammpsTadah(int narg, char **arg):
            c(arg[2]),
            norm(c.get<bool>("NORM") && !(c("MODEL")[1]=="BF_Linear" || c("MODEL")[1]=="Kern_Linear")),
            bias(c.get<bool>("BIAS")),
            cutoff_max_sq(pow(c.get<double>("RCUTMAX"),2)),
            cutoff_max(c.get<double>("RCUTMAX")),
            linear(c("MODEL")[1]=="BF_Linear" || c("MODEL")[1]=="Kern_Linear"),
            dimer(c.get<bool>("DIMER",0)),
            dimer_bond_bool(c.get<bool>("DIMER",2)),
            dimer_r(c.get<double>("DIMER",1)),
            init2b(c.get<bool>("INIT2B")),
            initmb(c.get<bool>("INITMB")),
            S(c),
            model(nullptr),
            fb(nullptr)
    {
        // check does pot file exists
        std::ifstream f(arg[2]);
        if(!f.good())
            throw std::runtime_error("Pot file does not exists.\n");

        // Map LAMMPS atom types to weight factors
        // e.g. pair_coaeff * * pot.tadah Ti Nb
        PeriodicTable pt;
        weights.resize(narg-3);
        for (int i=3; i<narg; ++i) {
            weights[i-3] = pt.find_by_symbol(arg[i]).Z;
        }

        // model does not need to know what calculators are being used
        if (c.size("MODEL")!=2)
            throw std::runtime_error("MODEL keyword requires two values\n\
                    e.g. MODEL M_BLR BF_Linear\n");

        fb = factory<Function_Base,Config&>(
                c.get<std::string>("MODEL",1),c);
        model = factory<M_MD_Base,Function_Base&,Config&>(
                c.get<std::string>("MODEL",0),*fb,c);

        // bias affects position of fidx for descriptors
        size_t bias=0;
        if (c.get<bool>("BIAS"))
            bias++;

        // add INTERNAL_KEYS
        size_t dsize=0;
        if (S.d2b->size()) {
            c.add("SIZE2B",S.d2b->size());
            dsize+=S.d2b->size();
        }
        else {
            c.add("SIZE2B",0);
        }

        if (S.dmb->size()) {
            c.add("SIZEMB",S.dmb->size());
            dsize+=S.dmb->size();
        }
        else {
            c.add("SIZEMB",0);
        }

        if (!dsize)
            throw std::runtime_error("The descriptor size is 0, check your config.");

        if (c.get<bool>("BIAS"))
            dsize++;

        c.add("DSIZE",dsize);
    }

        ~LammpsTadah() {

            if (model)
                delete model;
            if (fb)
                delete fb;
        }


int pack_reverse_linear(int n, int first, double *buf)
{
    int m = 0;
    int last = first + n;
    for (int i = first; i < last; i++) {
        for (size_t x = 0; x < S.dmb->rhoi_size(); x++)
            buf[m++] = rhos(x,i);
        // pack two-body only, as many-body will be computed from rho
        // which is being packed hence we do now comm mb-aeds
        for (size_t x = bias; x < S.d2b->size()+bias; x++)
            buf[m++] = aeds(x,i);
    }
    return m;
}
void unpack_reverse_linear(int n, int *list, double *buf)
{
    int m = 0;
    for (int i = 0; i < n; i++) {
        int j = list[i];
        for (size_t x = 0; x < S.dmb->rhoi_size(); x++)
            rhos(x,j) += buf[m++];
        for (size_t x = bias; x < S.d2b->size()+bias; x++)
            aeds(x,j) += buf[m++];
    }
}
// Pack the derivative of the embedding func
int pack_forward_linear(int n, int *list, double *buf)
{
    int m = 0;
    for (int i = 0; i < n; i++) {
        int j = list[i];
        for (size_t x = S.dmb->rhoi_size(); x < S.dmb->rhoi_size()+S.dmb->rhoip_size(); x++)
            buf[m++] = rhos(x,j);
    }
    return m;
}
// Unpack the derivative of the embedding func
void unpack_forward_linear(int n, int first, double *buf) {
    int m = 0;
    int last = first + n;
    for (int i = first; i < last; i++) {
        for (size_t x = S.dmb->rhoi_size(); x < S.dmb->rhoi_size()+S.dmb->rhoip_size(); x++)
            rhos(x,i) = buf[m++];
    }
}
// NONLINEAR
int pack_reverse_nonlinear(int n, int first, double *buf)
{
    int i,m,last;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
        for (size_t x = 0; x < S.dmb->rhoi_size(); x++)
            buf[m++] = rhos(x,i);
        // pack two-body only, as many-body will be computed from rho
        // which is packed above hence we do now comm mb-aeds
        for (size_t x = bias; x < S.d2b->size()+bias; x++)
            buf[m++] = aeds(x,i);
    }
    return m;
}
void unpack_reverse_nonlinear(int n, int *list, double *buf)
{
    int i,j,m;

    m = 0;
    for (i = 0; i < n; i++) {
        j = list[i];
        for (size_t x = 0; x < S.dmb->rhoi_size(); x++)
            rhos(x,j) += buf[m++];
        for (size_t x = bias; x < S.d2b->size()+bias; x++)
            aeds(x,j) += buf[m++];
    }
}
int pack_forward_nonlinear(int n, int *list, double *buf)
{
    int i,j,m;

    m = 0;
    for (i = 0; i < n; i++) {
        j = list[i];
        for (size_t x = S.dmb->rhoi_size(); x < S.dmb->rhoi_size()+S.dmb->rhoip_size(); x++)
            buf[m++] = rhos(x,j);
        for (size_t x = bias; x < S.d2b->size()+S.dmb->size()+bias; x++)
            buf[m++] = aeds(x,j);
    }
    return m;
}
void unpack_forward_nonlinear(int n, int first, double *buf) {
    int i,m,last;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
        for (size_t x = S.dmb->rhoi_size(); x < S.dmb->rhoi_size()+S.dmb->rhoip_size(); x++)
            rhos(x,i) = buf[m++];
        for (size_t x = bias; x < S.d2b->size()+S.dmb->size()+bias; x++)
            aeds(x,i) = buf[m++];
    }
}

};
}
#endif

