#pragma once
#ifndef LAMMPS_TADAH_H
#define LAMMPS_TADAH_H

#include <memory>
#include <vector>
#include <tadah/core/context.h>
#include <tadah/models/dc_selector.h>
#include <tadah/models/functions/function_base.h>
#include <tadah/md/m_md_base.h>

namespace tadah {

class LammpsTadah
{
public:
    explicit LammpsTadah(const  core::Context& ctx);   // copy
    explicit LammpsTadah(      core::Context&& ctx);   // move
    LammpsTadah(int argc, char** argv);                // from CLI

    void initialize(int argc, char** argv);

    /* … pack/unpack functions unchanged … */

    ~LammpsTadah() = default;                          // unique_ptr ⇒ OK
    LammpsTadah(const LammpsTadah&)            = delete;
    LammpsTadah& operator=(const LammpsTadah&) = delete;

  int pack_reverse_linear(int n, int first, double *buf);
  void unpack_reverse_linear(int n, int *list, double *buf);
  // Pack the derivative of the embedding func
  int pack_forward_linear(int n, int *list, double *buf);
  // Unpack the derivative of the embedding func
  void unpack_forward_linear(int n, int first, double *buf);
  // NONLINEAR
  int pack_reverse_nonlinear(int n, int first, double *buf);
  void unpack_reverse_nonlinear(int n, int *list, double *buf);
  int pack_forward_nonlinear(int n, int *list, double *buf);
  void unpack_forward_nonlinear(int n, int first, double *buf);

public:
  bool pot_is_loaded=false;  // true if potential is loaded

    /* ------------ plain data ------------ */
    core::Context  c;

    bool   norm   {};
    size_t bias   {0};
    size_t dsize  {0};
    size_t rsize  {0};

    double cutoff_max_sq {};
    double cutoff_max    {};

    bool   linear {};
    bool   dimer  {};
    double dimer_r{};
    bool   dimer_bond_bool{};

    bool init2b{};
    bool initmb{};

    std::vector<int> Z;

    /* ------------ large sub-objects -------- */
    models::DC_Selector                         S;
    std::unique_ptr<models::Function_Base>      fb;
    std::unique_ptr<md::M_MD_Base>              model;

    core::aeds_type2 aeds;
    core::rhos_type  rhos;
};
} // namespace tadah
#endif
