#ifndef LAMMPS_TADAH_H
#define LAMMPS_TADAH_H

#include <tadah/core/config.h>
#include <tadah/core/core_types.h>
#include <tadah/core/registry.h>
#include <tadah/core/periodic_table.h>
#include <tadah/models/descriptors/d_all.h>
#include <tadah/models/cutoffs.h>
#include <tadah/models/functions/function_base.h>
#include <tadah/models/dc_selector.h>
#include <tadah/md/m_md_base.h>

#include <iostream>
#include <string>

namespace TADAH {
class LammpsTadah {

public:
  Config c;
  bool norm;
  size_t bias=0;
  size_t dsize=0; // descriptor dim inc bias
  size_t rsize=0; // rho dim inc derivative
  double cutoff_max_sq;
  double cutoff_max;
  bool linear;
  bool dimer;
  double dimer_r;
  bool dimer_bond_bool;

  bool init2b;
  bool initmb;
  std::vector<int> Z;

  DC_Selector S;
  Function_Base *fb=nullptr;
  M_MD_Base *model=nullptr;

  aeds_type2 aeds;
  rhos_type rhos;

  LammpsTadah(int narg, char **arg);
  ~LammpsTadah();
  LammpsTadah(const LammpsTadah&) = delete;
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

};
}
#endif
