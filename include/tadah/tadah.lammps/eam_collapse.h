#pragma once
#ifndef TADAH_LAMMPS_EAM_COLLAPSE_H
#define TADAH_LAMMPS_EAM_COLLAPSE_H

/* Table collapse of a linear Tadah model.
 *
 * For a linear model the total energy is E_i = w . aed_i, so the 2b block
 * collapses to a scalar pair function phi(r) = sum_k w_k B_k(r) and an
 * EAM-form MB block to F(rho_bar) with rho_bar_i = sum_j rho(r_ij). The MD
 * hot path then needs only phi, phi', rho, rho' per pair (tabulated here as
 * cubic Hermite splines with analytic node derivatives, so energy and force
 * come from the same interpolant) plus the exact analytic F, F' per atom
 * (DM_Base::eam_embed).
 *
 * Two collapses are supported, selected by the model's descriptor set:
 *  - build_eam_collapse:    radial D2 + EAM-form DM (e.g. the `tadah refit`
 *    family D2_Knot5 + DM_REAM) -> phi + rho tables, analytic F.
 *  - build_pair2b_collapse: pure two-body (INITMB false) -> phi table only.
 *    This is the degenerate F==0, rho==0 case; no capability opt-in is
 *    needed because the phi sampling below is descriptor-agnostic.
 *
 * The pair tables are built by sampling the descriptors' OWN public calc
 * routines (zeroed-scratch trick for the 2b block), so they are node-exact
 * against the generic compute path for ANY radial D2 descriptor without
 * knowing its basis internals. */

#include <cmath>
#include <string>
#include <vector>

#include <tadah/core/core_types.h>
#include <tadah/models/dc_selector.h>

namespace tadah {

struct EamRadialTable {
  int nr = 0;
  double rmax = 0.0;       // table upper edge = descriptor cutoff
  double dr = 0.0;
  double inv_dr = 0.0;
  std::vector<double> v;   // node values   v_k   = f(k*dr)
  std::vector<double> d;   // node derivs   d_k   = f'(k*dr)

  /* Value and derivative from the same cubic Hermite segment, so that the
   * interpolated force is exactly minus the gradient of the interpolated
   * energy (matters for minimizer line searches). */
  inline void eval(const double r, double &val, double &der) const {
    double t = r * inv_dr;
    int k = static_cast<int>(t);
    if (k >= nr - 1) k = nr - 2;
    t -= k;
    const double p0 = v[k],   p1 = v[k + 1];
    const double m0 = d[k] * dr, m1 = d[k + 1] * dr;
    const double t2 = t * t;
    const double t3 = t2 * t;
    val = ( 2.0 * t3 - 3.0 * t2 + 1.0) * p0
        + (       t3 - 2.0 * t2 + t  ) * m0
        + (-2.0 * t3 + 3.0 * t2      ) * p1
        + (       t3 -       t2      ) * m1;
    der = (( 6.0 * t2 - 6.0 * t      ) * p0
        +  ( 3.0 * t2 - 4.0 * t + 1.0) * m0
        +  (-6.0 * t2 + 6.0 * t      ) * p1
        +  ( 3.0 * t2 - 2.0 * t      ) * m1) * inv_dr;
  }
  inline double eval_val(const double r) const {
    double val, der;
    eval(r, val, der);
    return val;
  }
};

struct EamCollapseTables {
  bool valid = false;
  bool has_mb = true;      // false for the pure-2b (phi-only) collapse
  EamRadialTable phi;      // pair function (per pair, total; halves tallied)
  EamRadialTable rho;      // scalar pair density (empty when !has_mb)
  tadah::core::t_type w;   // full model weight vector (for eam_embed)
  int Z = 0;               // single species atomic number
  double parity_err_phi = 0.0;  // max |table-analytic| mid-node self-check
  double parity_err_rho = 0.0;
};

/* phi(r) sampling: zeroed-scratch trick against the generic 2b routines.
 * (Zi,Zj) is passed explicitly so a future multi-species extension can key
 * one table per pair; today callers pass the single species twice. */
inline void sample_phi_table(const models::DC_Selector &S,
                             const tadah::core::t_type &w,
                             const int Zi, const int Zj, const int nr,
                             EamRadialTable &tab) {
  const std::size_t ndesc = S.descriptor_size();
  const double rc2b = S.d2b->get_rcut();
  tab.nr = nr;
  tab.rmax = rc2b;
  tab.dr = rc2b / (nr - 1);
  tab.inv_dr = 1.0 / tab.dr;
  tab.v.assign(nr, 0.0);
  tab.d.assign(nr, 0.0);
  tadah::core::aed_type aed(ndesc);          // zero-initialised
  tadah::core::fd_type fd(ndesc);
  for (int k = 1; k < nr; ++k) {             // r=0 kept at phi=0; never queried
    const double r = k * tab.dr;
    aed.set_zero();
    S.d2b->calc_aed(Zi, Zj, r, r * r, aed, 1.0);
    fd.set_zero();
    S.d2b->calc_dXijdri(Zi, Zj, r, r * r, fd, 1.0);
    double pv = 0.0, pd = 0.0;
    for (std::size_t c = 0; c < ndesc; ++c) {
      pv += w(c) * aed(c);
      pd += w(c) * fd(c, 0);
    }
    tab.v[k] = pv;
    tab.d[k] = pd;
  }
  /* r -> 0 head: extrapolate node 0 linearly from node 1 so the segment
   * [0,dr) is well-defined; pair distances that small are unphysical. */
  tab.d[0] = tab.d[1];
  tab.v[0] = tab.v[1] - tab.d[1] * tab.dr;
}

/* Mid-node phi parity self-check against the analytic sampling. */
inline double phi_parity_check(const models::DC_Selector &S,
                               const tadah::core::t_type &w,
                               const int Zi, const int Zj,
                               const EamRadialTable &tab) {
  const std::size_t ndesc = S.descriptor_size();
  tadah::core::aed_type aed(ndesc);
  double err = 0.0;
  for (int k = 1; k < tab.nr - 1; k += std::max(1, tab.nr / 64)) {
    const double r = (k + 0.5) * tab.dr;
    double tv, td;
    tab.eval(r, tv, td);
    aed.set_zero();
    S.d2b->calc_aed(Zi, Zj, r, r * r, aed, 1.0);
    double av = 0.0;
    for (std::size_t c = 0; c < ndesc; ++c) av += w(c) * aed(c);
    err = std::max(err, std::fabs(tv - av));
  }
  return err;
}

/* Build the pure two-body collapse (phi table only; INITMB false).
 * Returns false (with `why` set) when this S/model combination does not
 * collapse; caller falls back to the generic compute path. `w` must be the
 * full linear weight vector (no bias). */
inline bool build_pair2b_collapse(const models::DC_Selector &S,
                                  const tadah::core::t_type &w,
                                  const int Z, const int nr,
                                  EamCollapseTables &T, std::string &why) {
  T.valid = false;
  if (!S.d2b) { why = "missing d2b"; return false; }
  const std::size_t ndesc = S.descriptor_size();
  if (w.size() != ndesc) {
    why = "weight count " + std::to_string(w.size()) +
          " != descriptor size " + std::to_string(ndesc);
    return false;
  }
  const double rc2b = S.d2b->get_rcut();
  if (!(rc2b > 0.0)) { why = "non-positive cutoff"; return false; }

  T.w = w;
  T.Z = Z;
  T.has_mb = false;
  T.rho = EamRadialTable{};      // no density in the pure-2b collapse
  T.parity_err_rho = 0.0;

  sample_phi_table(S, w, Z, Z, nr, T.phi);
  T.parity_err_phi = phi_parity_check(S, w, Z, Z, T.phi);

  T.valid = true;
  return true;
}

/* Build the EAM collapse tables (radial D2 + EAM-form DM). Returns false
 * (with `why` set) when this S/model combination does not collapse; caller
 * falls back to the generic compute path. `w` must be the full linear
 * weight vector (no bias). */
inline bool build_eam_collapse(const models::DC_Selector &S,
                               const tadah::core::t_type &w,
                               const int Z, const int nr,
                               EamCollapseTables &T, std::string &why) {
  T.valid = false;
  if (!S.d2b || !S.dmb) { why = "missing d2b/dmb"; return false; }
  if (!S.dmb->eam_form()) {
    why = "DM descriptor " + S.dmb->label() + " is not EAM-form";
    return false;
  }
  const std::size_t ndesc = S.descriptor_size();
  if (w.size() != ndesc) {
    why = "weight count " + std::to_string(w.size()) +
          " != descriptor size " + std::to_string(ndesc);
    return false;
  }

  const double rc2b = S.d2b->get_rcut();
  const double rcmb = S.dmb->get_rcut();
  if (!(rc2b > 0.0) || !(rcmb > 0.0)) { why = "non-positive cutoff"; return false; }

  T.w = w;
  T.Z = Z;
  T.has_mb = true;

  sample_phi_table(S, w, Z, Z, nr, T.phi);

  /* ---- rho(r): exact scalar density from the DM descriptor ---- */
  T.rho.nr = nr;
  T.rho.rmax = rcmb;
  T.rho.dr = rcmb / (nr - 1);
  T.rho.inv_dr = 1.0 / T.rho.dr;
  T.rho.v.assign(nr, 0.0);
  T.rho.d.assign(nr, 0.0);
  for (int k = 1; k < nr; ++k) {
    const double r = k * T.rho.dr;
    S.dmb->eam_density(Z, Z, r, T.rho.v[k], T.rho.d[k]);
  }
  T.rho.d[0] = T.rho.d[1];
  T.rho.v[0] = T.rho.v[1] - T.rho.d[1] * T.rho.dr;

  /* ---- mid-node parity self-check against the analytic functions ---- */
  {
    tadah::core::aed_type aed(ndesc);
    tadah::core::fd_type fd(ndesc);
    double err_phi = 0.0, err_rho = 0.0;
    for (int k = 1; k < nr - 1; k += std::max(1, nr / 64)) {
      const double r = (k + 0.5) * T.phi.dr;
      double tv, td;
      T.phi.eval(r, tv, td);
      aed.set_zero();
      S.d2b->calc_aed(Z, Z, r, r * r, aed, 1.0);
      double av = 0.0;
      for (std::size_t c = 0; c < ndesc; ++c) av += w(c) * aed(c);
      err_phi = std::max(err_phi, std::fabs(tv - av));

      const double rr = (k + 0.5) * T.rho.dr;
      double rv, rd, arv, ard;
      T.rho.eval(rr, rv, rd);
      S.dmb->eam_density(Z, Z, rr, arv, ard);
      err_rho = std::max(err_rho, std::fabs(rv - arv));
    }
    T.parity_err_phi = err_phi;
    T.parity_err_rho = err_rho;
    (void)fd;
  }

  T.valid = true;
  return true;
}

} // namespace tadah
#endif
