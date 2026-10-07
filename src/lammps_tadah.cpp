#include <tadah/tadah.lammps/lammps_tadah.h>
#include <tadah/core/core_user_config_parser.h>
#include <tadah/core/intercept_layout.h>
#include <tadah/core/periodic_table.h>
#include <tadah/core/registry.h>
#include <tadah/models/descriptors/dimer_element_pairs.h>
#include <tadah/models/functions/basis_functions/bf_polynomial2.h>
#include <tadah/models/functions/kernels/kern_base.h>
#include <tadah/models/functions/weights_fit.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

using namespace tadah;

namespace {

/* With BIAS true on several species pair_style tadah's descriptor has one
 * bias component per species (S.bias of them), while tadah (MLIP) trains
 * with one bias component shared by the species, S.bias - 1 components
 * fewer: the descriptor size of a potential trained so, 0 when the two
 * layouts agree (BIAS false, or one species). */
std::size_t shared_bias_descriptor_size(const LammpsTadah& self)
{
    const std::size_t nb = static_cast<std::size_t>(self.S.bias);
    return nb < 2 ? 0 : self.S.descriptor_size() - (nb - 1);
}

/* The refusal of a non-linear potential laid out for that descriptor. */
std::string shared_bias_refusal(const LammpsTadah& self,
                                const std::string& model)
{
    return model + ": the potential has BIAS true and " +
           std::to_string(self.S.bias) + " species (ATOM). It is laid out "
           "for a descriptor of " +
           std::to_string(shared_bias_descriptor_size(self)) + " components, "
           "with one bias component shared by the species, as tadah trains "
           "it, but pair_style tadah's descriptor has one bias component per "
           "species, " + std::to_string(self.S.descriptor_size()) +
           " components. pair_style tadah cannot evaluate a non-linear "
           "potential trained with BIAS true on several species; one trained "
           "with BIAS false runs.";
}

}  // namespace

/*-------------------------------------------*/
/* helper that does the heavy initialisation */
/*-------------------------------------------*/
void build_sub_objects(LammpsTadah& self)
{
    self.norm = self.c.get<bool>("NORM") &&
                !( self.c("MODEL")[1] == "BF_Linear" ||
                   self.c("MODEL")[1] == "Kern_Linear" );

    self.cutoff_max     = self.c.get<double>("RCUTMAX");
    self.cutoff_max_sq  = std::pow(self.cutoff_max, 2);

    self.linear = self.c("MODEL")[1] == "BF_Linear" ||
                  self.c("MODEL")[1] == "Kern_Linear";

    self.dimer           = self.c.get<bool  >( "DIMER", 0 );
    self.dimer_r         = self.c.get<double>( "DIMER", 1 );
    self.dimer_bond_bool = self.c.get<bool  >( "DIMER", 2 );

    self.init2b = self.c.get<bool>("INIT2B");
    self.initmb = self.c.get<bool>("INITMB");

    /* PairTadah::compute_dimers evaluates every atom as Z=1: a descriptor
     * declared for specific element pairs would apply to no pair. */
    models::require_dimer_element_pairs(self.c, *self.S.d2b, *self.S.dmb);

    /* ---- create Function_Base ---- */
    auto fb_ptr = core::factory<models::Function_Base, core::Context&>(
                      self.c.get<std::string>("MODEL", 1), self.c);
    if (!fb_ptr)
        throw std::runtime_error("cannot build Function_Base for model "
                                 + self.c.get<std::string>("MODEL", 1));
    self.fb.reset(fb_ptr);

    /* ---- create MD model ---- */
    auto model_ptr = tadah::core::factory<
                         md::M_MD_Base,
                         models::Function_Base&, core::Context&>(
                         self.c.get<std::string>("MODEL", 0),
                         *self.fb, self.c);
    if (!model_ptr)
        throw std::runtime_error("cannot build MD model "
                                 + self.c.get<std::string>("MODEL", 0));
    self.model.reset(model_ptr);

    /* ---- a non-linear potential fits the descriptor ----
     * Only here is the descriptor size known. A non-linear kernel compares
     * each descriptor with its basis vectors (SBASIS vectors stored in
     * BASIS; M_KRR_Core refuses WEIGHTS without a basis or with another
     * number of vectors): vectors of another length read past one of the
     * two in every K(b, x). BF_Polynomial2 has one weight per product
     * x_i x_j, i <= j, of the D descriptor components, D(D+1)/2: other
     * WEIGHTS read past them, or leave some out. */
    const std::size_t dsize = self.S.descriptor_size();
    const std::string model = self.c.get<std::string>("MODEL", 0) + " " +
                              self.fb->get_label();
    if (const auto* kb = dynamic_cast<const models::Kern_Base*>(self.fb.get())) {
        const core::Matrix& b = kb->get_basis();
        if (b.cols() > 0 && b.rows() != dsize) {
            if (b.rows() == shared_bias_descriptor_size(self))
                throw std::runtime_error(shared_bias_refusal(self, model));
            throw std::runtime_error(
                model + ": BASIS holds SBASIS " +
                std::to_string(b.cols()) + " basis vectors of " +
                std::to_string(b.rows()) + " values each, but the descriptor "
                "of the potential has " + std::to_string(dsize) +
                " components: BASIS must have SBASIS x " +
                std::to_string(dsize) + " values.");
        }
    }
    if (dynamic_cast<const models::BF_Polynomial2*>(self.fb.get())) {
        const std::size_t nw =
            self.c.contains("WEIGHTS") ? self.c.size("WEIGHTS") : 0;
        if (nw != dsize * (dsize + 1) / 2) {
            const std::size_t d = shared_bias_descriptor_size(self);
            if (d > 0 && nw == d * (d + 1) / 2)
                throw std::runtime_error(shared_bias_refusal(self, model));
            throw std::runtime_error(
                model + ": the potential has " + std::to_string(nw) +
                " WEIGHTS, but the descriptor of the potential has " +
                std::to_string(dsize) + " components: BF_Polynomial2 has one "
                "weight per product x_i x_j (i <= j) of them, " +
                std::to_string(dsize * (dsize + 1) / 2) + " weights.");
        }
    }
}

/*-------------------------------------------*/
/* main constructor – takes Context by move  */
/*-------------------------------------------*/
LammpsTadah::LammpsTadah(core::Context&& ctx)
    : c(std::move(ctx)),
      embf0_declared(c.has_stored("EMBF0")),   // before S materialises it
      S(c)                     // DC_Selector may need the Context early
{
    build_sub_objects(*this);
}

/*-------------------------------------------*/
/* copy-in constructor                       */
/*-------------------------------------------*/
LammpsTadah::LammpsTadah(const core::Context& ctx)
    : c(ctx),                 // one copy
      embf0_declared(c.has_stored("EMBF0")),   // before S materialises it
      S(c)
{
    build_sub_objects(*this);
}

/*-------------------------------------------*/
/* command-line constructor                  */
/*-------------------------------------------*/
LammpsTadah::LammpsTadah(int argc, char** argv)
    : LammpsTadah(
          core::CoreUserConfigParser()
              .parse(core::to_string(argv[2])) )
{
}


void LammpsTadah::initialize(int narg, char **arg)
{
  // Map LAMMPS atom types to weight factors
  // e.g. pair_coaeff * * pot.tadah Ti Nb NULL
  // watom.resize(narg-3); // this includes NULL values
  Z.resize(narg-3); // this includes NULL values
                        // which want be set
  for (int i=3; i<narg; ++i) {
    std::string symbol=arg[i];
    if (symbol=="NULL") continue;
    size_t j=0;
    for (j=0; j<c.size("ATOM"); ++j) {
      if (c.get<std::string>("ATOM", j)==symbol) {
        Z[i-3] = core::PeriodicTable::find_by_symbol(symbol).Z;
        break;
      }
    }
    if (j==c.size("ATOM"))
      throw std::runtime_error("Symbol not in potential file: "+symbol);
  }

  c.addInternal("DSIZE",S.descriptor_size());
  bias = S.bias;
  rsize = S.rsize; // size of the density array (rhoi + rhoip)

  // The predictors' dot products run over the weights: refuse WEIGHTS,
  // NMEAN/NSTDEV (NORM true) or basis vectors that do not fit this
  // descriptor (e.g. a multi-species BIAS potential trained before 1.6.0,
  // with one intercept column shared by all species) instead of computing
  // with misaligned weights. pair_coeff reports it through error->all.
  models::require_weights_fit(c, S.descriptor_size(), *fb);

  // Precompute per-LAMMPS-atom-type intercept lookup tables. With per-species
  // BIAS, every atom of LAMMPS type t contributes a fixed constant
  // bias_per_atom_type[t] to its energy; the descriptor's intercept column for
  // that species is species_intercept_idx_per_type[t]. Probing the model with
  // a one-hot AED is identical to looking up w[k] in the linear case but works
  // for any model that implements epredict.
  const size_t ntypes = Z.size();
  species_intercept_idx_per_type.assign(ntypes, -1);
  bias_per_atom_type.assign(ntypes, 0.0);

  if (bias > 0) {
    const auto col = core::intercept_columns_by_Z(c);
    core::aed_type probe(S.descriptor_size());
    for (size_t t = 0; t < ntypes; ++t) {
      const int z = Z[t];
      if (z == 0) continue; // NULL slot — leave defaults
      const int k = col[z];
      if (k < 0) continue;
      species_intercept_idx_per_type[t] = k;
      probe.set_zero();
      probe(k) = 1.0;
      bias_per_atom_type[t] = model->epredict(probe);
    }
  }
}


int LammpsTadah::pack_reverse_linear(int n, int first, double *buf)
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
void LammpsTadah::unpack_reverse_linear(int n, int *list, double *buf)
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
int LammpsTadah::pack_forward_linear(int n, int *list, double *buf)
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
void LammpsTadah::unpack_forward_linear(int n, int first, double *buf) 
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    for (size_t x = S.dmb->rhoi_size(); x < S.dmb->rhoi_size()+S.dmb->rhoip_size(); x++)
      rhos(x,i) = buf[m++];
  }
}
// NONLINEAR
int LammpsTadah::pack_reverse_nonlinear(int n, int first, double *buf)
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
void LammpsTadah::unpack_reverse_nonlinear(int n, int *list, double *buf)
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
int LammpsTadah::pack_forward_nonlinear(int n, int *list, double *buf)
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
void LammpsTadah::unpack_forward_nonlinear(int n, int first, double *buf) 
{
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
