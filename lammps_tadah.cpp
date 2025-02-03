#include <lammps_tadah.h>
#include <tadah/core/periodic_table.h>

using namespace TADAH;

LammpsTadah::LammpsTadah(int narg, char **arg):
  c(arg[2]),
  norm(c.get<bool>("NORM") && !(c("MODEL")[1]=="BF_Linear" || c("MODEL")[1]=="Kern_Linear")),
  cutoff_max_sq(pow(c.get<double>("RCUTMAX"),2)),
  cutoff_max(c.get<double>("RCUTMAX")),
  linear(c("MODEL")[1]=="BF_Linear" || c("MODEL")[1]=="Kern_Linear"),
  dimer(c.get<bool>("DIMER",0)),
  dimer_r(c.get<double>("DIMER",1)),
  dimer_bond_bool(c.get<bool>("DIMER",2)),
  init2b(c.get<bool>("INIT2B")),
  initmb(c.get<bool>("INITMB")),
  S(c),
  fb(CONFIG::factory<Function_Base,Config&>(
    c.get<std::string>("MODEL",1),c)),
  model(CONFIG::factory<M_MD_Base,Function_Base&,Config&>(
    c.get<std::string>("MODEL",0),*fb,c))
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
    for (j=0; j<c.size("ATOMS"); ++j) {
      if (c.get<std::string>("ATOMS", j)==symbol) {
        Z[i-3] = PeriodicTable::find_by_symbol(symbol).Z;
        break;
      }
    }
    if (j==c.size("ATOMS"))
      throw std::runtime_error("Symbol not in potential file: "+symbol);
  }

  // bias affects position of fidx for descriptors
  if (c.get<bool>("BIAS"))
    bias++;

  // add INTERNAL_KEYS, this is normally done by the descriptor_calc::common_constructor()
  // mind that it cannot be performed by the Config
  if (S.d2b->size()) {
    c.add("SIZE2B",S.d2b->size());
    dsize+=S.d2b->size();
    S.d2b->set_fidx(bias);
  }
  else {
    c.add("SIZE2B",0);
  }

  if (S.dmb->size()) {
    c.add("SIZEMB",S.dmb->size());
    dsize+=S.dmb->size();
    S.dmb->set_fidx(bias+S.d2b->size());
    rsize = S.dmb->rhoi_size()+ S.dmb->rhoip_size();
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

LammpsTadah::~LammpsTadah() {

  if (model)
    delete model;
  if (fb)
    delete fb;
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
