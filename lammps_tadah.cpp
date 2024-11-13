#include "lammps_tadah.h"

using namespace TADAH;

// Const Assignment operator
LammpsTadah &LammpsTadah::operator=(const LammpsTadah& lt)
{
  c=lt.c;
  norm=lt.norm;
  bias=lt.bias;
  dsize=lt.dsize;
  cutoff_max_sq=lt.cutoff_max_sq;
  cutoff_max=lt.cutoff_max;
  linear=lt.linear;
  dimer=lt.dimer;
  dimer_r=lt.dimer_r;
  dimer_bond_bool=lt.dimer_bond_bool;
  init2b=lt.init2b;
  initmb=lt.initmb;
  watom=lt.watom;
  S=lt.S;
  // We need deep copy for fb and model.
  // For now we just init them the same way
  // as constructor does
  fb = CONFIG::factory<Function_Base,Config&>(
      c.get<std::string>("MODEL",1),c);
  model = CONFIG::factory<M_MD_Base,Function_Base&,Config&>(
      c.get<std::string>("MODEL",0),*fb,c);
  aeds=lt.aeds;
  rhos=lt.rhos;
  return *this;
}
LammpsTadah &LammpsTadah::operator=(LammpsTadah& lt)
{
  std::swap(c,lt.c);
  std::swap(norm , lt.norm);
  std::swap(bias,lt.bias);
  std::swap(dsize,lt.dsize);
  std::swap(cutoff_max_sq,lt.cutoff_max_sq);
  std::swap(cutoff_max,lt.cutoff_max);
  std::swap(linear,lt.linear);
  std::swap(dimer,lt.dimer);
  std::swap(dimer_r,lt.dimer_r);
  std::swap(dimer_bond_bool,lt.dimer_bond_bool);
  std::swap(init2b,lt.init2b);
  std::swap(initmb,lt.initmb);
  std::swap(watom,lt.watom);
  S=lt.S;
  std::swap(fb,lt.fb);
  std::swap(model,lt.model);
  aeds=lt.aeds;
  rhos=lt.rhos;
  return *this;
}
LammpsTadah::LammpsTadah():
  fb(nullptr),
  model(nullptr)
{}

LammpsTadah::LammpsTadah(int narg, char **arg):
  c(arg[2]),
  norm(c.get<bool>("NORM") && !(c("MODEL")[1]=="BF_Linear" || c("MODEL")[1]=="Kern_Linear")),
  //bias(c.get<bool>("BIAS")),
  cutoff_max_sq(pow(c.get<double>("RCUTMAX"),2)),
  cutoff_max(c.get<double>("RCUTMAX")),
  linear(c("MODEL")[1]=="BF_Linear" || c("MODEL")[1]=="Kern_Linear"),
  dimer(c.get<bool>("DIMER",0)),
  dimer_r(c.get<double>("DIMER",1)),
  dimer_bond_bool(c.get<bool>("DIMER",2)),
  init2b(c.get<bool>("INIT2B")),
  initmb(c.get<bool>("INITMB")),
  S(c),
  fb(nullptr),
  model(nullptr)
{
  // check does pot file exists
  std::ifstream f(arg[2]);
  if(!f.good())
    throw std::runtime_error("Pot file does not exists.\n");

  // Map LAMMPS atom types to weight factors
  // e.g. pair_coaeff * * pot.tadah Ti Nb NULL
  watom.resize(narg-3); // this includes NULL values
                        // which want be set

  for (int i=3; i<narg; ++i) {
    std::string symbol=arg[i];
    if (symbol=="NULL") continue;
    size_t j=0;
    for (j=0; j<c.size("ATOMS"); ++j) {
      if (c.get<std::string>("ATOMS", j)==symbol) {
        watom[i-3] = c.get<double>("WATOMS",j);
        break;
      }
    }
    if (j==c.size("ATOMS"))
      throw std::runtime_error("Symbol not in potential file: "+symbol);
  }

  fb = CONFIG::factory<Function_Base,Config&>(
      c.get<std::string>("MODEL",1),c);
  model = CONFIG::factory<M_MD_Base,Function_Base&,Config&>(
      c.get<std::string>("MODEL",0),*fb,c);

  // bias affects position of fidx for descriptors
  if (c.get<bool>("BIAS"))
    bias++;

  // add INTERNAL_KEYS, this is normally done by the descriptor_calc::common_constructor()
  // mind that it cannot be performed by the Config
  if (S.d2b->size()) {
    c.add("SIZE2B",S.d2b->size());
    dsize+=S.d2b->size();
    S.d2b->fidx=bias;
  }
  else {
    c.add("SIZE2B",0);
  }

  if (S.dmb->size()) {
    c.add("SIZEMB",S.dmb->size());
    dsize+=S.dmb->size();
    S.dmb->fidx=bias+S.d2b->size();
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
      buf[m++] = rhos(i,x);
    // pack two-body only, as many-body will be computed from rho
    // which is being packed hence we do now comm mb-aeds
    for (size_t x = bias; x < S.d2b->size()+bias; x++)
      buf[m++] = aeds(i,x);
  }
  return m;
}
void LammpsTadah::unpack_reverse_linear(int n, int *list, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    for (size_t x = 0; x < S.dmb->rhoi_size(); x++)
      rhos(j,x) += buf[m++];
    for (size_t x = bias; x < S.d2b->size()+bias; x++)
      aeds(j,x) += buf[m++];
  }
}
// Pack the derivative of the embedding func
int LammpsTadah::pack_forward_linear(int n, int *list, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    for (size_t x = S.dmb->rhoi_size(); x < S.dmb->rhoi_size()+S.dmb->rhoip_size(); x++)
      buf[m++] = rhos(j,x);
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
      rhos(i,x) = buf[m++];
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
      buf[m++] = rhos(i,x);
    // pack two-body only, as many-body will be computed from rho
    // which is packed above hence we do now comm mb-aeds
    for (size_t x = bias; x < S.d2b->size()+bias; x++)
      buf[m++] = aeds(i,x);
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
      rhos(j,x) += buf[m++];
    for (size_t x = bias; x < S.d2b->size()+bias; x++)
      aeds(j,x) += buf[m++];
  }
}
int LammpsTadah::pack_forward_nonlinear(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (size_t x = S.dmb->rhoi_size(); x < S.dmb->rhoi_size()+S.dmb->rhoip_size(); x++)
      buf[m++] = rhos(j,x);
    for (size_t x = bias; x < S.d2b->size()+S.dmb->size()+bias; x++)
      buf[m++] = aeds(j,x);
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
      rhos(i,x) = buf[m++];
    for (size_t x = bias; x < S.d2b->size()+S.dmb->size()+bias; x++)
      aeds(i,x) = buf[m++];
  }
}
