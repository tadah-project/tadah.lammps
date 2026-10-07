// pair_coeff loads a potential whose weights fit pair_style tadah's
// descriptor, and refuses one whose weights do not.
//
// pair_style tadah builds the descriptor with DC_Selector: with BIAS true
// one one-hot intercept column per ATOM species, then the two-body and the
// many-body descriptors. Its predictors' dot products run over the weights,
// so a potential whose WEIGHTS were written for another layout gave
// silently wrong energies and forces. A multi-species BIAS potential trained
// by tadah before 1.6.0 has such weights: one intercept column shared by all
// species (Ta-Nb D2_LJ: -622 against -7.65 eV/atom).

#include "catch2/catch.hpp"

#include <tadah/core/context.h>
#include <tadah/core/core_user_config_parser.h>
#include <tadah/core/periodic_table.h>
#include <tadah/tadah.lammps/lammps_tadah.h>

#include <memory>
#include <string>
#include <vector>

using namespace tadah::core;
using tadah::LammpsTadah;

namespace {

/* ATOM Nb Ta, BIAS true, D2_BP (3 functions) for every pair: DSIZE 5 for
 * pair_style tadah (2 intercepts + 3). `model` lines, then `extra`. */
std::vector<std::string> nbta(const std::vector<std::string> &extra,
                              const std::string &model =
                                  "MODEL M_BLR BF_Linear") {
  std::vector<std::string> lines = {
      "ATOM Nb Ta",      "BIAS true",           "NORM false",
      "DIMER false 0 false", "EAMFAST false",   "INIT2B true",
      "TYPE2B D2_BP 3 3 * *", "CGRID2B 2.4 3.1 3.9", "SGRID2B 1.0 1.2 0.9",
      "RCTYPE2B Cut_Cos", "RCUT2B 5.0",         "INITMB false",
      model};
  lines.insert(lines.end(), extra.begin(), extra.end());
  return lines;
}

/* Load a potential the way pair_coeff * * <file> Ta Nb does (LAMMPS type 1
 * Ta, type 2 Nb). */
std::unique_ptr<LammpsTadah> load(const std::vector<std::string> &lines) {
  PeriodicTable::initialize();
  std::unique_ptr<LammpsTadah> lt(
      new LammpsTadah(CoreUserConfigParser().parse(lines)));
  std::vector<std::string> args = {"*", "*", "pot.tadah", "Ta", "Nb"};
  std::vector<char *> argv;
  for (auto &a : args) argv.push_back(&a[0]);
  lt->initialize(static_cast<int>(argv.size()), argv.data());
  return lt;
}

std::string refusal(const std::vector<std::string> &lines) {
  try {
    load(lines);
  } catch (const std::exception &e) {
    return e.what();
  }
  return "<accepted>";
}

const char *kLegacy = "one intercept column shared by all species";

}  // namespace

// Pins: pair_coeff refuses, with a message naming WEIGHTS, DSIZE, BIAS and
// ATOM and the retrain hint, an ATOM Nb Ta BIAS true potential with 4
// WEIGHTS (one intercept shared by the two species, the layout tadah wrote
// before 1.6.0) for its 5-entry descriptor; a 6-weight one without the
// hint; an M_KRR Kern_RBF potential whose basis vectors have 4 entries
// (the one-intercept layout) naming BASIS, with the hint; and a NORM true
// potential with 5 WEIGHTS whose NMEAN and NSTDEV have 4 values naming
// NMEAN.
// Before the fix: every one of them was accepted (the predictors ran over
// the weights: wrong energies and forces).
TEST_CASE("pair_coeff refuses a potential whose weights do not fit the "
          "descriptor of pair_style tadah",
          "[lammps_tadah][weights_layout]") {
  SECTION("4 WEIGHTS: one shared intercept") {
    const std::string m = refusal(nbta({"WEIGHTS -8.2 0.1 -0.05 0.02"}));
    INFO(m);
    CHECK(m.find("WEIGHTS has 4 values") != std::string::npos);
    CHECK(m.find("DSIZE 5") != std::string::npos);
    CHECK(m.find("BIAS true") != std::string::npos);
    CHECK(m.find("ATOM Nb Ta") != std::string::npos);
    CHECK(m.find(kLegacy) != std::string::npos);
    CHECK(m.find("retrain") != std::string::npos);
  }
  SECTION("6 WEIGHTS") {
    const std::string m =
        refusal(nbta({"WEIGHTS -8.2 -7.1 0.1 -0.05 0.02 0.3"}));
    INFO(m);
    CHECK(m.find("WEIGHTS has 6 values") != std::string::npos);
    CHECK(m.find(kLegacy) == std::string::npos);
  }
  SECTION("Kern_RBF with basis vectors of 4 entries") {
    const std::string m = refusal(nbta(
        {"MPARAMS 1.5", "SBASIS 3",
         "BASIS 1 0.1 0.2 0.3 1 0.4 0.5 0.6 1 0.7 0.8 0.9",
         "WEIGHTS 0.5 -0.2 0.1"},
        "MODEL M_KRR Kern_RBF 1"));
    INFO(m);
    CHECK(m.find("BASIS") != std::string::npos);
    CHECK(m.find("4 entries") != std::string::npos);
    CHECK(m.find(kLegacy) != std::string::npos);
  }
  SECTION("NORM true with 4 NMEAN and NSTDEV values") {
    std::vector<std::string> lines =
        nbta({"WEIGHTS -8.2 -7.1 0.1 -0.05 0.02", "NMEAN 1 0.5 0.4 0.3",
              "NSTDEV 0 0.2 0.3 0.4"});
    for (auto &l : lines)
      if (l == "NORM false") l = "NORM true";
    const std::string m = refusal(lines);
    INFO(m);
    CHECK(m.find("NMEAN has 4 values") != std::string::npos);
    CHECK(m.find(kLegacy) != std::string::npos);
  }
}

// Pins: a potential whose weights fit loads: ATOM Nb Ta, BIAS true with 5
// WEIGHTS (DSIZE 5 = the WEIGHTS count; LAMMPS type 1 = Ta has intercept
// column 1, type 2 = Nb column 0, and their per-atom constants are WEIGHTS
// 2 and 1), and an M_KRR Kern_RBF potential with 3 basis vectors of 5
// entries.
// Before the fix: passes as well (a pin: the refusal above must not refuse
// a potential in the per-species layout).
TEST_CASE("pair_coeff loads a potential whose weights fit the descriptor",
          "[lammps_tadah][weights_layout]") {
  std::unique_ptr<LammpsTadah> lt;
  REQUIRE_NOTHROW(lt = load(nbta({"WEIGHTS -8.2 -7.1 0.1 -0.05 0.02"})));
  CHECK(lt->S.descriptor_size() == 5);
  CHECK(lt->S.descriptor_size() == lt->c.size("WEIGHTS"));
  REQUIRE(lt->species_intercept_idx_per_type.size() == 2);
  CHECK(lt->species_intercept_idx_per_type[0] == 1);   // Ta
  CHECK(lt->species_intercept_idx_per_type[1] == 0);   // Nb
  CHECK(lt->bias_per_atom_type[0] == -7.1);
  CHECK(lt->bias_per_atom_type[1] == -8.2);
  CHECK_NOTHROW(load(nbta(
      {"MPARAMS 1.5", "SBASIS 3",
       "BASIS 1 0 0.1 0.2 0.3 0 1 0.4 0.5 0.6 1 0 0.7 0.8 0.9",
       "WEIGHTS 0.5 -0.2 0.1"},
      "MODEL M_KRR Kern_RBF 1")));
}
