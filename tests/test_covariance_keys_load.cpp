// pair_style tadah loads (tadah::LammpsTadah, the object pair_coeff builds
// from the potential file) a potential that carries a weight covariance:
// WCOV and WNOISE, which `tadah train` writes for an OALGO 3 fit, and
// SIGMA, which Tadah! 1.0.0 to 1.5.0 and the 1.6.0 development builds wrote
// instead (a deprecated key now: one notice, not stored). pair_style tadah
// does not use the covariance.
//
// The potential is the D2_LJ pair potential of
// tests_data/pots/ta_lj_cutcos.tadah.

#include "catch2/catch.hpp"

#include <tadah/core/key_notice.h>
#include <tadah/core/periodic_table.h>
#include <tadah/tadah.lammps/lammps_tadah.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

const char *kPairPotential =
    "ATOM     Ta\n"
    "BIAS     false\n"
    "DIMER     false     0     false\n"
    "EWEIGHT     1.0\n"
    "FWEIGHT     1.0\n"
    "INIT2B     true\n"
    "INITMB     false\n"
    "LAMBDA     0\n"
    "MODEL     M_KRR     Kern_Linear     1\n"
    "NORM     false\n"
    "OALGO     3\n"
    "RCTYPE2B     Cut_Cos\n"
    "RCUT2B     6.5\n"
    "SWEIGHT     1.0\n"
    "TYPE2B     D2_LJ     Ta     Ta\n"
    "WATOM     73\n"
    "WEIGHTS     14.204213349318765225     2604.396965191186154698\n";

/* A potential file unique to this process, removed with the object. */
struct PotFile {
  std::string path;
  explicit PotFile(const std::string &extra) {
    path = "covariance_keys_" + std::to_string(::getpid()) + ".tadah";
    std::ofstream(path) << kPairPotential << extra;
  }
  ~PotFile() { std::remove(path.c_str()); }
};

/* pair_coeff * * <file> Ta: LammpsTadah from the file, then the type map. */
std::unique_ptr<tadah::LammpsTadah> load(const std::string &path) {
  std::vector<std::string> args = {"*", "*", path, "Ta"};
  std::vector<char *> argv;
  for (auto &a : args) argv.push_back(&a[0]);
  std::unique_ptr<tadah::LammpsTadah> lt(
      new tadah::LammpsTadah(static_cast<int>(argv.size()), argv.data()));
  lt->initialize(static_cast<int>(argv.size()), argv.data());
  return lt;
}

std::string load_error(const std::string &path) {
  try {
    load(path);
  } catch (const std::exception &e) {
    return e.what();
  }
  return "";
}

}  // namespace

/* Pins: a potential with WCOV (N = 2, then the 2 x 2 matrix) and WNOISE (E
 * and F) loads as pair_coeff loads it, and its Context keeps both keys
 * (5 and 4 values).
 * Before the fix: the load threw "Invalid key: WCOV" (pair_coeff failed). */
TEST_CASE("pair_style tadah loads a potential with WCOV and WNOISE",
          "[lammps_tadah][weight_covariance]") {
  tadah::core::PeriodicTable::initialize();
  const PotFile pot("WCOV     2     0.5 0.1 0.1 0.25\n"
                    "WNOISE     E 0.04     F 0.0025\n");
  std::unique_ptr<tadah::LammpsTadah> lt;
  std::string msg;
  try {
    lt = load(pot.path);
  } catch (const std::exception &e) {
    msg = e.what();
  }
  INFO("load: " << msg);
  REQUIRE(lt != nullptr);
  CHECK(lt->c.size("WCOV") == 5);
  CHECK(lt->c.size("WNOISE") == 4);
}

/* Pins: a potential with a SIGMA line (an OALGO 3 potential of Tadah!
 * 1.0.0 to 1.5.0 or a 1.6.0 development build) loads as pair_coeff loads
 * it, writes the deprecation
 * notice of SIGMA ("WARNING: Deprecated key SIGMA is ignored: ...") once,
 * on standard error (pair_style tadah routes it to LAMMPS' warning, see
 * HPO's test_pair_tadah_notices.cpp), and its Context has no SIGMA.
 * Before the fix: SIGMA was stored, with no notice. */
TEST_CASE("pair_style tadah loads a potential with SIGMA and ignores it",
          "[lammps_tadah][weight_covariance]") {
  tadah::core::PeriodicTable::initialize();
  const PotFile pot("SIGMA     2     1 0 0 1\n");
  const std::size_t expected =
      tadah::core::ignored_key_reported("SIGMA") ? 0 : 1;
  std::ostringstream err;
  std::streambuf *old = std::cerr.rdbuf(err.rdbuf());
  std::unique_ptr<tadah::LammpsTadah> lt;
  std::string msg;
  try {
    lt = load(pot.path);
  } catch (const std::exception &e) {
    msg = e.what();
  }
  std::cerr.rdbuf(old);
  INFO("load: " << msg << "\nstderr: " << err.str());
  REQUIRE(lt != nullptr);
  CHECK_FALSE(lt->c.has_stored("SIGMA"));
  std::size_t n = 0;
  const std::string text = err.str();
  for (std::size_t p = text.find("Deprecated key SIGMA is ignored");
       p != std::string::npos;
       p = text.find("Deprecated key SIGMA is ignored", p + 1))
    ++n;
  CHECK(n == expected);
}
