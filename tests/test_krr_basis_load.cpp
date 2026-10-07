// pair_style tadah and the basis of a non-linear KRR potential (M_KRR with
// Kern_RBF, Kern_Polynomial, ...): every K(b, x) compares a descriptor x
// with a basis vector b of the potential's BASIS (SBASIS vectors stored
// one after another), so the vectors must be there, one per weight, and
// have the length of the descriptor. Only the pair style knows that
// length (DC_Selector::descriptor_size()); M_KRR_Core checks the rest.
//
// Before the fix: a potential with basis vectors of another length
// loaded, and every K(b, x) read past the shorter of the two; one without
// BASIS loaded with an empty basis and predicted E = 0 for every atom.

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

/* An M_KRR Kern_RBF potential on D2_Blip with 3 functions (descriptor
 * length 3, BIAS false), SBASIS 2, with the given BASIS and WEIGHTS lines
 * ("" leaves the line out). */
std::vector<std::string> rbf_potential(const std::string& basis,
                                       const std::string& weights) {
    std::vector<std::string> lines = {
        "ATOM Ta",
        "WATOM 1",
        "INIT2B true",
        "INITMB false",
        "TYPE2B D2_Blip 3 3 Ta Ta",
        "CGRID2B 2.8 3.5 4.2",
        "SGRID2B 2.0 2.0 2.0",
        "RCTYPE2B Cut_Cos",
        "RCUT2B 5.0",
        "BIAS false",
        "NORM false",
        "VERBOSE 0",
        "MODEL M_KRR Kern_RBF 1",
        "MPARAMS 1.0",
        "SBASIS 2",
    };
    if (!basis.empty()) lines.push_back("BASIS " + basis);
    if (!weights.empty()) lines.push_back("WEIGHTS " + weights);
    return lines;
}

/* The message of the exception loading `lines` throws, "" when none. */
std::string load_error(const std::vector<std::string>& lines) {
    PeriodicTable::initialize();
    try {
        LammpsTadah lt(CoreUserConfigParser().parse(lines));
    } catch (const std::exception& e) {
        return e.what();
    }
    return "";
}

bool contains(const std::string& s, const std::string& part) {
    return s.find(part) != std::string::npos;
}

}  // namespace

/* Pins: the potential with 2 basis vectors of 3 values and 2 weights
 * loads (a pin); with 2 vectors of 4 values (BASIS of 8 values, a multiple
 * of SBASIS 2) it is refused, naming BASIS, the vector length and the
 * descriptor's ("SBASIS 2 basis vectors of 4 values each", "has 3
 * components").
 * Before the fix: the potential with 4-value vectors loaded (no throw). */
TEST_CASE("pair_style tadah refuses KRR basis vectors that do not have the "
          "descriptor's length",
          "[lammps_tadah][krr_basis]") {
    CHECK(load_error(rbf_potential("1 0 0 0.75 1.06 0.49", "0.5 -0.25")) == "");
    const std::string msg =
        load_error(rbf_potential("1 0 0 0 0.75 1.06 0.49 0", "0.5 -0.25"));
    INFO("message: " << msg);
    CHECK(contains(msg, "BASIS holds SBASIS 2 basis vectors of 4 values each"));
    CHECK(contains(msg, "has 3 components"));
}

/* Pins: a KRR potential of a non-linear kernel with WEIGHTS and no BASIS
 * is refused, naming SBASIS and BASIS; one with 3 WEIGHTS for 2 basis
 * vectors is refused, naming WEIGHTS and SBASIS.
 * Before the fix: both loaded; the first predicted E = 0 for every atom
 * (no basis vector, nothing to sum), the second read past its basis. */
TEST_CASE("pair_style tadah refuses a KRR potential whose basis does not "
          "match its WEIGHTS",
          "[lammps_tadah][krr_basis]") {
    {
        const std::string msg = load_error(rbf_potential("", "0.5 -0.25"));
        INFO("message: " << msg);
        CHECK(contains(msg, "no basis vectors (SBASIS and BASIS)"));
    }
    {
        const std::string msg =
            load_error(rbf_potential("1 0 0 0.75 1.06 0.49", "0.5 -0.25 0.1"));
        INFO("message: " << msg);
        CHECK(contains(msg, "3 WEIGHTS but SBASIS 2"));
    }
}
