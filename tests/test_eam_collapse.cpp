#include "catch2/catch.hpp"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <tadah/core/periodic_table.h>
#include <tadah/tadah.lammps/eam_collapse.h>
#include <tadah/tadah.lammps/lammps_tadah.h>

using tadah::EamCollapseTables;
using tadah::LammpsTadah;

namespace {

/* Build the fixture potential exactly the way pair_coeff does.
 * (unique_ptr: LammpsTadah is neither copyable nor movable.) */
std::unique_ptr<LammpsTadah> make_kr_pure2b() {
    tadah::core::PeriodicTable::initialize();
    const char* argv[] = {"*", "*", "tests_data/pots/kr_pure2b.tadah", "Kr"};
    auto lt = std::make_unique<LammpsTadah>(4, const_cast<char**>(argv));
    lt->initialize(4, const_cast<char**>(argv));
    return lt;
}

tadah::core::t_type weights_of(const LammpsTadah& lt) {
    tadah::core::t_type w;
    w.resize(lt.c.size("WEIGHTS"));
    lt.c.get<tadah::core::t_type>("WEIGHTS", w);
    return w;
}

} // namespace

TEST_CASE("build_pair2b_collapse: node-exact phi for pure-2B D2_Knot5",
          "[eam_collapse][pair2b]") {
    auto lt_ptr = make_kr_pure2b();
    LammpsTadah& lt = *lt_ptr;
    REQUIRE(lt.init2b);
    REQUIRE_FALSE(lt.initmb);

    const auto w = weights_of(lt);
    const int Z = lt.Z[0];
    constexpr int NR = 2048;

    EamCollapseTables T;
    std::string why;
    REQUIRE(tadah::build_pair2b_collapse(lt.S, w, Z, NR, T, why));
    REQUIRE(T.valid);
    REQUIRE_FALSE(T.has_mb);
    REQUIRE(T.phi.nr == NR);
    REQUIRE(T.phi.rmax == Approx(lt.S.d2b->get_rcut()));
    REQUIRE(T.rho.nr == 0);           // no density table in the 2b collapse

    /* Node exactness: table nodes must equal w . calc_aed / w . calc_dXijdri
     * bit-for-bit up to summation order (same loop as the builder). */
    const std::size_t ndesc = lt.S.descriptor_size();
    REQUIRE(w.size() == ndesc);
    double max_v = 0.0;
    for (int k = 1; k < NR; k += NR / 37) {
        const double r = k * T.phi.dr;
        tadah::core::aed_type aed(ndesc);
        tadah::core::fd_type  fd(ndesc);
        aed.set_zero();
        lt.S.d2b->calc_aed(Z, Z, r, r * r, aed, 1.0);
        fd.set_zero();
        lt.S.d2b->calc_dXijdri(Z, Z, r, r * r, fd, 1.0);
        double pv = 0.0, pd = 0.0;
        for (std::size_t c = 0; c < ndesc; ++c) {
            pv += w(c) * aed(c);
            pd += w(c) * fd(c, 0);
        }
        REQUIRE(T.phi.v[k] == Approx(pv).margin(1e-14));
        REQUIRE(T.phi.d[k] == Approx(pd).margin(1e-14));
        max_v = std::max(max_v, std::fabs(pv));
    }

    /* Mid-node parity self-check recorded by the builder. */
    REQUIRE(T.parity_err_phi < 1e-8 * std::max(1.0, max_v));

    /* Hermite self-consistency: derivative returned by eval() must match a
     * central difference of the interpolated value (same segment). */
    const double r0 = 0.43 * T.phi.rmax;
    const double h = 1e-6;
    double v_lo, v_hi, v_mid, d_mid, dummy;
    T.phi.eval(r0 - h, v_lo, dummy);
    T.phi.eval(r0 + h, v_hi, dummy);
    T.phi.eval(r0, v_mid, d_mid);
    REQUIRE(d_mid == Approx((v_hi - v_lo) / (2.0 * h)).epsilon(1e-6).margin(1e-8));
}

TEST_CASE("build_eam_collapse refuses a pure-2B model (DM_Dummy is not EAM-form)",
          "[eam_collapse][pair2b]") {
    auto lt_ptr = make_kr_pure2b();
    LammpsTadah& lt = *lt_ptr;
    const auto w = weights_of(lt);

    EamCollapseTables T;
    std::string why;
    REQUIRE_FALSE(tadah::build_eam_collapse(lt.S, w, lt.Z[0], 512, T, why));
    REQUIRE_FALSE(T.valid);
    REQUIRE(why.find("not EAM-form") != std::string::npos);
}

TEST_CASE("build_pair2b_collapse refuses on weight-size mismatch",
          "[eam_collapse][pair2b]") {
    auto lt_ptr = make_kr_pure2b();
    LammpsTadah& lt = *lt_ptr;

    tadah::core::t_type w_bad;
    w_bad.resize(lt.S.descriptor_size() + 3);

    EamCollapseTables T;
    std::string why;
    REQUIRE_FALSE(tadah::build_pair2b_collapse(lt.S, w_bad, lt.Z[0], 512, T, why));
    REQUIRE_FALSE(T.valid);
    REQUIRE(why.find("weight count") != std::string::npos);
}
