#include "catch2/catch.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
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

/* ---------------------------------------------------------------------- */
/* DM_REAM (EAM-form many-body block) and the EMBF0 embedding reference.   */

namespace {

/* Load a pot the way pair_coeff does; with `free` the same file is loaded
 * with an extra `EMBF0 FREE` line (written to a temporary copy). */
std::unique_ptr<LammpsTadah> load_pot(const std::string& file, bool free) {
    tadah::core::PeriodicTable::initialize();
    std::string path = file;
    if (free) {
        path = file + ".free.tmp";
        std::ifstream in(file);
        std::ofstream out(path);
        out << in.rdbuf() << "\nEMBF0 FREE\n";
    }
    const char* argv[] = {"*", "*", path.c_str(), "Ta"};
    auto lt = std::make_unique<LammpsTadah>(4, const_cast<char**>(argv));
    lt->initialize(4, const_cast<char**>(argv));
    if (free) std::remove(path.c_str());
    return lt;
}

/* w . aed over the DM block, with aed from DM_Base::calc_aed at rho_bar;
 * `abs_sum` (optional) receives sum |w . aed|, the scale of its rounding
 * (this test file is compiled with the preset's fast-math flags, so its
 * reassociated sum need not round like the library's loop). */
double embedding_via_aed(LammpsTadah& lt, const tadah::core::t_type& w,
                         double rho_bar, double* abs_sum = nullptr) {
    const std::size_t ndesc = lt.S.descriptor_size();
    tadah::core::rho_type rho;
    lt.S.dmb->init_rhoi(rho);
    rho(0) = rho_bar;
    tadah::core::aed_type aed(ndesc);
    aed.set_zero();
    lt.S.dmb->calc_aed(rho, aed);
    double F = 0.0, A = 0.0;
    const std::size_t off = ndesc - lt.S.dmb->size();
    for (std::size_t c = off; c < ndesc; ++c) {
        F += w(c) * aed(c);
        A += std::fabs(w(c) * aed(c));
    }
    if (abs_sum) *abs_sum = A;
    return F;
}

} // namespace

TEST_CASE("DM_REAM collapse: eam_embed equals w.calc_aed; EMBF0 ZERO gives F(0)=0",
          "[eam_collapse][ream][embf0]") {
    for (const std::string file : {"tests_data/pots/ta_ream_blip5.tadah",
                                   "tests_data/pots/ta_ream_knot5.tadah"}) {
        INFO(file);
        auto z_ptr = load_pot(file, false);
        auto f_ptr = load_pot(file, true);
        LammpsTadah& lz = *z_ptr;
        LammpsTadah& lf = *f_ptr;
        REQUIRE(lz.initmb);
        REQUIRE(lz.S.dmb->eam_form());
        REQUIRE(lz.S.dmb->embf0_mode() == tadah::models::EmbF0Mode::Zero);
        REQUIRE(lf.S.dmb->embf0_mode() == tadah::models::EmbF0Mode::Free);

        const auto w = weights_of(lz);
        REQUIRE(w.size() == lz.S.descriptor_size());

        EamCollapseTables T;
        std::string why;
        REQUIRE(tadah::build_eam_collapse(lz.S, w, lz.Z[0], 1024, T, why));
        REQUIRE(T.valid);
        REQUIRE(T.has_mb);

        const double raw_f0 = *lf.S.dmb->raw_embedding_f0(w);
        REQUIRE(raw_f0 != 0.0);   // the fixture is sensitive to the reference
        // Rounding floor of sums cancelling down from S0 = sum |w g(0)|
        // (every embedding family has g >= 0, so S0 = raw F(0) of |w|).
        tadah::core::t_type wabs = w;
        for (std::size_t i = 0; i < wabs.size(); ++i) wabs(i) = std::fabs(w(i));
        const double s0 = *lf.S.dmb->raw_embedding_f0(wabs);
        const double eps = std::numeric_limits<double>::epsilon();

        for (const double rb : {0.0, 1e-9, 0.05, 0.3, 0.9, 1.7, 3.1, 7.5}) {
            INFO("rho_bar = " << rb);
            double Fz, Fpz, Ff, Fpf;
            lz.S.dmb->eam_embed(w, rb, Fz, Fpz);
            lf.S.dmb->eam_embed(w, rb, Ff, Fpf);
            // F' is identical in both references (forces unchanged).
            REQUIRE(Fpz == Fpf);
            // The EAMFAST embedding equals the generic w . aed (both modes).
            double az = 0.0, af = 0.0;
            const double gz = embedding_via_aed(lz, w, rb, &az);
            const double gf = embedding_via_aed(lf, w, rb, &af);
            const double tol = 128.0 * eps * (s0 + az + af) + 1e-12;
            REQUIRE(std::fabs(Fz - gz) <= tol);
            REQUIRE(std::fabs(Ff - gf) <= tol);
            // ZERO is FREE shifted by the raw F(0).
            REQUIRE(std::fabs(Fz - (Ff - raw_f0)) <= tol);
        }
        // Exactly zero at rho = 0 under ZERO, on both evaluation paths.
        double F0, Fp0;
        lz.S.dmb->eam_embed(w, 0.0, F0, Fp0);
        REQUIRE(F0 == 0.0);
        REQUIRE(embedding_via_aed(lz, w, 0.0) == 0.0);
        // FREE: the isolated-atom embedding energy is the raw F(0).
        lf.S.dmb->eam_embed(w, 0.0, F0, Fp0);
        REQUIRE(std::fabs(F0 - raw_f0) <= 128.0 * eps * s0);
    }
}
