#include "catch2/catch.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <unistd.h>

#include <tadah/core/core_user_config_parser.h>
#include <tadah/core/periodic_table.h>
#include <tadah/models/descriptors/dm/dm_ream.h>
#include <tadah/models/descriptors/dm/rhoform.h>
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

TEST_CASE("DM_REAM collapse: eam_embed equals w.calc_aed, EMBF0 ZERO gives F(0)=0",
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

/* ---------------------------------------------------------------------- */
/* DM_REAM density forms (RHOFORM) in the EAM collapse.                    */

namespace {

/* Copy of @p file with its RHOFORM line replaced by "RHOFORM <form>" (form
 * empty: the line is dropped), written next to it under @p tag and the
 * process id (ctest runs the cases as concurrent processes). */
std::string pot_with_rhoform(const std::string& file, const std::string& form,
                             const std::string& tag) {
    std::ifstream in(file);
    REQUIRE(in.good());
    const std::string path = file + "." + tag + "." +
                             std::to_string(::getpid()) + ".tmp";
    std::ofstream out(path);
    std::string line;
    int replaced = 0;
    while (std::getline(in, line)) {
        if (line.rfind("RHOFORM", 0) == 0) {
            ++replaced;
            if (!form.empty()) out << "RHOFORM           " << form << "\n";
            continue;
        }
        out << line << "\n";
    }
    REQUIRE(replaced == 1);
    return path;
}

std::unique_ptr<LammpsTadah> load_path(const std::string& path) {
    tadah::core::PeriodicTable::initialize();
    const char* argv[] = {"*", "*", path.c_str(), "Ta"};
    auto lt = std::make_unique<LammpsTadah>(4, const_cast<char**>(argv));
    lt->initialize(4, const_cast<char**>(argv));
    return lt;
}

std::unique_ptr<LammpsTadah> load_with_rhoform(const std::string& file,
                                               const std::string& form) {
    const std::string path = pot_with_rhoform(file, form, form);
    auto lt = load_path(path);
    std::remove(path.c_str());
    return lt;
}

tadah::models::RhoForm density_form(const LammpsTadah& lt) {
    const auto* ream =
        dynamic_cast<const tadah::models::DM_REAM*>(lt.S.dmb);
    REQUIRE(ream != nullptr);
    return ream->density().form;
}

} // namespace

TEST_CASE("DM_REAM collapse under RHOFORM MONOTONE: the rho table follows "
          "eam_density (relative tolerance)",
          "[eam_collapse][ream][rhoform]") {
    /* The Blip5 fixture has every density amplitude > 0, SGRIDMB > 0,
     * WATOM 1 and Cut_Dummy, so it is also a valid MONOTONE potential:
     * rho = sum_n a_n S_n(r), S_n = 1 inside the bump's inner support edge
     * (c_n - 3/s_n), 0 beyond its outer edge (c_n + 3/s_n). This file is
     * compiled with the preset's fast-math flags, so every comparison with
     * the library (built without them) is relative, never bitwise. */
    const std::string file = "tests_data/pots/ta_ream_blip5.tadah";
    auto m_ptr = load_with_rhoform(file, "MONOTONE");
    auto s_ptr = load_path(file);   // RHOFORM SIGNED, as written
    LammpsTadah& lm = *m_ptr;
    LammpsTadah& ls = *s_ptr;
    REQUIRE(tadah::models::stored_rhoform(lm.c) ==
            tadah::models::RhoForm::Monotone);
    REQUIRE(density_form(lm) == tadah::models::RhoForm::Monotone);
    REQUIRE(density_form(ls) == tadah::models::RhoForm::Signed);

    const auto w = weights_of(lm);
    constexpr int NR = 8192;   // the table size pair_tadah uses
    EamCollapseTables T, TS;
    std::string why;
    REQUIRE(tadah::build_eam_collapse(lm.S, w, lm.Z[0], NR, T, why));
    REQUIRE(tadah::build_eam_collapse(ls.S, w, ls.Z[0], NR, TS, why));
    REQUIRE(T.valid);
    REQUIRE(T.has_mb);
    const int Z = lm.Z[0];

    double rho_max = 0.0, rhop_max = 0.0;
    for (int k = 1; k < NR; ++k) {
        rho_max = std::max(rho_max, std::fabs(T.rho.v[k]));
        rhop_max = std::max(rhop_max, std::fabs(T.rho.d[k]));
    }
    REQUIRE(rho_max > 0.0);
    REQUIRE(rhop_max > 0.0);
    const double eps = std::numeric_limits<double>::epsilon();

    SECTION("nodes are eam_density; rho >= 0 and non-increasing") {
        for (int k = 1; k < NR; ++k) {
            const double r = k * T.rho.dr;
            double v, d;
            lm.S.dmb->eam_density(Z, Z, r, v, d);
            REQUIRE(std::fabs(T.rho.v[k] - v) <= 4.0 * eps * rho_max);
            REQUIRE(std::fabs(T.rho.d[k] - d) <= 4.0 * eps * rhop_max);
            // library values: the form's guarantees, to rounding
            REQUIRE(T.rho.v[k] >= 0.0);
            REQUIRE(T.rho.d[k] <= 16.0 * eps * rhop_max);
            if (k + 1 < NR)
                REQUIRE(T.rho.v[k + 1] - T.rho.v[k] <= 16.0 * eps * rho_max);
        }
    }

    SECTION("off-node interpolation: value and derivative, relative") {
        REQUIRE(T.parity_err_rho <= 1e-10 * rho_max);
        // 4001 points strictly inside (0, rcut), none on a node
        const double rmax = T.rho.rmax;
        double worst_v = 0.0, worst_d = 0.0;
        for (int i = 1; i <= 4001; ++i) {
            const double r = rmax * (i - 0.37) / 4002.0;
            double tv, td, av, ad;
            T.rho.eval(r, tv, td);
            lm.S.dmb->eam_density(Z, Z, r, av, ad);
            worst_v = std::max(worst_v, std::fabs(tv - av));
            worst_d = std::max(worst_d, std::fabs(td - ad));
        }
        INFO("worst |table - eam_density| = " << worst_v << " (rho_max "
             << rho_max << "), worst |d| = " << worst_d << " (rhop_max "
             << rhop_max << ")");
        REQUIRE(worst_v <= 1e-10 * rho_max);
        REQUIRE(worst_d <= 1e-6 * rhop_max);
    }

    SECTION("the step form, not the bump form, is tabulated") {
        // Inside every inner support edge (1.8 A): S_n = 1, rho = sum a_n.
        // Beyond every outer edge (5.4 A) and below RCUTMB: rho = 0; the
        // Hermite segments there interpolate zeros (undershoot ~1e-22).
        const double sum_a = 0.5 + 0.3 + 0.1;
        for (const double r : {0.3, 0.9, 1.5, 1.75}) {
            INFO("r = " << r);
            REQUIRE(T.rho.eval_val(r) == Approx(sum_a).epsilon(1e-12));
            // SIGNED: every bump is 0 there
            REQUIRE(std::fabs(TS.rho.eval_val(r)) <= 1e-12 * rho_max);
        }
        for (const double r : {5.41, 5.45, 5.49}) {
            INFO("r = " << r);
            REQUIRE(std::fabs(T.rho.eval_val(r)) <= 1e-12 * rho_max);
        }
        // The pair function and the embedding do not depend on the form.
        for (int k = 1; k < NR; k += 97) {
            REQUIRE(T.phi.v[k] == Approx(TS.phi.v[k]).epsilon(4.0 * eps)
                                      .margin(4.0 * eps));
        }
        for (const double rb : {0.0, 0.4, 0.9, 2.5}) {
            double Fm, Fpm, Fs, Fps;
            lm.S.dmb->eam_embed(w, rb, Fm, Fpm);
            ls.S.dmb->eam_embed(w, rb, Fs, Fps);
            REQUIRE(Fm == Approx(Fs).epsilon(4.0 * eps).margin(4.0 * eps));
            REQUIRE(Fpm == Approx(Fps).epsilon(4.0 * eps).margin(4.0 * eps));
        }
    }

    SECTION("NONNEG tabulates the SIGNED density of the same amplitudes") {
        // Every amplitude >= 0 and Cut_Dummy >= 0: the NONNEG density is the
        // SIGNED one (bit for bit in the library; relative here).
        auto n_ptr = load_with_rhoform(file, "NONNEG");
        REQUIRE(density_form(*n_ptr) == tadah::models::RhoForm::NonNeg);
        EamCollapseTables TN;
        REQUIRE(tadah::build_eam_collapse(n_ptr->S, w, Z, NR, TN, why));
        double smax = 0.0;
        for (int k = 1; k < NR; ++k)
            smax = std::max(smax, std::fabs(TS.rho.v[k]));
        for (int k = 1; k < NR; ++k)
            REQUIRE(std::fabs(TN.rho.v[k] - TS.rho.v[k]) <= 4.0 * eps * smax);
    }
}

TEST_CASE("DM_REAM potential without RHOFORM: the load fails with the "
          "migration text",
          "[ream][rhoform]") {
    /* PairTadah::load keeps e.what() of exactly these exceptions and
     * pair_coeff puts it in its error->all message (checked through LAMMPS
     * in the HPO test_ream_fastpath). */
    tadah::core::PeriodicTable::initialize();
    const std::string file = "tests_data/pots/ta_ream_knot5.tadah";

    SECTION("pair_coeff * * <file> (the CLI constructor)") {
        const std::string path = pot_with_rhoform(file, "", "absent");
        const char* argv[] = {"*", "*", path.c_str(), "Ta"};
        REQUIRE_THROWS_AS(LammpsTadah(4, const_cast<char**>(argv)),
                          tadah::models::DensityFormError);
        REQUIRE_THROWS_WITH(LammpsTadah(4, const_cast<char**>(argv)),
                            Catch::Contains(
                                tadah::models::rhoform_absent_message()) &&
                            Catch::Contains("add 'RHOFORM SIGNED' to use an "
                                            "existing potential unchanged"));
        std::remove(path.c_str());
    }
    SECTION("a driver-built potential (the Context constructor)") {
        const std::string path = pot_with_rhoform(file, "", "absent_ctx");
        tadah::core::Context ctx =
            tadah::core::CoreUserConfigParser().parse(path);
        std::remove(path.c_str());
        REQUIRE_FALSE(ctx.has_stored("RHOFORM"));
        REQUIRE_THROWS_WITH(LammpsTadah(ctx),
                            Catch::Contains("DM_REAM: RHOFORM is required"));
        REQUIRE_FALSE(ctx.has_stored("RHOFORM"));   // nothing materialised
    }
    SECTION("an invalid form") {
        const std::string path = pot_with_rhoform(file, "SQUARE", "square");
        const char* argv[] = {"*", "*", path.c_str(), "Ta"};
        REQUIRE_THROWS_WITH(LammpsTadah(4, const_cast<char**>(argv)),
                            Catch::Contains("RHOFORM 'SQUARE' is not a "
                                            "density form"));
        std::remove(path.c_str());
    }
    SECTION("the fixture as written (RHOFORM SIGNED) loads") {
        auto lt = load_path(file);
        REQUIRE(density_form(*lt) == tadah::models::RhoForm::Signed);
    }
}
