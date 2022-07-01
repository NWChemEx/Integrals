#include "make_engine.hpp"

const double eps  = 100000000 * std::numeric_limits<double>::epsilon();
const double marg = 1000000 * std::numeric_limits<double>::epsilon();

TEST_CASE("make_engine") {
    /// Libint Basis Set
    auto bset = libint2::BasisSet();
    bset.push_back(
      libint2::Shell{{130.7093200, 23.8088610, 6.4436083},
                     {{0, true, {0.15432897, 0.53532814, 0.44463454}}},
                     {{0.0, -0.143222342980786, 0.0}}});
    bset.push_back(
      libint2::Shell{{5.0331513, 1.1695961, 0.3803890},
                     {{0, true, {-0.09996723, 0.39951283, 0.70011547}}},
                     {{0.0, -0.143222342980786, 0.0}}});
    bset.push_back(
      libint2::Shell{{5.0331513, 1.1695961, 0.3803890},
                     {{1, true, {0.15591627, 0.60768372, 0.39195739}}},
                     {{0.0, -0.143222342980786, 0.0}}});
    bset.push_back(
      libint2::Shell{{3.425250914, 0.6239137298, 0.1688554040},
                     {{0, true, {0.1543289673, 0.5353281423, 0.4446345422}}},
                     {{1.638033502034240, 1.136556880358410, 0.0}}});
    bset.push_back(
      libint2::Shell{{3.425250914, 0.6239137298, 0.1688554040},
                     {{0, true, {0.1543289673, 0.5353281423, 0.4446345422}}},
                     {{-1.638033502034240, 1.136556880358410, 0.0}}});

    /// Basis set inputs
    std::vector<libint2::BasisSet> two_sets{bset, bset};
    std::vector<libint2::BasisSet> three_sets{bset, bset, bset};
    std::vector<libint2::BasisSet> four_sets{bset, bset, bset, bset};

    /// Threshold input
    double t = 1.0E-16;

    /// Particles for operators
    std::size_t one{1}, eight{8};
    std::array<double, 3> coords1{0.0, -0.143222342980786, 0.0};
    std::array<double, 3> coords2{1.638033502034240, 1.136556880358410, 0.0};
    std::array<double, 3> coords3{-1.638033502034240, 1.136556880358410, 0.0};

    chemist::Electron e{};
    chemist::operators::STG stg{1.0, 1.0};
    chemist::Atom o1{eight, coords1};
    chemist::Atom h1{one, coords2};
    chemist::Atom h2{one, coords3};
    chemist::Nuclei water{o1, h1, h2};

    SECTION("el_el_coulomb") {
        using op_t = simde::type::el_el_coulomb;
        op_t op;

        SECTION("Two basis sets") {
            auto engine = integrals::detail_::make_engine(two_sets, op, t);
            testing::test_engine_standard<op_t>(engine);
            REQUIRE(engine.braket() == libint2::BraKet::xs_xs);
        }

        SECTION("Three basis sets") {
            auto engine = integrals::detail_::make_engine(three_sets, op, t);
            testing::test_engine_standard<op_t>(engine);
            REQUIRE(engine.braket() == libint2::BraKet::xs_xx);
        }

        SECTION("Four basis sets") {
            auto engine = integrals::detail_::make_engine(four_sets, op, t);
            testing::test_engine_standard<op_t>(engine);
            REQUIRE(engine.braket() == libint2::BraKet::xx_xx);
        }
    }

    /// For the rest of these, the correctness of the parameter setting
    /// is tested indirectly by checking some computed values from the
    /// engine against the expected result.
    SECTION("el_nuc_coulomb") {
        using op_t = simde::type::el_nuc_coulomb;
        op_t op(e, water);
        auto engine = integrals::detail_::make_engine(two_sets, op, t);
        testing::test_engine_standard<op_t>(engine);
        const auto& buf = engine.compute(two_sets[0][0], two_sets[1][0]);
        REQUIRE(buf[0][0] ==
                Approx(-61.5805952694322).epsilon(eps).margin(marg));
    }

    SECTION("el_el_stg") {
        using op_t = simde::type::el_el_stg;
        op_t op(stg);
        auto engine = integrals::detail_::make_engine(two_sets, op, t);
        testing::test_engine_standard<op_t>(engine);
        const auto& buf = engine.compute(two_sets[0][0], two_sets[1][0]);
        REQUIRE(buf[0][0] ==
                Approx(0.2537942162587338).epsilon(eps).margin(marg));
    }

    SECTION("el_el_yukawa") {
        using op_t = simde::type::el_el_yukawa;
        op_t op(stg, e, e);
        auto engine = integrals::detail_::make_engine(two_sets, op, t);
        testing::test_engine_standard<op_t>(engine);
        const auto& buf = engine.compute(two_sets[0][0], two_sets[1][0]);
        REQUIRE(buf[0][0] ==
                Approx(0.7174258302062281).epsilon(eps).margin(marg));
    }

    SECTION("el_el_f12_commutator") {
        using op_t = simde::type::el_el_f12_commutator;
        op_t op(stg, e, e);
        auto engine = integrals::detail_::make_engine(two_sets, op, t);
        testing::test_engine_standard<op_t>(engine);
        const auto& buf = engine.compute(two_sets[0][0], two_sets[1][0]);
        REQUIRE(buf[0][0] == Approx(0.1612010046).epsilon(eps).margin(marg));
    }

    SECTION("el_dipole") {
        using op_t = simde::type::el_dipole;
        op_t op;
        auto engine = integrals::detail_::make_engine(two_sets, op, t);
        testing::test_engine_standard<op_t>(engine);
        const auto& buf = engine.results();
        engine.compute(two_sets[0][0], two_sets[1][0]);
        REQUIRE(buf[0][0] ==
                Approx(1.0000000000000004).epsilon(eps).margin(marg));
        REQUIRE(buf[2][0] ==
                Approx(-0.1432223429807861).epsilon(eps).margin(marg));
    }

    SECTION("el_quadrupole") {
        using op_t = simde::type::el_quadrupole;
        op_t op;
        auto engine = integrals::detail_::make_engine(two_sets, op, t);
        testing::test_engine_standard<op_t>(engine);
        const auto& buf = engine.results();
        engine.compute(two_sets[0][0], two_sets[1][0]);
        REQUIRE(buf[0][0] ==
                Approx(1.0000000000000004).epsilon(eps).margin(marg));
        REQUIRE(buf[2][0] ==
                Approx(-0.1432223429807861).epsilon(eps).margin(marg));
        REQUIRE(buf[4][0] ==
                Approx(0.0170208494985677).epsilon(eps).margin(marg));
    }

    SECTION("el_octupole") {
        using op_t = simde::type::el_octupole;
        op_t op;
        auto engine = integrals::detail_::make_engine(two_sets, op, t);
        testing::test_engine_standard<op_t>(engine);
        const auto& buf = engine.results();
        engine.compute(two_sets[0][0], two_sets[1][0]);
        REQUIRE(buf[0][0] ==
                Approx(1.0000000000000004).epsilon(eps).margin(marg));
        REQUIRE(buf[2][0] ==
                Approx(-0.1432223429807861).epsilon(eps).margin(marg));
        REQUIRE(buf[4][0] ==
                Approx(0.0170208494985677).epsilon(eps).margin(marg));
        REQUIRE(buf[11][0] ==
                Approx(-0.0024377659447082).epsilon(eps).margin(marg));
    }
}