/*
 * Copyright 2022 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "libint_basis_set_water.hpp"
#include "make_engine.hpp"

const double eps  = 100000000 * std::numeric_limits<double>::epsilon();
const double marg = 1000000 * std::numeric_limits<double>::epsilon();

using integrals::libint::detail_::make_engine;

TEST_CASE("make_engine") {
    /// Libint Basis Set
    auto bset = testing::water_basis_set();

    /// Basis set inputs
    std::vector two_sets{bset, bset};
    std::vector three_sets{bset, bset, bset};
    std::vector four_sets{bset, bset, bset, bset};

    /// Threshold input
    double t = 1.0E-16;

    /// Components for operators
    chemist::Electron e{};
    chemist::operators::STG stg{1.0, 1.0};
    chemist::Atom o1{"O", 8ul, 0.0, 0.0, -0.143222342980786, 0.0};
    chemist::Atom h1{"H", 1ul, 0.0, 1.638033502034240, 1.136556880358410, 0.0};
    chemist::Atom h2{"H", 1ul, 0.0, -1.638033502034240, 1.136556880358410, 0.0};
    chemist::Nuclei water{o1, h1, h2};

    SECTION("el_el_coulomb") {
        using op_t = simde::type::el_el_coulomb;
        op_t op;

        SECTION("Two basis sets") {
            auto engine = make_engine(two_sets, op, t);
            testing::test_engine_standard<op_t>(engine);
            REQUIRE(engine.braket() == libint2::BraKet::xs_xs);
        }

        SECTION("Three basis sets") {
            auto engine = make_engine(three_sets, op, t);
            testing::test_engine_standard<op_t>(engine);
            REQUIRE(engine.braket() == libint2::BraKet::xs_xx);
        }

        SECTION("Four basis sets") {
            auto engine = make_engine(four_sets, op, t);
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
        auto engine = make_engine(two_sets, op, t);
        testing::test_engine_standard<op_t>(engine);
        const auto& buf = engine.compute(two_sets[0][0], two_sets[1][0]);
        REQUIRE(buf[0][0] ==
                Approx(-61.5805952694322).epsilon(eps).margin(marg));
    }

    SECTION("el_el_stg") {
        using op_t = simde::type::el_el_stg;
        op_t op(stg);
        auto engine = make_engine(two_sets, op, t);
        testing::test_engine_standard<op_t>(engine);
        const auto& buf = engine.compute(two_sets[0][0], two_sets[1][0]);
        REQUIRE(buf[0][0] ==
                Approx(0.2537942162587338).epsilon(eps).margin(marg));
    }

    SECTION("el_el_yukawa") {
        using op_t = simde::type::el_el_yukawa;
        op_t op(stg, e, e);
        auto engine = make_engine(two_sets, op, t);
        testing::test_engine_standard<op_t>(engine);
        const auto& buf = engine.compute(two_sets[0][0], two_sets[1][0]);
        REQUIRE(buf[0][0] ==
                Approx(0.7174258302062281).epsilon(eps).margin(marg));
    }

    SECTION("el_el_f12_commutator") {
        using op_t = simde::type::el_el_f12_commutator;
        op_t op(stg, e, e);
        auto engine = make_engine(two_sets, op, t);
        testing::test_engine_standard<op_t>(engine);
        const auto& buf = engine.compute(two_sets[0][0], two_sets[1][0]);
        REQUIRE(buf[0][0] == Approx(0.1612010046).epsilon(eps).margin(marg));
    }

    SECTION("el_dipole") {
        using op_t = simde::type::el_dipole;
        op_t op;
        auto engine = make_engine(two_sets, op, t);
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
        auto engine = make_engine(two_sets, op, t);
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
        auto engine = make_engine(two_sets, op, t);
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