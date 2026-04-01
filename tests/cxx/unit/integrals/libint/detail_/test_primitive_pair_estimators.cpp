/*
 * Copyright 2026 NWChemEx-Project
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

#include "integrals/libint/detail_/primitive_pair_estimators.hpp"
#undef DEPRECATED
#include "../../testing/testing.hpp"

// Expected values below were recorded from a single run of each estimator on
// water_sto3g_basis_set() (same basis for both arguments). If the test basis or
// libint embedding changes, re-record by temporarily printing the outputs.

using namespace integrals::testing;
using integrals::libint::detail_::coarse_k_ij;
using integrals::libint::detail_::fine_k_ij;
using integrals::libint::detail_::gamma_ij;
using integrals::libint::detail_::k_ij;

namespace {

constexpr double tol = 1.0e-11;

} // namespace

TEST_CASE("primitive_pair_estimators: gamma_ij") {
    const auto aobs = water_sto3g_basis_set();
    REQUIRE(aobs.n_primitives() == 15);
    const auto g = gamma_ij(aobs, aobs);
    REQUIRE(g.size() == 15);
    REQUIRE(g[0].size() == 15);
    REQUIRE(g[0][0] == Catch::Approx(2.61418639999999982e+02).margin(tol));
    REQUIRE(g[0][14] == Catch::Approx(1.30878175403999990e+02).margin(tol));
    REQUIRE(g[14][0] == Catch::Approx(1.30878175403999990e+02).margin(tol));
    REQUIRE(g[14][14] == Catch::Approx(3.37710807999999973e-01).margin(tol));
    REQUIRE(g[1][2] == Catch::Approx(3.02524693000000013e+01).margin(tol));
}

TEST_CASE("primitive_pair_estimators: k_ij") {
    const auto aobs = water_sto3g_basis_set();
    REQUIRE(aobs.n_primitives() == 15);
    const auto K = k_ij(aobs, aobs);
    REQUIRE(K.size() == 15);
    REQUIRE(K[0].size() == 15);
    for(std::size_t i = 0; i < 15; ++i) {
        REQUIRE(K[i][i] == Catch::Approx(1.0).margin(1.0e-14));
    }
    REQUIRE(K[0][9] == Catch::Approx(5.44974713611047489e-07).margin(1.0e-18));
}

TEST_CASE("primitive_pair_estimators: coarse_k_ij") {
    const auto aobs = water_sto3g_basis_set();
    REQUIRE(aobs.n_primitives() == 15);
    const auto Kc = coarse_k_ij(aobs, aobs);
    REQUIRE(Kc.size() == 15);
    REQUIRE(Kc[0].size() == 15);
    REQUIRE(Kc[0][0] == Catch::Approx(1.80790216813702180e+01).margin(tol));
    REQUIRE(Kc[0][14] == Catch::Approx(1.71267463425960637e-01).margin(tol));
    REQUIRE(Kc[14][0] == Catch::Approx(1.71267463425960637e-01).margin(tol));
    REQUIRE(Kc[14][14] == Catch::Approx(6.96785377189288076e-03).margin(tol));
    REQUIRE(Kc[1][2] == Catch::Approx(5.27040829013411383e+00).margin(tol));
}

TEST_CASE("primitive_pair_estimators: fine_k_ij") {
    const auto aobs = water_sto3g_basis_set();
    REQUIRE(aobs.n_primitives() == 15);
    const auto Kf = fine_k_ij(aobs, aobs);
    REQUIRE(Kf.size() == 15);
    REQUIRE(Kf[0].size() == 15);
    REQUIRE(Kf[0][0] == Catch::Approx(4.09063484384912246e-01).margin(tol));
    REQUIRE(Kf[0][14] == Catch::Approx(7.74033883652055620e-03).margin(tol));
    REQUIRE(Kf[14][0] == Catch::Approx(7.74033883652055620e-03).margin(tol));
    REQUIRE(Kf[14][14] == Catch::Approx(1.22041182423709954e-01).margin(tol));
    REQUIRE(Kf[1][2] == Catch::Approx(1.03047099111916607e+00).margin(tol));
}
