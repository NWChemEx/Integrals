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

#include "detail_/libint_basis_set_water.hpp"
#include "integrals/libint/libint_factory.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>
#include <simde/types.hpp>

using namespace mokup;

const double eps  = 100000000 * std::numeric_limits<double>::epsilon();
const double marg = 1000000 * std::numeric_limits<double>::epsilon();

using op1_type = simde::type::el_identity;
using op2_type = simde::type::el_el_stg;
using test1_t  = integrals::libint::LibintFactory<2, op1_type>;
using test2_t  = integrals::libint::LibintFactory<2, op2_type>;

TEST_CASE("LibintFactory") {
    auto basis = testing::water_basis_set();

    /// Inputs
    std::vector sets1{basis, basis};
    std::vector sets2{basis, basis, basis};
    op1_type I;
    op2_type stg1(chemist::operators::STG(1.0, 1.0));
    op2_type stg2(chemist::operators::STG(2.0, 2.0));
    double thresh = 1.0E-16;
    int deriv     = 0;

    /// Construct
    test1_t fac1(sets1, I, thresh, deriv);
    test2_t fac2(sets1, stg1, thresh, deriv);
    test2_t fac3(sets1, stg2, thresh, deriv);

    SECTION("compute") {
        const auto& results = fac1.compute({2, 2});
        for(auto i = 0; i < 3; ++i) {
            for(auto j = 0; j < 3; ++j) {
                auto val = results[0][i * 3 + j];
                if(i == j) {
                    REQUIRE(val == Approx(1.0).epsilon(eps).margin(marg));
                } else {
                    REQUIRE(val == Approx(0.0).epsilon(eps).margin(marg));
                }
            }
        }
    }

    SECTION("clone") {
        auto clone = fac1.clone();
        REQUIRE(clone->are_equal(fac1));
    }

    SECTION("are_equal") {
        SECTION("Same") {
            REQUIRE(fac1.are_equal(test1_t(sets1, I, thresh, deriv)));
        }
        SECTION("Different OpType") { REQUIRE_FALSE(fac1.are_equal(fac2)); }
        SECTION("Differnt Basis sets") {
            REQUIRE_FALSE(fac1.are_equal(test1_t(sets2, I, thresh, deriv)));
        }
        SECTION("Differnt Operator") { REQUIRE_FALSE(fac2.are_equal(fac3)); }
        SECTION("Differnt Threshold") {
            REQUIRE_FALSE(fac1.are_equal(test1_t(sets1, I, 1.0E-12, deriv)));
        }
        /// Libint installation doesn't have derivatives, so this segfaults ATM.
        // SECTION("Differnt Derivatives") {
        //     REQUIRE_FALSE(fac1.are_equal(test1_t(sets1, I, thresh, 1)));
        // }
    }
}