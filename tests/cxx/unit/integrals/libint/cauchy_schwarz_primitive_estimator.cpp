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
#include "../testing/testing.hpp"
#include <integrals/property_types.hpp>

using pt = integrals::property_types::PrimitivePairEstimator;
using namespace integrals::testing;
using tensorwrapper::operations::approximately_equal;

namespace {

/* Corr values come from a companion Python script that uses the dumped ERIs.
 * The script uses Numpy and is much cleaner than the C++ code, which decreases
 * the likelihood of errors in the corr values. That said, it is possible that
 * there is an error in how I approached this. If such an error occurred it
 * would likely show up in both places and NOT be caught by these tests.
 */

// Assumes all s-type functions
template<typename FloatType>
auto Q_ab_corr() {
    std::vector<FloatType> data{0.00000000e+00, 0.00000000e+00, 3.96602449e-10,
                                0.00000000e+00, 3.75905705e-15, 2.01602647e-08,
                                3.96602449e-10, 2.01602647e-08, 8.48436812e-07};
    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(data), shape);
    return simde::type::tensor(std::move(cont), shape);
}

// Assumes p for basis 0 and s for basis 1 functions
template<typename FloatType>
auto Q_ab_psps_corr() {
    std::vector<FloatType> data{0.00000000e+00, 0.00000000e+00, 2.88455917e-09,
                                0.00000000e+00, 1.93317006e-13, 1.96314576e-07,
                                1.03768892e-08, 3.60782670e-07, 6.23308753e-06};
    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(data), shape);
    return simde::type::tensor(std::move(cont), shape);
}

// Assumes all p-type functions
template<typename FloatType>
auto Q_ab_pppp_corr() {
    std::vector<FloatType> data{0.00000000e+00, 0.00000000e+00, 6.85302613e-08,
                                0.00000000e+00, 9.31619988e-12, 3.04372388e-06,
                                6.85302613e-08, 3.04372388e-06, 3.65812603e-05};
    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(data), shape);
    return simde::type::tensor(std::move(cont), shape);
}

template<typename FloatType>
auto Q_cd_corr() {
    std::vector<FloatType> data{0.00000000,     2.08398137e-08, 3.10753955e-05,
                                2.08398137e-08, 1.15521140e-05, 2.21832728e-04,
                                3.10753955e-05, 2.21832728e-04, 3.13532082e-04};
    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(data), shape);
    return simde::type::tensor(std::move(cont), shape);
}

} // namespace

TEST_CASE("CauchySchwarzPrimitiveEstimator") {
    auto mm                       = initialize_integrals();
    auto& mod                     = mm.at("CauchySchwarz Estimator");
    auto [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();

    SECTION("Q(bra0, bra1)") {
        auto Q_ab = mod.run_as<pt>(bra0, bra1);
        auto corr = Q_ab_corr<double>();
        REQUIRE(approximately_equal(Q_ab, corr, 1E-8));
    }

    SECTION("Q(ket0, ket1)") {
        auto Q_cd = mod.run_as<pt>(ket0, ket1);
        auto corr = Q_cd_corr<double>();
        REQUIRE(approximately_equal(Q_cd, corr, 1E-8));
    }

    SECTION("Q_ab (ps|ps)") {
        bra0.shell(0).l() = 1;
        auto Q_ab         = mod.run_as<pt>(bra0, bra1);
        auto corr         = Q_ab_psps_corr<double>();
        REQUIRE(approximately_equal(Q_ab, corr, 1E-8));
    }

    SECTION("Q_ab (pp|pp)") {
        bra0.shell(0).l() = 1;
        bra1.shell(0).l() = 1;
        auto Q_ab         = mod.run_as<pt>(bra0, bra1);
        auto corr         = Q_ab_pppp_corr<double>();
        REQUIRE(approximately_equal(Q_ab, corr, 1E-8));
    }
}
