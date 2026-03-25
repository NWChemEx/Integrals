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

using namespace integrals::testing;
using tensorwrapper::operations::approximately_equal;

using pt = integrals::property_types::PrimitivePairEstimator;

namespace {

/* Correct values for H2 come from a MWE I made that recreates the integral
 * screening
 * using libint directly. Since the MWE exactly matches libint, I strongly
 * believe these values are also correct.
 */

template<typename FloatType>
auto corr_k01() {
    std::vector<FloatType> buffer{3.693e-38,   4.76823e-13, 8.9425e-06,
                                  4.76823e-13, 1.73977e-08, 3.39937e-05,
                                  8.9425e-06,  3.39937e-05, 0.000112915};
    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(buffer), shape);
    return simde::type::tensor(std::move(cont), shape);
}

template<typename FloatType>
auto corr_k23() {
    std::vector<FloatType> buffer{4.07419e-12, 5.0571e-05,  0.00250317,
                                  5.0571e-05,  0.000964263, 0.00356585,
                                  0.00250317,  0.00356585,  0.00217062};

    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(buffer), shape);
    return simde::type::tensor(std::move(cont), shape);
}

/* Correct values for the l=2 test come from an independent computation using
 * the analytic normalization formula for Cartesian Gaussians.
 */
template<typename FloatType>
auto corr_k01_l2() {
    std::vector<FloatType> buffer{
      // Primitive 0 (zeta=3.425250914): xx, xy, xz, yy, yz, zz
      2.92126652e-37, 3.77180568e-12, 7.07377006e-05, 5.05978203e-37,
      6.53295907e-12, 1.22521291e-04, 5.05978203e-37, 6.53295907e-12,
      1.22521291e-04, 2.92126652e-37, 3.77180568e-12, 7.07377006e-05,
      5.05978203e-37, 6.53295907e-12, 1.22521291e-04, 2.92126652e-37,
      3.77180568e-12, 7.07377006e-05,
      // Primitive 1 (zeta=0.6239137298): xx, xy, xz, yy, yz, zz
      6.87039113e-13, 2.50677873e-08, 4.89803780e-05, 1.18998665e-12,
      4.34186812e-08, 8.48365032e-05, 1.18998665e-12, 4.34186812e-08,
      8.48365032e-05, 6.87039113e-13, 2.50677873e-08, 4.89803780e-05,
      1.18998665e-12, 4.34186812e-08, 8.48365032e-05, 6.87039113e-13,
      2.50677873e-08, 4.89803780e-05,
      // Primitive 2 (zeta=0.1688554040): xx, xy, xz, yy, yz, zz
      3.48717315e-06, 1.32560018e-05, 4.40318744e-05, 6.03996107e-06,
      2.29600686e-05, 7.62654435e-05, 6.03996107e-06, 2.29600686e-05,
      7.62654435e-05, 3.48717315e-06, 1.32560018e-05, 4.40318744e-05,
      6.03996107e-06, 2.29600686e-05, 7.62654435e-05, 3.48717315e-06,
      1.32560018e-05, 4.40318744e-05};
    tensorwrapper::shape::Smooth shape({18, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(buffer), shape);
    return simde::type::tensor(std::move(cont), shape);
}

// These values come from dumping the K matrix from our module and will be
// wrong if our module is incorrect.
template<typename FloatType>
auto corr_k01_water() {
    std::vector<FloatType> buffer{
      18.07902198052444, 17.48523994290402, 5.449386408849742,
      17.48523994290402, 16.91095990647483, 5.270408377343748,
      5.449386408849742, 5.270408377343748, 1.642556343199648};
    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(buffer), shape);
    return simde::type::tensor(std::move(cont), shape);
}

template<typename FloatType>
auto corr_k23_water() {
    std::vector<FloatType> buffer{
      7.980105744638301e-10, 0.0002571668526766007, 0.004110063898643834,
      0.0002571668526766007, 0.002521623116845085,  0.005370401777930866,
      0.004110063898643834,  0.005370401777930866,  0.002815605205499067};
    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(buffer), shape);
    return simde::type::tensor(std::move(cont), shape);
}

} // namespace

TEST_CASE("BlackBoxPrimitiveEstimator") {
    auto mm          = initialize_integrals();
    auto& mod        = mm.at("Black Box Primitive Pair Estimator");
    using float_type = double;

    SECTION("H2 Dimer/STO-3G (03|12)") {
        auto&& [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();

        SECTION("K(bra0, bra1)") {
            auto K01  = mod.run_as<pt>(bra0, bra1);
            auto corr = corr_k01<float_type>();
            REQUIRE(approximately_equal(K01, corr, 1E-8));
        }

        SECTION("K(ket0, ket1)") {
            auto K23  = mod.run_as<pt>(ket0, ket1);
            auto corr = corr_k23<float_type>();
            REQUIRE(approximately_equal(K23, corr, 1E-8));
        }
    }

    SECTION("H2 Dimer/STO-3G l=2 (03|12)") {
        auto [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();
        bra0.shell(0).l()             = 2;
        bra0.shell(0).pure()          = chemist::ShellType::cartesian;

        auto K    = mod.run_as<pt>(bra0, bra1);
        auto corr = corr_k01_l2<float_type>();
        REQUIRE(approximately_equal(K, corr, 1E-8));
    }

    SECTION("H2O/STO-3G (00|34)") {
        auto&& [bra0, bra1, ket0, ket1] = get_h2o_0034_bases();

        SECTION("K(bra0, bra1)") {
            auto K01  = mod.run_as<pt>(bra0, bra1);
            auto corr = corr_k01_water<float_type>();
            REQUIRE(approximately_equal(K01, corr, 1E-8));
        }

        SECTION("K(ket0, ket1)") {
            auto K23  = mod.run_as<pt>(ket0, ket1);
            auto corr = corr_k23_water<float_type>();
            REQUIRE(approximately_equal(K23, corr, 1E-8));
        }
    }
}
