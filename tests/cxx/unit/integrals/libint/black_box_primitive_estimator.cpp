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

/* Correct values for water comes from a companion python script that I am less
 * confident in than the H2 values. The Python script is much simpler than the
 * C++ implementation, and I didn't see any obvious errors, so it's quite
 * plausible that it's correct, but if I made a logic error in the C++ code I
 * likely made the same logic error int the Python code.
 */
template<typename FloatType>
auto corr_k01_water() {
    std::vector<FloatType> buffer{
      18.0790216813702180, 17.4852396535751247, 5.4493863186785507,
      17.4852396535751247, 16.9109596266485731, 5.2704082901341138,
      5.4493863186785507,  5.2704082901341138,  1.6425563160202101};
    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(buffer), shape);
    return simde::type::tensor(std::move(cont), shape);
}

template<typename FloatType>
auto corr_k23_water() {
    std::vector<FloatType> buffer{
      0.0000000007980106, 0.0002571668526586, 0.0041100638983562,
      0.0002571668526586, 0.0025216231166686, 0.0053704017775550,
      0.0041100638983562, 0.0053704017775550, 0.0028156052053020};
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