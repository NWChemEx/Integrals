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

} // namespace

TEST_CASE("BlackBoxPrimitiveEstimator") {
    auto mm   = initialize_integrals();
    auto& mod = mm.at("Black Box Primitive Pair Estimator");
    auto&& [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();

    using float_type = double;

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