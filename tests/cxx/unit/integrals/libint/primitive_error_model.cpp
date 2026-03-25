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
#include <integrals/integrals.hpp>

using eri4_pt = simde::ERI4;
using pt      = integrals::property_types::Uncertainty<eri4_pt>;

using namespace integrals::testing;

namespace {

template<typename T0, typename T1>
bool compare_magnitudes(T0&& lhs, T1&& corr, double tol) {
    using float_type = double;
    using tensorwrapper::buffer::get_raw_data;
    const auto lhs_buffer  = get_raw_data<float_type>(lhs.buffer());
    const auto corr_buffer = get_raw_data<float_type>(corr.buffer());
    assert(lhs_buffer.size() == corr_buffer.size());

    for(std::size_t i = 0; i < lhs_buffer.size(); ++i) {
        auto lhs_val  = lhs_buffer[i];
        auto corr_val = corr_buffer[i];
        if(lhs_val == 0.0 && corr_val == 0.0) { // Both didn't screen
            continue;
        } else if(lhs_val == 0.0) {
            std::cout << "We didn't screen, but libint did";
            return false;

        } else if(corr_val == 0.0) {
            std::cout << "We screened, but libint didn't";
            return false;
        }
        auto lhs_mag  = std::log10(std::fabs(lhs_val));
        auto corr_mag = std::log10(std::fabs(corr_val));
        if(lhs_mag != Catch::Approx(corr_mag).epsilon(0.2)) {
            auto pct_error = (lhs_mag - corr_mag) / corr_mag;
            std::cout << lhs_mag << " " << corr_mag << " " << pct_error
                      << std::endl;
            return false;
        }
    }
    return true;
}

} // namespace

TEST_CASE("PrimitiveErrorModel") {
    auto mm          = initialize_integrals();
    auto& mod        = mm.at("Primitive Error Model");
    auto& anal_error = mm.at("Analytic Error");
    simde::type::v_ee_type v_ee{};

    SECTION("Single (ss|ss) quartet") {
        auto [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();
        simde::type::aos_squared bra(bra0, bra1);
        simde::type::aos_squared ket(ket0, ket1);
        chemist::braket::BraKet mnls(bra, v_ee, ket);

        double tol = 1E-10;
        auto error = mod.run_as<pt>(mnls, tol);
        auto corr  = anal_error.run_as<pt>(mnls, tol);
        REQUIRE(compare_magnitudes(error, corr, tol));
    }

    SECTION("H2 molecule") {
        auto aobs = h2_sto3g_basis_set();
        simde::type::aos_squared bra(aobs, aobs);
        simde::type::aos_squared ket(aobs, aobs);
        chemist::braket::BraKet mnls(bra, v_ee, ket);

        double tol = 1E-6;
        auto error = mod.run_as<pt>(mnls, tol);
        auto corr  = anal_error.run_as<pt>(mnls, tol);
        REQUIRE(compare_magnitudes(error, corr, tol));
    }

    SECTION("Water/STO-3G(00|34)") {
        auto [bra0, bra1, ket0, ket1] = get_h2o_0034_bases();
        simde::type::aos_squared bra(bra0, bra1);
        simde::type::aos_squared ket(ket0, ket1);
        chemist::braket::BraKet mnls(bra, v_ee, ket);

        double tol = 1E-10;
        auto error = mod.run_as<pt>(mnls, tol);
        auto corr  = anal_error.run_as<pt>(mnls, tol);
        REQUIRE(compare_magnitudes(error, corr, tol));
    }

    SECTION("H2O molecule") {
        auto aobs = water_sto3g_basis_set();
        simde::type::aos_squared bra(aobs, aobs);
        simde::type::aos_squared ket(aobs, aobs);
        chemist::braket::BraKet mnls(bra, v_ee, ket);

        double tol = 1E-10;
        auto error = mod.run_as<pt>(mnls, tol);
        auto corr  = anal_error.run_as<pt>(mnls, tol);
        REQUIRE(compare_magnitudes(error, corr, tol));
    }
}
