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

#include "../testing/testing.hpp"

using namespace integrals::testing;

using namespace tensorwrapper;

namespace {

// N.b. The "means" of the correct values are validated by comparing to Libint's
// results. With the exception of the "No Mean" values, the
// "standard deviations" are not validated, but seem reasonable (the "No Mean"
// values come from the unit tests of the underlying error module).

template<typename FloatType>
auto corr_answer(const simde::type::tensor& T) {
    if constexpr(std::is_same_v<FloatType, double>) {
        return T;
    } else {
        simde::type::tensor T_corr(T);
        auto& corr_buffer = buffer::make_contiguous(T_corr.buffer());
        corr_buffer.set_elem({0, 0, 0, 0}, FloatType{0.774606, 0});
        corr_buffer.set_elem({0, 0, 0, 1},
                             FloatType{0.265558, 0.0000010000000000});
        corr_buffer.set_elem({0, 0, 1, 0},
                             FloatType{0.265558, 0.0000010000000000});
        corr_buffer.set_elem({0, 0, 1, 1}, FloatType{0.446701, 0});
        corr_buffer.set_elem({0, 1, 0, 0},
                             FloatType{0.265558, 0.0000010000000000});
        corr_buffer.set_elem({0, 1, 0, 1},
                             FloatType{0.120666, 0.0000170000000000});
        corr_buffer.set_elem({0, 1, 1, 0},
                             FloatType{0.120666, 0.0000170000000000});
        corr_buffer.set_elem({0, 1, 1, 1},
                             FloatType{0.265558, 0.0000010000000000});
        corr_buffer.set_elem({1, 0, 0, 0},
                             FloatType{0.265558, 0.0000010000000000});
        corr_buffer.set_elem({1, 0, 0, 1},
                             FloatType{0.120666, 0.0000170000000000});
        corr_buffer.set_elem({1, 0, 1, 0},
                             FloatType{0.120666, 0.0000170000000000});
        corr_buffer.set_elem({1, 0, 1, 1},
                             FloatType{0.265558, 0.0000010000000000});
        corr_buffer.set_elem({1, 1, 0, 0}, FloatType{0.446701, 0});
        corr_buffer.set_elem({1, 1, 0, 1},
                             FloatType{0.265558, 0.0000010000000000});
        corr_buffer.set_elem({1, 1, 1, 0},
                             FloatType{0.265558, 0.0000010000000000});
        corr_buffer.set_elem({1, 1, 1, 1}, FloatType{0.774606, 0});
        return T_corr;
    }
}

} // namespace

TEST_CASE("UQ Atom Blocked Driver") {
    using float_type = tensorwrapper::types::udouble;
    using test_pt    = simde::ERI4;

    if constexpr(tensorwrapper::types::is_uncertain_v<float_type>) {
        auto rt = std::make_unique<parallelzone::runtime::RuntimeView>();
        pluginplay::ModuleManager mm(std::move(rt), nullptr);
        integrals::load_modules(mm);
        integrals::set_defaults(mm);
        REQUIRE(mm.count("UQ Atom Blocked Driver"));
        mm.change_input("ERI4", "Threshold", 1.0e-6);
        auto mod = mm.at("UQ Atom Blocked Driver");

        // Get basis set
        auto aobs = h2_sto3g_basis_set();

        // Make AOS object
        simde::type::aos aos(aobs);
        simde::type::aos_squared aos_squared(aos, aos);

        // Make Operator
        simde::type::v_ee_type op{};

        // Make BraKet Input
        chemist::braket::BraKet braket(aos_squared, op, aos_squared);

        using tensorwrapper::operations::approximately_equal;
        SECTION("No Mean") { REQUIRE_THROWS(mod.run_as<test_pt>(braket)); }
        SECTION("Max Error") {
            mod.change_input("Mean Type", "max");
            auto T = mod.run_as<test_pt>(braket);

            auto T_corr = corr_answer<float_type>(T);
            REQUIRE(approximately_equal(T_corr, T, 1E-6));
        }
        SECTION("Geometric Mean") {
            mod.change_input("Mean Type", "geometric");
            auto T = mod.run_as<test_pt>(braket);

            auto T_corr = corr_answer<float_type>(T);
            REQUIRE(approximately_equal(T_corr, T, 1E-6));
        }
        SECTION("Harmonic Mean") {
            mod.change_input("Mean Type", "harmonic");
            auto T = mm.at("UQ Driver").run_as<test_pt>(braket);

            auto T_corr = corr_answer<float_type>(T);
            REQUIRE(approximately_equal(T_corr, T, 1E-6));
        }
    }
}
