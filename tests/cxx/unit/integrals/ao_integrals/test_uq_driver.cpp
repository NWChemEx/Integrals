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

#include "../testing.hpp"

TEST_CASE("UQ Driver") {
    using float_type = tensorwrapper::types::udouble;
    if constexpr(tensorwrapper::types::is_uncertain_v<float_type>) {
        using test_pt = simde::ERI4;

        auto rt = std::make_unique<parallelzone::runtime::RuntimeView>();
        pluginplay::ModuleManager mm(std::move(rt), nullptr);
        integrals::load_modules(mm);
        REQUIRE(mm.count("UQ Driver"));

        mm.change_input("UQ Driver", "precision", 1.0e-6);

        // Get basis set
        auto mol  = test::h2_molecule();
        auto aobs = test::h2_sto3g_basis_set();

        // Make AOS object
        simde::type::aos aos(aobs);
        simde::type::aos_squared aos_squared(aos, aos);

        // Make Operator
        simde::type::v_ee_type op{};

        // Make BraKet Input
        chemist::braket::BraKet braket(aos_squared, op, aos_squared);

        // Call module
        auto T = mm.at("UQ Driver").run_as<test_pt>(braket);

        simde::type::tensor T_corr(T);

        using alloc_type  = tensorwrapper::allocator::Eigen<float_type>;
        auto& corr_buffer = alloc_type::rebind(T_corr.buffer());
        corr_buffer.at(0, 0, 0, 0) = float_type{0.774606, 0};
        corr_buffer.at(0, 0, 0, 1) = float_type{0.265558, 2.49687e-06};
        corr_buffer.at(0, 0, 1, 0) = float_type{0.265558, 2.49687e-06};
        corr_buffer.at(0, 0, 1, 1) = float_type{0.446701, 0};
        corr_buffer.at(0, 1, 0, 0) = float_type{0.265558, 2.49687e-06};
        corr_buffer.at(0, 1, 0, 1) = float_type{0.120666, 1.10748e-05};
        corr_buffer.at(0, 1, 1, 0) = float_type{0.120666, 1.10748e-05};
        corr_buffer.at(0, 1, 1, 1) = float_type{0.265558, 2.49687e-06};
        corr_buffer.at(1, 0, 0, 0) = float_type{0.265558, 2.49687e-06};
        corr_buffer.at(1, 0, 0, 1) = float_type{0.120666, 1.10748e-05};
        corr_buffer.at(1, 0, 1, 0) = float_type{0.120666, 1.10748e-05};
        corr_buffer.at(1, 0, 1, 1) = float_type{0.265558, 2.49687e-06};
        corr_buffer.at(1, 1, 0, 0) = float_type{0.446701, 0};
        corr_buffer.at(1, 1, 0, 1) = float_type{0.265558, 2.49687e-06};
        corr_buffer.at(1, 1, 1, 0) = float_type{0.265558, 2.49687e-06};
        corr_buffer.at(1, 1, 1, 1) = float_type{0.774606, 0};

        using tensorwrapper::operations::approximately_equal;
        REQUIRE(approximately_equal(T_corr, T, 1E-6));
    }
}
