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

template<typename FloatType>
auto corr_answer_no_mean(const simde::type::tensor& T,
                         const simde::type::tensor& T_eri,
                         const simde::type::tensor& T_err) {
    if constexpr(std::is_same_v<FloatType, double>) {
        return true;
    } else {
        auto t_uq  = eigen_tensor<4, FloatType>(T.buffer());
        auto t_eri = eigen_tensor<4, double>(T_eri.buffer());
        auto t_err = eigen_tensor<4, double>(T_err.buffer());

        for(Eigen::Index i = 0; i < 2; ++i)
            for(Eigen::Index j = 0; j < 2; ++j)
                for(Eigen::Index k = 0; k < 2; ++k)
                    for(Eigen::Index l = 0; l < 2; ++l) {
                        REQUIRE(t_uq(i, j, k, l).median() ==
                                Catch::Approx(t_eri(i, j, k, l)).margin(1E-6));
                        REQUIRE(t_uq(i, j, k, l).radius() ==
                                Catch::Approx(t_err(i, j, k, l)).margin(1E-6));
                    }
        return true;
    }
}

} // namespace

TEST_CASE("UQ Atom Symm Blocked Driver") {
    using float_type = tensorwrapper::types::idouble;
    using test_pt    = simde::ERI4;
    using error_pt   = integrals::property_types::Uncertainty<test_pt>;
    using tensorwrapper::operations::approximately_equal;

    if constexpr(tensorwrapper::types::is_uncertain_v<float_type>) {
        auto rt = std::make_unique<parallelzone::runtime::RuntimeView>();
        pluginplay::ModuleManager mm(std::move(rt));
        integrals::load_modules(mm);
        integrals::set_defaults(mm);
        REQUIRE(mm.count("UQ Atom Symm Blocked Driver"));
        mm.change_input("ERI4", "Threshold", 1.0e-6);
        auto mod = mm.at("UQ Atom Symm Blocked Driver");

        // Make Operator
        simde::type::v_ee_type op{};

        SECTION("H2") {
            // Get basis set
            auto aobs = h2_sto3g_basis_set();

            // Make AOS object
            simde::type::aos aos(aobs);
            simde::type::aos_squared aos_squared(aos, aos);

            // Make BraKet Input
            chemist::braket::BraKet braket(aos_squared, op, aos_squared);

            SECTION("No Mean") {
                auto T     = mod.run_as<test_pt>(braket);
                auto tol   = 1.0e-6;
                auto T_eri = mm.at("ERI4").run_as<test_pt>(braket);
                auto T_err =
                  mm.at("Primitive Error Model").run_as<error_pt>(braket, tol);
                REQUIRE(corr_answer_no_mean<float_type>(T, T_eri, T_err));
            }
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
                auto T = mod.run_as<test_pt>(braket);

                auto T_corr = corr_answer<float_type>(T);
                REQUIRE(approximately_equal(T_corr, T, 1E-6));
            }
        }

        SECTION("H2 Dimer") {
            simde::type::nucleus h0("H", 1ul, 1836.15, 0.0, 0.0, 0.0);
            simde::type::nucleus h1("H", 1ul, 1836.15, 0.0, 0.0, 1.39839);
            simde::type::nucleus h2("H", 1ul, 1836.15, 0.0, 0.0, 4.39839);
            simde::type::nucleus h3("H", 1ul, 1836.15, 0.0, 0.0, 5.79678);
            auto aobs = apply_h2_sto3g_basis_set(std::vector{h0, h1, h2, h3});

            // Make AOS object
            simde::type::aos aos(aobs);
            simde::type::aos_squared aos_squared(aos, aos);

            // Make BraKet Input
            chemist::braket::BraKet braket(aos_squared, op, aos_squared);
            const auto corr_key = "UQ Atom Blocked Driver";
            auto corr_mod       = mm.at(corr_key);
            SECTION("Max Error") {
                mod.change_input("Mean Type", "max");
                corr_mod.change_input("Mean Type", "max");
                auto eris      = mod.run_as<test_pt>(braket);
                auto eris_corr = corr_mod.run_as<test_pt>(braket);
                REQUIRE(approximately_equal(eris_corr, eris, 1E-6));
            }
        }
    }
}
