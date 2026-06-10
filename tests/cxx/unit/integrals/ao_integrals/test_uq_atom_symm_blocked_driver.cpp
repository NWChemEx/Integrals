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

template<typename UQType, typename FloatType>
auto elem(FloatType value, FloatType error) {
    if constexpr(tensorwrapper::types::is_uncertain_v<UQType>) {
        return UQType{value, error};
    } else if constexpr(tensorwrapper::types::is_interval_v<UQType>) {
        return UQType{value - error, value + error};
    } else if constexpr(tensorwrapper::types::is_affine_v<UQType>) {
        return UQType{value - error, value + error};
    } else {
        ::utilities::printing::Demangler demangler;
        auto type0 = demangler.template demangle<UQType>();
        throw std::runtime_error(
          "UQ Atom Symm Blocked Driver Test: Invalid UQ type " + type0);
    }
}

// N.b. The "means" of the correct values are validated by comparing to Libint's
// results. With the exception of the "No Mean" values, the
// "standard deviations" are not validated, but seem reasonable (the "No Mean"
// values come from the unit tests of the underlying error module).

template<typename UQType>
auto corr_answer(const simde::type::tensor& T) {
    if constexpr(std::is_same_v<UQType, double>) {
        return T;
    } else {
        using value_t = typename UQType::value_t;
        simde::type::tensor T_corr(T);
        auto& corr_buffer = buffer::make_contiguous(T_corr.buffer());
        corr_buffer.set_elem({0, 0, 0, 0}, elem<UQType>(0.774606, value_t(0)));
        corr_buffer.set_elem({0, 0, 0, 1},
                             elem<UQType>(0.265558, 0.0000010000000000));
        corr_buffer.set_elem({0, 0, 1, 0},
                             elem<UQType>(0.265558, 0.0000010000000000));
        corr_buffer.set_elem({0, 0, 1, 1}, elem<UQType>(0.446701, value_t(0)));
        corr_buffer.set_elem({0, 1, 0, 0},
                             elem<UQType>(0.265558, 0.0000010000000000));
        corr_buffer.set_elem({0, 1, 0, 1},
                             elem<UQType>(0.120666, 0.0000170000000000));
        corr_buffer.set_elem({0, 1, 1, 0},
                             elem<UQType>(0.120666, 0.0000170000000000));
        corr_buffer.set_elem({0, 1, 1, 1},
                             elem<UQType>(0.265558, 0.0000010000000000));
        corr_buffer.set_elem({1, 0, 0, 0},
                             elem<UQType>(0.265558, 0.0000010000000000));
        corr_buffer.set_elem({1, 0, 0, 1},
                             elem<UQType>(0.120666, 0.0000170000000000));
        corr_buffer.set_elem({1, 0, 1, 0},
                             elem<UQType>(0.120666, 0.0000170000000000));
        corr_buffer.set_elem({1, 0, 1, 1},
                             elem<UQType>(0.265558, 0.0000010000000000));
        corr_buffer.set_elem({1, 1, 0, 0}, elem<UQType>(0.446701, value_t(0)));
        corr_buffer.set_elem({1, 1, 0, 1},
                             elem<UQType>(0.265558, 0.0000010000000000));
        corr_buffer.set_elem({1, 1, 1, 0},
                             elem<UQType>(0.265558, 0.0000010000000000));
        corr_buffer.set_elem({1, 1, 1, 1}, elem<UQType>(0.774606, value_t(0)));
        return T_corr;
    }
}

template<typename UQType>
auto corr_answer_no_mean(const simde::type::tensor& T,
                         const simde::type::tensor& T_eri,
                         const simde::type::tensor& T_err) {
    if constexpr(std::is_same_v<UQType, double>) {
        return true;
    } else {
        auto t_uq  = eigen_tensor<4, UQType>(T.buffer());
        auto t_eri = eigen_tensor<4, double>(T_eri.buffer());
        auto t_err = eigen_tensor<4, double>(T_err.buffer());

        for(Eigen::Index i = 0; i < 2; ++i)
            for(Eigen::Index j = 0; j < 2; ++j)
                for(Eigen::Index k = 0; k < 2; ++k)
                    for(Eigen::Index l = 0; l < 2; ++l) {
                        if constexpr(tensorwrapper::types::is_uncertain_v<
                                       UQType>) {
                            REQUIRE(
                              t_uq(i, j, k, l).mean() ==
                              Catch::Approx(t_eri(i, j, k, l)).margin(1E-6));
                            REQUIRE(
                              t_uq(i, j, k, l).sd() ==
                              Catch::Approx(t_err(i, j, k, l)).margin(1E-6));
                        }
                    }
        return true;
    }
}

} // namespace

TEST_CASE("UQ Atom Symm Blocked Driver") {
    using uncertain_type = tensorwrapper::types::udouble;
    using interval_type  = tensorwrapper::types::interval_type<double>;
    using test_pt        = simde::ERI4;
    using error_pt       = integrals::property_types::Uncertainty<test_pt>;
    using tensorwrapper::operations::approximately_equal;

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

        if constexpr(tensorwrapper::types::is_uncertain_v<uncertain_type>) {
            SECTION("No Mean") {
                auto T     = mod.run_as<test_pt>(braket);
                auto tol   = 1.0e-6;
                auto T_eri = mm.at("ERI4").run_as<test_pt>(braket);
                auto T_err =
                  mm.at("Primitive Error Model").run_as<error_pt>(braket, tol);
                REQUIRE(corr_answer_no_mean<uncertain_type>(T, T_eri, T_err));
            }
            SECTION("Max Error") {
                mod.change_input("Mean Type", "max");
                auto T = mod.run_as<test_pt>(braket);

                auto T_corr = corr_answer<uncertain_type>(T);
                REQUIRE(approximately_equal(T_corr, T, 1E-6));
            }
            SECTION("Geometric Mean") {
                mod.change_input("Mean Type", "geometric");
                auto T = mod.run_as<test_pt>(braket);

                auto T_corr = corr_answer<uncertain_type>(T);
                REQUIRE(approximately_equal(T_corr, T, 1E-6));
            }
        }
        if constexpr(tensorwrapper::types::is_interval_v<interval_type> ||
                     tensorwrapper::types::is_affine_v<interval_type>) {
            // The errors for the integrals are on the order of 1e-5.
            // Subtracting
            // the intervals in approximately_equal treats the intervals as
            // independent and doubles them instead of cancelling them. This in
            // turn means center + radius on the difference inside
            // approximately_equal will be between 1e-5 and 1e-4.
            SECTION("No Mean") {
                mod.change_input("UQ Type", "interval");
                auto T     = mod.run_as<test_pt>(braket);
                auto tol   = 1.0e-6;
                auto T_eri = mm.at("ERI4").run_as<test_pt>(braket);
                auto T_err =
                  mm.at("Primitive Error Model").run_as<error_pt>(braket, tol);
                REQUIRE(corr_answer_no_mean<interval_type>(T, T_eri, T_err));
            }
            SECTION("Max Error") {
                mod.change_input("UQ Type", "interval");
                mod.change_input("Mean Type", "max");
                auto T = mod.run_as<test_pt>(braket);

                auto T_corr = corr_answer<interval_type>(T);
                REQUIRE(approximately_equal(T_corr, T, 1E-4));
            }
            SECTION("Geometric Mean") {
                mod.change_input("UQ Type", "interval");
                mod.change_input("Mean Type", "geometric");
                auto T = mod.run_as<test_pt>(braket);

                auto T_corr = corr_answer<interval_type>(T);
                REQUIRE(approximately_equal(T_corr, T, 1E-4));
            }
        }
    }
}
