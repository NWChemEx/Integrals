/*
 * Copyright 2024 NWChemEx-Project
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

using simde::type::tensor;

namespace {

void compare_matrices(const tensor& A, const tensor& A_corr) {
    using alloc_type          = tensorwrapper::allocator::Eigen<double>;
    const auto& A_buffer      = alloc_type::rebind(A.buffer());
    const auto& A_corr_buffer = alloc_type::rebind(A_corr.buffer());

    const auto tol = 1E-6;
    auto A00       = A_buffer.get_elem({0, 0});
    auto A01       = A_buffer.get_elem({0, 1});
    auto A10       = A_buffer.get_elem({1, 0});
    auto A11       = A_buffer.get_elem({1, 1});

    REQUIRE(A00 == Catch::Approx(A_corr_buffer.get_elem({0, 0})).margin(tol));
    REQUIRE(A01 == Catch::Approx(A_corr_buffer.get_elem({0, 1})).margin(tol));
    REQUIRE(A10 == Catch::Approx(A_corr_buffer.get_elem({1, 0})).margin(tol));
    REQUIRE(A11 == Catch::Approx(A_corr_buffer.get_elem({1, 1})).margin(tol));
}

} // namespace

TEST_CASE("AOIntegralsDriver") {
    using pt = simde::aos_op_base_aos;
    using erased_type =
      chemist::braket::BraKet<simde::type::aos, simde::type::op_base_type,
                              simde::type::aos>;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    REQUIRE(mm.count("AO integral driver"));
    auto& mod = mm.at("AO integral driver");

    // Get basis set
    auto mol  = test::h2_molecule();
    auto aobs = test::h2_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);

    // Operator Inputs
    simde::type::electron e;
    auto rho = test::h2_density();

    SECTION("Calling Kinetic") {
        auto& tmod = mm.at("Kinetic");
        simde::type::t_e_type t_e(e);
        chemist::braket::BraKet braket(aos, t_e, aos);
        erased_type copy_braket(braket);
        const auto& T      = mod.run_as<pt>(copy_braket);
        const auto& T_corr = tmod.run_as<simde::aos_t_e_aos>(braket);
        compare_matrices(T, T_corr);
    }

    SECTION("Calling Electron-Nuclear Attraction") {
        auto& tmod = mm.at("Nuclear");
        simde::type::v_en_type v_en(e, mol.nuclei().as_nuclei());
        chemist::braket::BraKet braket(aos, v_en, aos);
        erased_type copy_braket(braket);
        const auto& V      = mod.run_as<pt>(copy_braket);
        const auto& V_corr = tmod.run_as<simde::aos_v_en_aos>(braket);
        compare_matrices(V, V_corr);
    }

    SECTION("Calling J Matrix") {
        auto& jmod = mm.at("Four center J builder");
        simde::type::j_e_type j_e(e, rho);
        chemist::braket::BraKet braket(aos, j_e, aos);
        erased_type copy_braket(braket);
        const auto& J      = mod.run_as<pt>(copy_braket);
        const auto& J_corr = jmod.run_as<simde::aos_j_e_aos>(braket);
        compare_matrices(J, J_corr);
    }

    SECTION("Calling K Matrix") {
        auto& kmod = mm.at("Four center K builder");
        simde::type::k_e_type k_e(e, rho);
        chemist::braket::BraKet braket(aos, k_e, aos);
        erased_type copy_braket(braket);
        const auto& K      = mod.run_as<pt>(copy_braket);
        const auto& K_corr = kmod.run_as<simde::aos_k_e_aos>(braket);
        compare_matrices(K, K_corr);
    }
}