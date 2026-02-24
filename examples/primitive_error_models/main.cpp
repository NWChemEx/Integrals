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

#include <integrals/integrals.hpp>
#include <simde/simde.hpp>

/* This example showcases how to:
 *
 *  1. Compute the analytic error in an ERI4 integral tensor owing to primitive
 *     pair screening.
 */

namespace {

// This makes a basis set for H2 (bond distance 1.40 a.u.) using STO-3G.
inline simde::type::ao_basis_set h2_sto3g_basis_set() {
    using ao_basis_t     = simde::type::ao_basis_set;
    using atomic_basis_t = simde::type::atomic_basis_set;
    using cg_t           = simde::type::contracted_gaussian;
    using point_t        = simde::type::point;
    using doubles_t      = std::vector<double>;

    point_t r0{0.0, 0.0, 0.0};
    point_t r1{0.0, 0.0, 1.40};

    doubles_t cs{0.1543289673, 0.5353281423, 0.4446345422};
    doubles_t es{3.425250914, 0.6239137298, 0.1688554040};
    cg_t cg0(cs.begin(), cs.end(), es.begin(), es.end(), r0);
    cg_t cg1(cs.begin(), cs.end(), es.begin(), es.end(), r1);
    atomic_basis_t h0("sto-3g", 1, r0);
    atomic_basis_t h1("sto-3g", 1, r1);
    h0.add_shell(chemist::ShellType::cartesian, 0, cg0);
    h1.add_shell(chemist::ShellType::cartesian, 0, cg1);

    ao_basis_t bs;
    bs.add_center(h0);
    bs.add_center(h1);
    return bs;
}

} // namespace

// Property types for the ERI4 and the error in the ERI4
using eri4_pt       = simde::ERI4;
using eri4_error_pt = integrals::property_types::Uncertainty<eri4_pt>;

int main(int argc, char* argv[]) {
    // Makes sure the environment doesn't go out of scope before the end.
    auto rt = std::make_unique<parallelzone::runtime::RuntimeView>();

    // Initializes a ModuleManager object with the integrals plugin
    pluginplay::ModuleManager mm(std::move(rt), nullptr);
    integrals::load_modules(mm);
    integrals::set_defaults(mm);

    // Modules for computing analytic error and estimating error
    auto& analytic_error_mod = mm.at("Analytic Error");
    auto& error_model        = mm.at("Primitive Error Model");

    // Makes: basis set, direct product of the basis set, and 1/r12 operator
    simde::type::aos aos(h2_sto3g_basis_set());
    simde::type::aos_squared aos2(aos, aos);
    simde::type::v_ee_type op{};

    // Make BraKet
    chemist::braket::BraKet mnls(aos2, op, aos2);

    // Compute the error by screening with tolerance "tol"
    double tol        = 1E-10;
    auto error        = analytic_error_mod.run_as<eri4_error_pt>(mnls, tol);
    auto approx_error = error_model.run_as<eri4_error_pt>(mnls, tol);

    std::cout << "Analytic error: " << error << std::endl;
    std::cout << "Estimated error: " << approx_error << std::endl;

    return 0;
}
