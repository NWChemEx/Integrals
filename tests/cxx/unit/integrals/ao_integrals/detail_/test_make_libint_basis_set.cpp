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

#include "integrals/ao_integrals/detail_/make_libint_basis_set.hpp"
#include <catch2/catch_test_macros.hpp>

TEST_CASE("make_libint_basis_set") {
    using atomic_basis_t = simde::type::atomic_basis_set;
    using cg_t           = simde::type::contracted_gaussian;
    using doubles_t      = std::vector<double>;

    doubles_t cs0{1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01};
    doubles_t es0{1.3070932140e+02, 2.3808866050e+01, 6.4436083130e+00};
    cg_t cg0(cs0.begin(), cs0.end(), es0.begin(), es0.end(), 0.0, 0.0, 0.0);
    doubles_t cs1{-9.9967229190e-02, 3.9951282610e-01, 7.0011546890e-01};
    doubles_t es1{5.0331513190e+00, 1.1695961250e+00, 3.8038896000e-01};
    cg_t cg1(cs1.begin(), cs1.end(), es1.begin(), es1.end(), 0.0, 0.0, 0.0);
    doubles_t cs2{1.5591627500e-01, 6.0768371860e-01, 3.9195739310e-01};
    doubles_t es2{5.0331513190e+00, 1.1695961250e+00, 3.8038896000e-01};
    cg_t cg2(cs2.begin(), cs2.end(), es2.begin(), es2.end(), 0.0, 0.0, 0.0);

    atomic_basis_t corr("sto-3g", 8, 0.0, 0.0, 0.0);
    corr.add_shell(chemist::ShellType::pure, 0, cg0);
    corr.add_shell(chemist::ShellType::pure, 0, cg1);
    corr.add_shell(chemist::ShellType::pure, 1, cg2);
}
