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

#include "detail_/libint_basis_set_water.hpp"
#include "integrals/integrals.hpp"
#include "integrals/libint/libint_factory.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>
#include <simde/integral_factory.hpp>
#include <simde/types.hpp>

using namespace mokup;

TEST_CASE("Make Overlap LibintFactory") {
    using op_t         = simde::type::el_identity;
    using factory_t    = simde::type::integral_factory;
    using libint_fac_t = integrals::libint::LibintFactory<2, op_t>;
    using test_pt      = simde::IntegralFactory<op_t>;

    /// Load up modules
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    /// Inputs
    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto aos        = get_bases(name, bs);
    auto basis      = aos.basis_set();
    std::vector sets{basis, basis};
    op_t op;

    /// Correct result
    auto libint_basis = testing::water_basis_set();
    std::vector libint_sets{libint_basis, libint_basis};
    auto p =
      std::make_unique<libint_fac_t>(std::move(libint_sets), op, 1.0E-16, 0);
    factory_t corr(std::move(p));

    /// Check module output
    auto fac = mm.at("Overlap Factory").run_as<test_pt>(sets, op);
    REQUIRE(fac == corr);
}
