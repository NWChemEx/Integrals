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

#include "../../testing.hpp"
#include "integrals/libint/detail_/make_libint_basis_set.hpp"
#include "libint_basis_set_water.hpp"

TEST_CASE("make_libint_basis_set") {
    using integrals::libint::detail_::make_libint_basis_set;
    auto aobs        = test::water_sto3g_basis_set();
    auto libint_bs   = make_libint_basis_set(aobs);
    auto libint_corr = test::water_basis_set();
    REQUIRE(libint_bs == libint_corr);
}
