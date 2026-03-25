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

#include "integrals/libint/detail_/make_libint_basis_set.hpp"
#include "libint_basis_set_water.hpp"
#undef DEPRECATED
// Must be last due to conflicting macros
#include "../../testing/testing.hpp"

using namespace integrals::testing;

TEST_CASE("make_libint_basis_set") {
    using integrals::libint::detail_::make_libint_basis_set;
    auto aobs = water_sto3g_basis_set();

    SECTION("embed_normalization=true (default)") {
        auto libint_bs   = make_libint_basis_set(aobs);
        auto libint_corr = test::water_basis_set();
        REQUIRE(libint_bs == libint_corr);
    }

    SECTION("embed_normalization=false") {
        auto libint_norm   = make_libint_basis_set(aobs, true);
        auto libint_nonorm = make_libint_basis_set(aobs, false);

        // The two basis sets should have the same number of shells
        REQUIRE(libint_norm.size() == libint_nonorm.size());

        // With normalization disabled the coefficients must differ from the
        // renormalized ones for at least one shell (the p shell on oxygen
        // is the clearest case since its normalization factor != 1)
        bool any_differ = false;
        for(std::size_t s = 0; s < libint_norm.size(); ++s) {
            const auto& c_norm   = libint_norm[s].contr[0].coeff;
            const auto& c_nonorm = libint_nonorm[s].contr[0].coeff;
            REQUIRE(c_norm.size() == c_nonorm.size());
            for(std::size_t p = 0; p < c_norm.size(); ++p) {
                if(c_norm[p] != c_nonorm[p]) {
                    any_differ = true;
                    break;
                }
            }
            if(any_differ) break;
        }
        REQUIRE(any_differ);
    }
}
