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

#include "integrals/libint/detail_/shells2ord.hpp"
#include "libint_basis_set_water.hpp"
#undef DEPRECATED
#include <catch2/catch_test_macros.hpp>

using namespace integrals;

TEST_CASE("shells2ord") {
    using libint::detail_::shells2ord;
    auto aobs = testing::water_basis_set();
    std::vector<libint2::BasisSet> basis_sets{aobs, aobs};
    std::vector<std::size_t> shells{2, 2};
    auto out                      = shells2ord(basis_sets, shells);
    std::vector<std::size_t> corr = {16, 17, 18, 23, 24, 25, 30, 31, 32};
    REQUIRE(out == corr);
}
