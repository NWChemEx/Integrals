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

#include "integrals/ao_integrals/detail_/aos2shells.hpp"
#include <catch2/catch.hpp>

using integrals::ao_integrals::detail_::aos2shells;

TEST_CASE("aos2shells") {
    std::vector<std::size_t> shell_sizes{1, 1, 3, 1, 1};

    // Run with different inputs
    auto all     = aos2shells(shell_sizes, 0, 7);
    auto only_o  = aos2shells(shell_sizes, 0, 5);
    auto only_h1 = aos2shells(shell_sizes, 5, 6);
    auto only_h2 = aos2shells(shell_sizes, 6, 7);
    auto both_hs = aos2shells(shell_sizes, 5, 7);

    /// Check outputs
    REQUIRE(all == std::vector<std::size_t>{0, 1, 2, 3, 4});
    REQUIRE(only_o == std::vector<std::size_t>{0, 1, 2});
    REQUIRE(only_h1 == std::vector<std::size_t>{3});
    REQUIRE(only_h2 == std::vector<std::size_t>{4});
    REQUIRE(both_hs == std::vector<std::size_t>{3, 4});
}