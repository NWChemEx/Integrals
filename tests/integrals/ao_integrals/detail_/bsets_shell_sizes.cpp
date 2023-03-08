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

#include "integrals/ao_integrals/detail_/bsets_shell_sizes.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>

using namespace mokup;

using integrals::ao_integrals::detail_::bsets_shell_sizes;

TEST_CASE("bsets_shell_sizes") {
    const auto name  = molecule::h2o;
    const auto bs    = basis_set::sto3g;
    auto aos         = get_bases(name, bs);
    const auto& bset = aos.basis_set();
    std::vector input{bset};

    std::vector<std::vector<std::size_t>> corr{{1, 1, 3, 1, 1}};

    auto result = bsets_shell_sizes(input);
    REQUIRE(result == corr);
}