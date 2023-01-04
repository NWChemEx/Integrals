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

#include "integrals/libint/detail_/shells2ord.hpp"
#include "libint_basis_set_water.hpp"
#include <catch2/catch.hpp>

TEST_CASE("shell2ord") {
    using size_vector_t  = std::vector<std::size_t>;
    using basis_vector_t = std::vector<libint2::BasisSet>;

    /// Common basis set
    auto bset = testing::water_basis_set();

    /// Check different dimensionalities
    SECTION("2D") {
        basis_vector_t sets{bset, bset};
        size_vector_t curr{2, 0}, lo{0, 0}, up{4, 4};
        auto ord_pos = integrals::detail_::shells2ord(sets, curr, lo, up);
        REQUIRE(ord_pos == size_vector_t{14, 21, 28});
    }

    SECTION("3D") {
        basis_vector_t sets{bset, bset, bset};
        size_vector_t curr{2, 0, 0}, lo{0, 0, 0}, up{4, 4, 4};
        auto ord_pos = integrals::detail_::shells2ord(sets, curr, lo, up);
        REQUIRE(ord_pos == size_vector_t{98, 147, 196});
    }

    SECTION("4D") {
        basis_vector_t sets{bset, bset, bset, bset};
        size_vector_t curr{2, 0, 0, 0}, lo{0, 0, 0, 0}, up{4, 4, 4, 4};
        auto ord_pos = integrals::detail_::shells2ord(sets, curr, lo, up);
        REQUIRE(ord_pos == size_vector_t{686, 1029, 1372});
    }

    SECTION("specific tile") {
        basis_vector_t sets{bset, bset};
        size_vector_t lo{2, 0}, up{3, 0};
        SECTION("lower") {
            auto ord_pos = integrals::detail_::shells2ord(sets, lo, lo, up);
            REQUIRE(ord_pos == size_vector_t{0, 1, 2});
        }
        SECTION("upper") {
            auto ord_pos = integrals::detail_::shells2ord(sets, up, lo, up);
            REQUIRE(ord_pos == size_vector_t{3});
        }
    }
}