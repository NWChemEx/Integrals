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

#include "integrals/libint/detail_/hash_inputs.hpp"
#include "libint_basis_set_water.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>

using namespace mokup;
using namespace integrals::detail_;

using bases_vector_t = std::vector<libint2::BasisSet>;

TEST_CASE("Combining hashes") {
    /// Values for input and comparison
    const double v     = 1.23456;
    std::size_t v_hash = std::hash<double>{}(v);
    std::size_t corr1  = 0 ^ (v_hash + 0x9e3779b9 + (0 << 6) + (0 >> 2));
    std::size_t corr2 =
      corr1 ^ (v_hash + 0x9e3779b9 + (corr1 << 6) + (corr1 >> 2));

    SECTION("hash_together") {
        std::size_t seed = 0;
        hash_together(seed, v_hash);
        REQUIRE(seed == corr1);
        hash_together(seed, v_hash);
        REQUIRE(seed == corr2);
    }

    SECTION("combine_hash") {
        std::size_t seed = 0;
        combine_hash(seed, v);
        REQUIRE(seed == corr1);
        combine_hash(seed, v);
        REQUIRE(seed == corr2);

        /// Reset seed and check if it works as one call
        seed = 0;
        combine_hash(seed, v, v);
        REQUIRE(seed == corr2);
    }
}

TEST_CASE("Hashing Operators") {
    simde::type::el_el_coulomb op;

    auto hash = hash_operator(op);
    /// TODO: Validate the output
}

TEST_CASE("Hashing Bases") {
    auto bset = testing::water_basis_set();
    bases_vector_t two_sets{bset, bset};
    auto hash = std::hash<bases_vector_t>{}(two_sets);
    /// TODO: Validate the output
}

TEST_CASE("hash_inputs") {
    auto bset = testing::water_basis_set();

    bases_vector_t two_sets{bset, bset};
    bases_vector_t three_sets{bset, bset, bset};
    bases_vector_t four_sets{bset, bset, bset, bset};

    double t = 1.23456;

    simde::type::el_el_coulomb op;

    auto hash1 = hash_inputs(two_sets, op, t);
    auto hash2 = hash_inputs(three_sets, op, t);
    auto hash3 = hash_inputs(four_sets, op, t);
    /// TODO: Validate the output
}