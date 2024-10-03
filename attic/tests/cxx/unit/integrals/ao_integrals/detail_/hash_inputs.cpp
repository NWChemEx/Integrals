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

#include "integrals/ao_integrals/detail_/hash_inputs.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>

using namespace mokup;
using namespace integrals::ao_integrals::detail_;

using bases_vector_t = std::vector<simde::type::ao_basis_set>;

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
    SECTION("el_el_coulomb") {
        simde::type::el_el_coulomb op;
        simde::type::el_el_delta other_op;
        auto op_hash = hash_operator(op);
        auto corr    = std::hash<std::string>{}("(r₁₂)⁻¹");

        /// Check correctness of hash
        REQUIRE(op_hash == corr);
        /// Check differentiation of hashes
        REQUIRE(op_hash != hash_operator(other_op));
    }

    SECTION("el_nuc_coulomb") {
        using e_t    = chemist::Electron;
        using atom_t = chemist::Atom;
        using nuc_t  = chemist::Nuclei;
        using op_t   = simde::type::el_nuc_coulomb;

        atom_t a1{"X", std::size_t{1}, 0.0, 0.0, 0.0, 0.0};
        atom_t a2{"X", std::size_t{2}, 0.0, 1.0, 0.0, 0.0};
        op_t op(e_t{}, nuc_t{a1.nucleus()});
        op_t other_op1(e_t{}, nuc_t{a2.nucleus()});
        op_t other_op2(e_t{}, nuc_t{a1.nucleus(), a2.nucleus()});
        auto op_hash = hash_operator(op);

        auto corr = std::hash<std::string>{}("(r₁₂)⁻¹");
        combine_hash(corr, std::size_t{1}, std::string("X"), 0.0, 1.0, 0.0, 0.0,
                     0.0);

        /// Check correctness of hash
        REQUIRE(op_hash == corr);
        /// Check differentiation of hashes
        REQUIRE(op_hash != hash_operator(other_op1));
        REQUIRE(op_hash != hash_operator(other_op2));
    }

    SECTION("el_el_stg") {
        using stg_t = chemist::operators::STG;
        using op_t  = simde::type::el_el_stg;

        op_t op(stg_t{1.0, 1.0});
        op_t other_op1(stg_t{2.0, 1.0});
        op_t other_op2(stg_t{1.0, 2.0});
        auto op_hash = hash_operator(op);

        auto corr = std::hash<std::string>{}("f₁₂");
        combine_hash(corr, 1.0, 1.0);

        /// Check correctness of hash
        REQUIRE(op_hash == corr);
        /// Check differentiation of hashes
        REQUIRE(op_hash != hash_operator(other_op1));
        REQUIRE(op_hash != hash_operator(other_op2));
    }
}

TEST_CASE("Hashing Bases") {
    const auto name = molecule::h2o;
    auto sto3g_nwx  = get_bases(name, basis_set::sto3g);
    auto ccpvdz_nwx = get_bases(name, basis_set::ccpvdz);

    auto sto3g  = sto3g_nwx.basis_set();
    auto ccpvdz = ccpvdz_nwx.basis_set();

    bases_vector_t set1{sto3g};
    bases_vector_t set2{ccpvdz};
    bases_vector_t set3{sto3g, sto3g};
    bases_vector_t set4{sto3g, ccpvdz};
    bases_vector_t set5{ccpvdz, sto3g};

    std::hash<bases_vector_t> hasher;
    auto hash1 = hasher(set1);
    auto hash2 = hasher(set2);
    auto hash3 = hasher(set3);
    auto hash4 = hasher(set4);
    auto hash5 = hasher(set5);

    /// Checking the absolute correctness of the hash is overly involved,
    /// so I'm just checking for differentiation between instances.
    REQUIRE(hash1 != hash2); // Same number, different contents
    REQUIRE(hash1 != hash3); // Different number of the same set
    REQUIRE(hash4 != hash5); // Check asymetry of hashing
}

/// Everything else is tested by this point, so just make sure that this
/// is consistent with expectation.
TEST_CASE("hash_inputs") {
    /// Inputs
    simde::type::el_el_coulomb op;
    const auto name = molecule::h2o;
    auto sto3g_nwx  = get_bases(name, basis_set::sto3g);
    auto sto3g      = sto3g_nwx.basis_set();
    bases_vector_t bases{sto3g};
    double t = 1.23456;

    /// Hash inputs
    auto hash = hash_inputs(bases, op, t);

    /// Correct value
    auto corr = hash_operator(op);
    combine_hash(corr, bases, t);

    /// Check correctness of hash
    REQUIRE(hash == std::to_string(corr));
}
