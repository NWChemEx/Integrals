/*
 * Copyright 2026 NWChemEx-Project
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
#include "../testing/testing.hpp"
#include <integrals/integrals.hpp>

using namespace integrals::testing;
using pair_type   = std::vector<std::size_t>;
using pair_vector = std::vector<pair_type>;

namespace {

auto corr_bra_pairs() {
    pair_vector rv;
    rv.push_back({0, 2});
    rv.push_back({1, 1});
    rv.push_back({1, 2});
    rv.push_back({2, 0});
    rv.push_back({2, 1});
    rv.push_back({2, 2});
    return rv;
}

auto corr_ket_pairs() {
    pair_vector rv;
    rv.push_back({0, 1});
    rv.push_back({0, 2});
    rv.push_back({1, 0});
    rv.push_back({1, 1});
    rv.push_back({1, 2});
    rv.push_back({2, 0});
    rv.push_back({2, 1});
    rv.push_back({2, 2});
    return rv;
}

} // namespace

using pt = integrals::property_types::PairScreener;

TEST_CASE("Screen Primitive Pairs") {
    auto mm = initialize_integrals();

    auto mod = mm.at("Screen Primitive Pairs");

    auto [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();

    SECTION("(Bra0, Bra1) pairs") {
        auto pairs = mod.run_as<pt>(bra0, bra1, 1E-10);
        REQUIRE(pairs == corr_bra_pairs());
    }

    SECTION("(Ket0, Ket1) pairs") {
        auto pairs = mod.run_as<pt>(ket0, ket1, 1E-10);
        REQUIRE(pairs == corr_ket_pairs());
    }
}