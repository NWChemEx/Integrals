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

using pt = integrals::property_types::PairScreener;

TEST_CASE("Screen Primitive Pairs") {
    auto mm = initialize_integrals();

    auto mod = mm.at("Screen Primitive Pairs");

    auto [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();

    SECTION("(Bra0, Bra1) pairs") {
        auto pairs = mod.run_as<pt>(bra0, bra1, 1E-10);
        for(const auto& pair : pairs) {
            for(const auto& prim : pair) { std::cout << prim << " "; }
            std::cout << std::endl;
        }
    }
}