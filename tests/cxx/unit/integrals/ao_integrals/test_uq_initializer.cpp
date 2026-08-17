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

using namespace integrals::testing;

using pt = integrals::property_types::UQInitializer;
using integrals::property_types::TaylorModelFactory;

TEST_CASE("UQInitializer") {
    auto mm   = initialize_integrals();
    auto& mod = mm.at("UQ Initializer");

    SECTION("Default order matches compiled-in default (2)") {
        auto factory = mod.run_as<pt>();
        REQUIRE(factory.order() == 2);

        auto elem = factory(0.774606, 0.0000010000000000);
        REQUIRE(elem.max_order() == 2);
    }

    SECTION("Custom order is honored") {
        auto mod4 = mod.unlocked_copy();
        mod4.change_input("Order", std::size_t(4));
        auto factory = mod4.run_as<pt>();
        REQUIRE(factory.order() == 4);

        auto elem = factory(0.774606, 0.0000010000000000);
        REQUIRE(elem.max_order() == 4);
    }
}

TEST_CASE("UQInitializer::TaylorModelFactory") {
    SECTION("operator== / operator< satisfy AnyField comparability") {
        TaylorModelFactory f2(2), f2b(2), f4(4);
        REQUIRE(f2 == f2b);
        REQUIRE_FALSE(f2 == f4);
        REQUIRE(f2 < f4);
        REQUIRE_FALSE(f4 < f2);
    }

    SECTION("order() reports the configured order") {
        TaylorModelFactory f; // default order (2)
        TaylorModelFactory f3(3);
        REQUIRE(f.order() == 2);
        REQUIRE(f3.order() == 3);
    }
}
