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

using pt = simde::UQInitializer;
using simde::UQFactory;
using simde::UQKind;

TEST_CASE("UQInitializer") {
    auto mm   = initialize_integrals();
    auto& mod = mm.at("UQ Initializer");

    SECTION("Default UQ Type/order matches compiled-in defaults") {
        auto factory = mod.run_as<pt>();
        REQUIRE(factory.kind() == UQKind::none);
    }

    SECTION("No UQ Initializer") {
        auto& no_uq = mm.at("No UQ Initializer");
        REQUIRE(no_uq.run_as<pt>().kind() == UQKind::none);
    }

    SECTION("Each UQ Type string produces a factory of the matching kind") {
        std::vector<std::pair<std::string, UQKind>> cases{
          {"uncertain", UQKind::uncertain},
          {"interval", UQKind::interval},
          {"affine", UQKind::affine},
          {"thresholded affine", UQKind::thresholded_affine}};

        for(const auto& [uq_type, kind] : cases) {
            auto copy = mod.unlocked_copy();
            copy.change_input("UQ Type", uq_type);
            auto factory = copy.run_as<pt>();
            REQUIRE(factory.kind() == kind);

            // factory(center, radius) should not throw and (when Sigma
            // supplies real UQ types) should reflect the requested
            // center/radius via the underlying UQ type's bounds.
            [[maybe_unused]] auto elem = factory(0.774606, 0.0000010000000000);
#ifdef ENABLE_SIGMA
            wtf::fp::visit_float<tensorwrapper::types::floating_point_types>(
              [](auto value) {
                  auto lo = tensorwrapper::types::uq_lower(value);
                  auto hi = tensorwrapper::types::uq_upper(value);
                  REQUIRE(lo <= 0.774606);
                  REQUIRE(hi >= 0.774606);
              },
              elem);
#endif
        }
    }

    SECTION("Invalid UQ Type throws") {
        auto copy = mod.unlocked_copy();
        copy.change_input("UQ Type", std::string("not a real uq type"));
        REQUIRE_THROWS_AS(copy.run_as<pt>(), std::runtime_error);
    }
}
