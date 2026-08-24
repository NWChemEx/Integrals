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
using integrals::property_types::uq_kind_from_string;
using integrals::property_types::UQFactory;
using integrals::property_types::UQKind;

TEST_CASE("UQInitializer") {
    auto mm   = initialize_integrals();
    auto& mod = mm.at("UQ Initializer");

    SECTION("Default UQ Type/order matches compiled-in defaults") {
        auto factory = mod.run_as<pt>();
        REQUIRE(factory.kind() == UQKind::uncertain);
    }

    SECTION("Each UQ Type string produces a factory of the matching kind") {
        std::vector<std::pair<std::string, UQKind>> cases{
          {"uncertain", UQKind::uncertain},
          {"interval", UQKind::interval},
          {"affine", UQKind::affine},
          {"thresholded affine", UQKind::thresholded_affine},
          {"taylor model", UQKind::taylor_model}};

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

    SECTION("Custom order is honored for taylor model") {
        auto mod4 = mod.unlocked_copy();
        mod4.change_input("UQ Type", std::string("taylor model"));
        mod4.change_input("Order", std::size_t(4));
        auto factory = mod4.run_as<pt>();
        REQUIRE(factory.kind() == UQKind::taylor_model);

        [[maybe_unused]] auto elem = factory(0.774606, 0.0000010000000000);
        // Without Sigma the "taylor model" type is a plain double, which has
        // no order to honor.
#ifdef ENABLE_SIGMA
        auto tm =
          wtf::fp::float_cast<tensorwrapper::types::taylor_model_type<double>>(
            elem);
        REQUIRE(tm.max_order() == 4);
#endif
    }

    SECTION("Invalid UQ Type throws") {
        auto copy = mod.unlocked_copy();
        copy.change_input("UQ Type", std::string("not a real uq type"));
        REQUIRE_THROWS_AS(copy.run_as<pt>(), std::runtime_error);
    }
}
