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
#include <integrals/utils/uncertainty_reductions.hpp>
using namespace integrals::utils;

TEST_CASE("max_mean") {
    SECTION("Basic Test") {
        std::vector<double> values{1.0, 10.0, 100.0};
        auto result = max_mean(values);
        REQUIRE(result == Catch::Approx(100.0));
    }
    SECTION("Contains zero") {
        std::vector<double> values{1.0, 10.0, 100.0, 0.0};
        auto result = max_mean(values);
        REQUIRE(result == Catch::Approx(100.0));
    }
    SECTION("All zero") {
        std::vector<double> values{0.0, 0.0, 0.0};
        auto result = max_mean(values);
        REQUIRE(result == Catch::Approx(0.0));
    }
    SECTION("Empty") {
        std::vector<double> values{};
        auto result = max_mean(values);
        REQUIRE(result == Catch::Approx(0.0));
    }
}

TEST_CASE("geometric_mean") {
    SECTION("Basic Test") {
        std::vector<double> values{1.0, 10.0, 100.0};
        auto result = geometric_mean(values);
        REQUIRE(result == Catch::Approx(10.0));
    }

    SECTION("Skips zero") {
        std::vector<double> values{1.0, 10.0, 100.0, 0.0};
        auto result = geometric_mean(values);
        REQUIRE(result == Catch::Approx(10.0));
    }
    SECTION("All zero") {
        std::vector<double> values{0.0, 0.0, 0.0};
        auto result = geometric_mean(values);
        REQUIRE(result == Catch::Approx(0.0));
    }
    SECTION("No elements") {
        std::vector<double> values{};
        auto result = geometric_mean(values);
        REQUIRE(result == Catch::Approx(0.0));
    }
}

TEST_CASE("harmonic_mean") {
    SECTION("Basic Test") {
        std::vector<double> values{1.0, 2.0, 3.0};
        auto result = harmonic_mean(values);
        REQUIRE(result == Catch::Approx(1.63636363636));
    }

    SECTION("Skips zero") {
        std::vector<double> values{1.0, 2.0, 3.0, 0.0};
        auto result = harmonic_mean(values);
        REQUIRE(result == Catch::Approx(1.63636363636));
    }

    SECTION("All zero") {
        std::vector<double> values{0.0, 0.0, 0.0};
        auto result = harmonic_mean(values);
        REQUIRE(result == Catch::Approx(0.0));
    }

    SECTION("No elements") {
        std::vector<double> values{};
        auto result = harmonic_mean(values);
        REQUIRE(result == Catch::Approx(0.0));
    }
}

TEST_CASE("compute_mean") {
    std::vector<double> values{1.0, 2.0, 3.0};
    SECTION("Max") {
        auto result = compute_mean(mean_type::max, values);
        REQUIRE(result == Catch::Approx(3.0));
    }
    SECTION("Geometric") {
        auto result = compute_mean(mean_type::geometric, values);
        REQUIRE(result == Catch::Approx(1.81712059283));
    }
    SECTION("Harmonic") {
        auto result = compute_mean(mean_type::harmonic, values);
        REQUIRE(result == Catch::Approx(1.63636363636));
    }
    SECTION("None") {
        using except_t = std::runtime_error;
        REQUIRE_THROWS_AS(compute_mean(mean_type::none, values), except_t);
    }
}

TEST_CASE("mean_from_string") {
    SECTION("Max") {
        auto result = mean_from_string("max");
        REQUIRE(result == mean_type::max);
    }
    SECTION("Geometric") {
        auto result = mean_from_string("geometric");
        REQUIRE(result == mean_type::geometric);
    }
    SECTION("Harmonic") {
        auto result = mean_from_string("harmonic");
        REQUIRE(result == mean_type::harmonic);
    }
    SECTION("None") {
        auto result = mean_from_string("none");
        REQUIRE(result == mean_type::none);
    }
    SECTION("Invalid") {
        using except_t = std::runtime_error;
        REQUIRE_THROWS_AS(mean_from_string("invalid"), except_t);
    }
}
