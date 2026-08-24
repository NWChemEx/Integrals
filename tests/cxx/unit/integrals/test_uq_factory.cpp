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

#include "testing/testing.hpp"
#include <integrals/uq_factory.hpp>

using namespace integrals::property_types;

namespace {

/// The center/radius used throughout; the radius is on the order of the error
/// the UQ-aware integral drivers actually see.
constexpr double center = 0.774606;
constexpr double radius = 0.0000010000000000;

/// All of the kinds, in the order they're declared (which is also the order
/// UQFactory::operator< sorts them in).
const std::vector<UQKind> all_kinds{UQKind::uncertain, UQKind::interval,
                                    UQKind::affine, UQKind::thresholded_affine,
                                    UQKind::taylor_model};

} // namespace

TEST_CASE("uq_kind_from_string") {
    SECTION("Valid names") {
        REQUIRE(uq_kind_from_string("uncertain") == UQKind::uncertain);
        REQUIRE(uq_kind_from_string("interval") == UQKind::interval);
        REQUIRE(uq_kind_from_string("affine") == UQKind::affine);
        REQUIRE(uq_kind_from_string("thresholded affine") ==
                UQKind::thresholded_affine);
        REQUIRE(uq_kind_from_string("taylor model") == UQKind::taylor_model);
    }

    SECTION("Invalid names throw") {
        REQUIRE_THROWS_AS(uq_kind_from_string("not a uq type"),
                          std::runtime_error);
        REQUIRE_THROWS_AS(uq_kind_from_string(""), std::runtime_error);

        // The mapping is case- and spelling-sensitive.
        REQUIRE_THROWS_AS(uq_kind_from_string("Uncertain"), std::runtime_error);
        REQUIRE_THROWS_AS(uq_kind_from_string("thresholded_affine"),
                          std::runtime_error);
        REQUIRE_THROWS_AS(uq_kind_from_string("taylor_model"),
                          std::runtime_error);
    }
}

TEST_CASE("to_string(UQKind)") {
    SECTION("Each kind") {
        REQUIRE(to_string(UQKind::uncertain) == "uncertain");
        REQUIRE(to_string(UQKind::interval) == "interval");
        REQUIRE(to_string(UQKind::affine) == "affine");
        REQUIRE(to_string(UQKind::thresholded_affine) == "thresholded affine");
        REQUIRE(to_string(UQKind::taylor_model) == "taylor model");
    }

    SECTION("Is the inverse of uq_kind_from_string") {
        for(auto kind : all_kinds)
            REQUIRE(uq_kind_from_string(to_string(kind)) == kind);
    }
}

TEST_CASE("UQFactory") {
    UQFactory defaulted;
    UQFactory uncertain(UQKind::uncertain);
    UQFactory interval(UQKind::interval);
    UQFactory affine(UQKind::affine);
    UQFactory taffine(UQKind::thresholded_affine);
    UQFactory tm(UQKind::taylor_model);
    UQFactory tm4(UQKind::taylor_model, 4);

    SECTION("Default ctor") { REQUIRE(defaulted.kind() == UQKind::uncertain); }

    SECTION("Value ctor") {
        for(auto kind : all_kinds) REQUIRE(UQFactory(kind).kind() == kind);
    }

    SECTION("Value ctor ignores order for non-TaylorModel kinds") {
        for(auto kind : all_kinds) {
            if(kind == UQKind::taylor_model) continue;
            REQUIRE(UQFactory(kind, 7) == UQFactory(kind));
        }
    }

    SECTION("Copy ctor") {
        UQFactory copy(tm4);
        REQUIRE(copy == tm4);
        REQUIRE(copy.kind() == UQKind::taylor_model);
    }

    SECTION("Copy assignment") {
        UQFactory copy;
        copy = interval;
        REQUIRE(copy == interval);
    }

// Sigma provides the UQ representations; without it TensorWrapper's UQ types
// are plain doubles that carry neither bounds nor a TaylorModel order, leaving
// nothing here to check.
#ifdef ENABLE_SIGMA
    SECTION("operator()") {
        using tensorwrapper::types::uq_center;
        using tensorwrapper::types::uq_lower;
        using tensorwrapper::types::uq_upper;

        SECTION("uncertain") {
            auto erased = uncertain(center, radius);
            auto value =
              wtf::fp::float_cast<tensorwrapper::types::uncertain_type<double>>(
                erased);
            REQUIRE(uq_center(value) == Catch::Approx(center));
            REQUIRE(value.sd() == Catch::Approx(radius));
        }

        SECTION("interval") {
            auto erased = interval(center, radius);
            auto value =
              wtf::fp::float_cast<tensorwrapper::types::interval_type<double>>(
                erased);
            REQUIRE(uq_lower(value) == Catch::Approx(center - radius));
            REQUIRE(uq_upper(value) == Catch::Approx(center + radius));
        }

        SECTION("affine") {
            auto erased = affine(center, radius);
            auto value =
              wtf::fp::float_cast<tensorwrapper::types::affine_type<double>>(
                erased);
            REQUIRE(uq_center(value) == Catch::Approx(center));
            REQUIRE(uq_lower(value) == Catch::Approx(center - radius));
            REQUIRE(uq_upper(value) == Catch::Approx(center + radius));
        }

        SECTION("thresholded affine") {
            auto erased = taffine(center, radius);
            auto value  = wtf::fp::float_cast<
               tensorwrapper::types::thresholded_affine_type<double>>(erased);
            REQUIRE(uq_center(value) == Catch::Approx(center));
            REQUIRE(uq_lower(value) == Catch::Approx(center - radius));
            REQUIRE(uq_upper(value) == Catch::Approx(center + radius));
        }

        SECTION("taylor model") {
            using tm_type = tensorwrapper::types::taylor_model_type<double>;
            auto erased   = tm(center, radius);
            auto value    = wtf::fp::float_cast<tm_type>(erased);
            REQUIRE(uq_center(value) == Catch::Approx(center));
            REQUIRE(uq_lower(value) == Catch::Approx(center - radius));
            REQUIRE(uq_upper(value) == Catch::Approx(center + radius));

            // The order defaults to 2 and is settable.
            REQUIRE(value.max_order() == 2);
            auto erased4 = tm4(center, radius);
            auto value4  = wtf::fp::float_cast<tm_type>(erased4);
            REQUIRE(value4.max_order() == 4);
        }
    }
#endif

    SECTION("kind") {
        REQUIRE(defaulted.kind() == UQKind::uncertain);
        REQUIRE(interval.kind() == UQKind::interval);
        REQUIRE(tm4.kind() == UQKind::taylor_model);
    }

    SECTION("operator==") {
        SECTION("Same kind") {
            REQUIRE(defaulted == uncertain);
            REQUIRE(interval == UQFactory(UQKind::interval));
            REQUIRE(affine == UQFactory(UQKind::affine));
            REQUIRE(taffine == UQFactory(UQKind::thresholded_affine));
        }

        SECTION("Different kinds") {
            for(auto lhs : all_kinds) {
                for(auto rhs : all_kinds) {
                    if(lhs == rhs) continue;
                    REQUIRE_FALSE(UQFactory(lhs) == UQFactory(rhs));
                }
            }
        }

        SECTION("TaylorModel additionally compares order") {
            REQUIRE(tm == UQFactory(UQKind::taylor_model, 2));
            REQUIRE_FALSE(tm == tm4);
            REQUIRE(tm4 == UQFactory(UQKind::taylor_model, 4));
        }

        SECTION("Comparison is symmetric") {
            REQUIRE_FALSE(uncertain == tm);
            REQUIRE_FALSE(tm == uncertain);
            REQUIRE_FALSE(tm4 == interval);
            REQUIRE_FALSE(interval == tm4);
        }
    }

    SECTION("operator<") {
        SECTION("Orders by kind") {
            for(std::size_t i = 0; i < all_kinds.size(); ++i) {
                for(std::size_t j = 0; j < all_kinds.size(); ++j) {
                    UQFactory lhs(all_kinds[i]), rhs(all_kinds[j]);
                    REQUIRE((lhs < rhs) == (i < j));
                }
            }
        }

        SECTION("Irreflexive") {
            for(auto kind : all_kinds) {
                UQFactory factory(kind);
                REQUIRE_FALSE(factory < factory);
            }
            REQUIRE_FALSE(tm4 < tm4);
        }

        SECTION("TaylorModel additionally orders by order") {
            REQUIRE(tm < tm4);
            REQUIRE_FALSE(tm4 < tm);
        }
    }
}

TEST_CASE("UQFactoryBase") {
    detail::UncertainFactoryImpl uncertain;
    detail::IntervalFactoryImpl interval;
    detail::TaylorModelFactoryImpl tm2;
    detail::TaylorModelFactoryImpl tm4(4);

    SECTION("kind") {
        REQUIRE(uncertain.kind() == UQKind::uncertain);
        REQUIRE(interval.kind() == UQKind::interval);
        REQUIRE(detail::AffineFactoryImpl{}.kind() == UQKind::affine);
        REQUIRE(detail::ThresholdedAffineFactoryImpl{}.kind() ==
                UQKind::thresholded_affine);
        REQUIRE(tm2.kind() == UQKind::taylor_model);
    }

    SECTION("clone") {
        auto pcopy = interval.clone();
        REQUIRE(pcopy.get() != &interval); // Deep, not shallow
        REQUIRE(pcopy->kind() == UQKind::interval);
        REQUIRE(pcopy->equal(interval));

        auto ptm = tm4.clone();
        REQUIRE(ptm.get() != &tm4);
        REQUIRE(ptm->equal(tm4));
        REQUIRE_FALSE(ptm->equal(tm2));
    }

    SECTION("equal") {
        REQUIRE(uncertain.equal(detail::UncertainFactoryImpl{}));
        REQUIRE_FALSE(uncertain.equal(interval));
        REQUIRE(tm2.equal(detail::TaylorModelFactoryImpl{}));
        REQUIRE_FALSE(tm2.equal(tm4));
        REQUIRE_FALSE(tm2.equal(uncertain));
    }

    SECTION("less") {
        REQUIRE(uncertain.less(interval));
        REQUIRE_FALSE(interval.less(uncertain));
        REQUIRE(tm2.less(tm4));
        REQUIRE_FALSE(tm4.less(tm2));
        REQUIRE_FALSE(tm2.less(tm2));
    }

    SECTION("order") {
        REQUIRE(tm2.order() == 2); // Default
        REQUIRE(tm4.order() == 4);
    }
}
