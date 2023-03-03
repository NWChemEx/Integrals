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

#include "integrals/ao_integrals/detail_/make_shape.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>

using namespace mokup;

using extents_t = typename simde::type::tensor::shape_type::extents_type;
using integrals::ao_integrals::detail_::make_shape;

TEST_CASE("make_shape") {
    /// Basis set inputs
    const auto name  = molecule::h2o;
    const auto bs    = basis_set::sto3g;
    auto aos         = get_bases(name, bs);
    const auto& bset = aos.basis_set();
    std::vector bsets{bset, bset};

    /// Check output
    SECTION("standard usage") {
        auto shape_ptr = make_shape(bsets);
        REQUIRE(shape_ptr->extents() == extents_t{7, 7});
    }

    /// Check leading extent option
    SECTION("with leading extent") {
        std::size_t extra = 3;
        auto shape_ptr    = make_shape(bsets, extra);
        REQUIRE(shape_ptr->extents() == extents_t{3, 7, 7});
    }
}