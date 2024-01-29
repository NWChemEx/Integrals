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

#pragma once
#include "integrals/libint/detail_/make_engine.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>

namespace testing {

/// Checks the common parts of the test engines
template<typename OpType>
void test_engine_standard(libint2::Engine& engine) {
    REQUIRE(engine.oper() == integrals::libint::op_v<OpType>);
    REQUIRE(engine.max_nprim() == 3);
    REQUIRE(engine.max_l() == 1);
    REQUIRE(engine.precision() == 1.0E-16);
}

} // namespace testing