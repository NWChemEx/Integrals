#pragma once
#include "integrals/libint/detail_/make_engine.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>

namespace testing {

/// Checks the common parts of the test engines
template<typename OpType>
void test_engine_standard(libint2::Engine& engine) {
    REQUIRE(engine.oper() == integrals::op_v<OpType>);
    REQUIRE(engine.max_nprim() == 3);
    REQUIRE(engine.max_l() == 1);
    REQUIRE(engine.precision() == 1.0E-16);
}

} // namespace testing