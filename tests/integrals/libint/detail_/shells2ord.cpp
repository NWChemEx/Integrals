#include "integrals/libint/detail_/shells2ord.hpp"
#include "libint_basis_set_water.hpp"
#include <catch2/catch.hpp>

TEST_CASE("shell2ord") {
    using size_vector_t  = std::vector<std::size_t>;
    using basis_vector_t = std::vector<libint2::BasisSet>;

    auto bset = testing::water_basis_set();

    SECTION("2D") {
        auto ord_pos = integrals::detail_::shells2ord(
          basis_vector_t{bset, bset}, size_vector_t{2, 0});
        REQUIRE(ord_pos == size_vector_t{14, 21, 28});
    }

    SECTION("3D") {
        auto ord_pos = integrals::detail_::shells2ord(
          basis_vector_t{bset, bset, bset}, size_vector_t{2, 0, 0});
        REQUIRE(ord_pos == size_vector_t{98, 147, 196});
    }

    SECTION("4D") {
        auto ord_pos = integrals::detail_::shells2ord(
          basis_vector_t{bset, bset, bset, bset}, size_vector_t{2, 0, 0, 0});
        REQUIRE(ord_pos == size_vector_t{686, 1029, 1372});
    }
}