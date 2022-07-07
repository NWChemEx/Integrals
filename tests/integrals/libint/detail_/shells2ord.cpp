#include "integrals/libint/detail_/shells2ord.hpp"
#include "libint_basis_set_water.hpp"
#include <catch2/catch.hpp>

TEST_CASE("shell2ord") {
    using size_vector_t  = std::vector<std::size_t>;
    using basis_vector_t = std::vector<libint2::BasisSet>;

    /// Common basis set
    auto bset = testing::water_basis_set();

    /// Check different dimensionalities
    SECTION("2D") {
        basis_vector_t sets{bset, bset};
        size_vector_t curr{2, 0}, lo{0, 0}, up{4, 4};
        auto ord_pos = integrals::detail_::shells2ord(sets, curr, lo, up);
        REQUIRE(ord_pos == size_vector_t{14, 21, 28});
    }

    SECTION("3D") {
        basis_vector_t sets{bset, bset, bset};
        size_vector_t curr{2, 0, 0}, lo{0, 0, 0}, up{4, 4, 4};
        auto ord_pos = integrals::detail_::shells2ord(sets, curr, lo, up);
        REQUIRE(ord_pos == size_vector_t{98, 147, 196});
    }

    SECTION("4D") {
        basis_vector_t sets{bset, bset, bset, bset};
        size_vector_t curr{2, 0, 0, 0}, lo{0, 0, 0, 0}, up{4, 4, 4, 4};
        auto ord_pos = integrals::detail_::shells2ord(sets, curr, lo, up);
        REQUIRE(ord_pos == size_vector_t{686, 1029, 1372});
    }

    SECTION("specific tile") {
        basis_vector_t sets{bset, bset};
        size_vector_t lo{2, 0}, up{3, 0};
        SECTION("lower") {
            auto ord_pos = integrals::detail_::shells2ord(sets, lo, lo, up);
            REQUIRE(ord_pos == size_vector_t{0, 1, 2});
        }
        SECTION("upper") {
            auto ord_pos = integrals::detail_::shells2ord(sets, up, lo, up);
            REQUIRE(ord_pos == size_vector_t{3});
        }
    }
}