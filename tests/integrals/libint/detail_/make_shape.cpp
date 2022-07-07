#include "integrals/libint/detail_/make_shape.hpp"
#include "libint_basis_set_water.hpp"
#include <catch2/catch.hpp>

TEST_CASE("make_shape") {
    using extents_t = typename simde::type::tensor::shape_type::extents_type;

    /// Basis set inputs
    auto bset = testing::water_basis_set();
    std::vector<libint2::BasisSet> bsets{bset, bset};

    /// Check output
    SECTION("standard usage") {
        auto shape_ptr = integrals::detail_::make_shape(bsets);
        REQUIRE(shape_ptr->extents() == extents_t{7, 7});
    }

    /// Check leading extent option
    SECTION("with leading extent") {
        std::size_t extra = 3;
        auto shape_ptr    = integrals::detail_::make_shape(bsets, extra);
        REQUIRE(shape_ptr->extents() == extents_t{3, 7, 7});
    }
}