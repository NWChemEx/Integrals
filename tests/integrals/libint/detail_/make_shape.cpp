#include "integrals/libint/detail_/make_shape.hpp"
#include "libint_basis_set_water.hpp"
#include <catch2/catch.hpp>

TEST_CASE("make_shape") {
    using extents_t = typename simde::type::tensor::shape_type::extents_type;

    auto bset = testing::water_basis_set();
    std::vector<libint2::BasisSet> bsets{bset, bset};

    auto shape_ptr = integrals::detail_::make_shape(bsets);

    REQUIRE(shape_ptr->extents() == extents_t{7, 7});
}