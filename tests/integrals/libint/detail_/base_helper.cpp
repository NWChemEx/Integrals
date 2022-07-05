#include "integrals/libint/detail_/bases_helper.hpp"
#include "libint_basis_set_water.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>

using namespace mokup;

TEST_CASE("make_libint_basis_set") {
    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto aos        = get_bases(name, bs);

    /// Expected Libint result
    auto libint_corr = testing::water_basis_set();

    /// Check output
    auto libint_bs = integrals::detail_::make_libint_basis_set(aos.basis_set());
    REQUIRE(libint_bs == libint_corr);
}