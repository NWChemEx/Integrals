#include "integrals/integralsmm.hpp"
#include "integrals/property_types.hpp"
#include <catch2/catch.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <testing/testing.hpp>

using namespace testing;

TEST_CASE("DOI") {
    using integral_type = integrals::pt::doi<double>;

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto mol        = testing::get_molecules().at(name);
    auto aos        = testing::get_bases().at(name).at(bs);
    std::vector bases{bs, bs};
    auto tensors = testing::get_ao_data(world).at(name).at(bases);
    auto [X]     = mm.at("DOI").run_as<integral_type>(aos, aos);

    REQUIRE(libchemist::ta_helpers::allclose(X, tensors.at(property::dois)));
}