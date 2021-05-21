#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <testing/testing.hpp>

using namespace testing;

TEST_CASE("Yukawa4C") {
    using integral_type = integrals::pt::yukawa4c<double>;
    using size_vector   = integrals::type::size_vector;
    const auto key1     = "Yukawa4";
    const auto key2     = "Yukawa4 CS";

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto mol        = testing::get_molecules().at(name);
    auto aos        = testing::get_bases().at(name).at(bs);
    std::vector bases{bs, bs, bs, bs};
    auto tensors      = testing::get_ao_data(world).at(name).at(bases);
    auto stg_exponent = 1.0;

    auto [X] = mm.run_as<integral_type>(key1, stg_exponent, aos, aos, aos, aos);
    auto corr_R = TA::retile(tensors.at(property::yukawa), X.trange());
    REQUIRE(libchemist::ta_helpers::allclose(X, corr_R));

    mm.change_input(key2, "Screening Threshold", 0.000001);
    auto [X2] =
      mm.run_as<integral_type>(key2, stg_exponent, aos, aos, aos, aos);
    REQUIRE(libchemist::ta_helpers::allclose(X2, corr_R));
}
