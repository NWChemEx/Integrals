#include "../../test_common_TA.hpp"
#include "integrals/integrals.hpp"
#include <testing/testing.hpp>

TEST_CASE("ERI4C CS") {
    using integral_type = integrals::pt::eri4c<double>;

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.change_input("ERI4 CS", "Screening Threshold", 0.005);

    const auto name = "h2o";
    const auto bs   = "sto-3g";
    auto mol        = testing::get_molecules().at(name);
    auto aos        = testing::get_bases().at(name).at(bs);
    auto tensors    = testing::get_data(world).at(name).at(bs);
    auto corr_S     = tensors.at("Screened ERIs");
    auto [X] = mm.at("ERI4 CS").run_as<integral_type>(aos, aos, aos, aos);

    REQUIRE(libchemist::ta_helpers::allclose(X, corr_S));
}