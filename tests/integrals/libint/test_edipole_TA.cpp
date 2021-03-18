#include "integrals/integralsmm.hpp"
#include "integrals/property_types.hpp"
#include <catch2/catch.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <testing/testing.hpp>

using namespace integrals;

TEST_CASE("Dipole") {
    using s_type = pt::overlap<double>;
    using d_type = pt::edipole<double>;

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = "h2o";
    const auto bs   = "sto-3g";
    auto mol        = testing::get_molecules().at(name);
    auto aos        = testing::get_bases().at(name).at(bs);
    auto tensors    = testing::get_data(world).at(name).at(bs);
    std::array<double, 3> origin{0, 0, 0};

    SECTION("overlap matrix") {
        mm.change_input("EDipole", "Origin", origin);
        auto [S] = mm.at("EDipole").run_as<s_type>(aos, aos);
        auto X   = TA::retile(tensors.at("Overlap"), S.trange());
        REQUIRE(libchemist::ta_helpers::allclose(S, X));
    }

    SECTION("dipole matrix") {
        auto [D] = mm.at("EDipole").run_as<d_type>(origin, aos, aos);
        auto X   = TA::retile(tensors.at("Dipole"), D.trange());
        REQUIRE(libchemist::ta_helpers::allclose(D, X));
    }
}
