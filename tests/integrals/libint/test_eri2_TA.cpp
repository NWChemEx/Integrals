#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

TEST_CASE("ERI2C") {
    using op            = simde::type::el_el_coulomb;
    using integral_type = simde::AOTensorRepresentation<2, op>;
    using size_vector   = std::vector<std::size_t>;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    auto name = mokup::molecule::h2o;
    auto bs   = mokup::basis_set::sto3g;
    auto mol  = mokup::get_molecules().at(name);
    auto aos  = mokup::get_bases().at(name).at(bs);
    std::vector bases{bs, bs};
    auto corr = mokup::get_ao_data(world).at(name).at(bases);

    // mm.at("ERI2").change_input("Tile size", size_vector{1});
    simde::type::el_el_coulomb r12;
    auto [X] = mm.at("ERI2").run_as<integral_type>(aos, r12, aos);
    REQUIRE(libchemist::tensor::allclose(X, corr.at(mokup::property::eris)));
}
