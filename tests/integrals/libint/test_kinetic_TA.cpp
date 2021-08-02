#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

TEST_CASE("Kinetic") {
    using op_type       = simde::type::el_kinetic;
    using integral_type = simde::AOTensorRepresentation<2, op_type>;
    using size_vector   = std::vector<std::size_t>;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    auto name = mokup::molecule::h2o;
    auto bs   = mokup::basis_set::sto3g;
    auto aos  = mokup::get_bases().at(name).at(bs);
    std::vector bases{bs, bs};
    auto corr = mokup::get_ao_data(world).at(name).at(bases);

    // mm.at("Kinetic").change_input("Tile size", size_vector{1, 2});
    op_type t;
    auto [T] = mm.at("Kinetic").run_as<integral_type>(aos, t, aos);
    REQUIRE(libchemist::tensor::allclose(T, corr.at(mokup::property::kinetic)));
}
