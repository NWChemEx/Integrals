#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

TEST_CASE("DOI") {
    using op_type       = simde::type::el_el_delta;
    using integral_type = simde::EDOI;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = mokup::molecule::h2o;
    const auto bs   = mokup::basis_set::sto3g;
    auto mol        = mokup::get_molecules().at(name);
    auto aos        = mokup::get_bases().at(name).at(bs);
    std::vector bases{bs, bs};
    auto tensors = mokup::get_ao_data(world).at(name).at(bases);
    op_type d;
    auto [X]  = mm.at("DOI").run_as<integral_type>(aos, d, aos);
    auto corr = tensors.at(mokup::property::dois);
    REQUIRE(libchemist::tensor::allclose(X, corr));
}
