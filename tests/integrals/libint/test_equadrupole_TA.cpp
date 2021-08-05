#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <libchemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>
#include <simde/tensor_representation/tensor_representation.hpp>

using namespace integrals;

TEST_CASE("Quadrupole") {
    using i_op   = simde::type::el_identity;
    using d_op   = simde::type::el_dipole;
    using q_op   = simde::type::el_quadrupole;
    using s_type = simde::AOTensorRepresentation<2, i_op>;
    using d_type = simde::AOTensorRepresentation<2, d_op>;
    using q_type = simde::AOTensorRepresentation<2, q_op>;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = mokup::molecule::h2o;
    const auto bs   = mokup::basis_set::sto3g;
    auto mol        = mokup::get_molecules().at(name);
    auto aos        = mokup::get_bases().at(name).at(bs);
    std::vector bases{bs, bs};
    auto tensors = mokup::get_ao_data(world).at(name).at(bases);
    d_op r;
    q_op r2;

    // SECTION("overlap matrix") {
    //     mm.change_input("EDipole", "Origin", origin);
    //     auto [S] = mm.at("EDipole").run_as<s_type>(aos, aos);
    //     REQUIRE(libchemist::ta_helpers::allclose(S, X));
    // }

    SECTION("dipole matrix") {
        // auto [D]  = mm.at("EQuadrupole").run_as<d_type>(aos, r, aos);
        // auto corr = tensors.at(mokup::property::dipole);
        // REQUIRE(libchemist::tensor::allclose(D, corr));
    }

    SECTION("Quadrupole") {
        // auto [Q]  = mm.at("EQuadrupole").run_as<q_type>(aos, r2, aos);
        // auto corr = tensors.at(mokup::property::quadrupole);
        // REQUIRE(libchemist::tensor::allclose(Q, corr));
    }
}
