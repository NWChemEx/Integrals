#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <chemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>
#include <simde/tensor_representation/tensor_representation.hpp>

using namespace integrals;
using namespace mokup;

TEST_CASE("Octupole") {
    using i_op   = simde::type::el_identity;
    using d_op   = simde::type::el_dipole;
    using q_op   = simde::type::el_quadrupole;
    using o_op   = simde::type::el_octupole;
    using s_type = simde::AOTensorRepresentation<2, i_op>;
    using d_type = simde::AOTensorRepresentation<2, d_op>;
    using q_type = simde::AOTensorRepresentation<2, q_op>;
    using o_type = simde::AOTensorRepresentation<2, o_op>;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto mol        = get_molecule(name);
    auto aos        = get_bases(name, bs);
    std::vector bases{bs, bs};
    // auto corr = get_ao_data(name, bases, property::octopole, world);
    d_op r;
    q_op r2;
    o_op r3;

    // SECTION("overlap matrix") {
    //     mm.change_input("EDipole", "Origin", origin);
    //     auto [S] = mm.at("EDipole").run_as<s_type>(aos, aos);
    //     REQUIRE(chemist::ta_helpers::allclose(S, X));
    // }

    SECTION("dipole matrix") {
        // auto [D]  = mm.at("EOctupole").run_as<d_type>(aos, r, aos);
        // auto corr = tensors.at(mokup::property::dipole);
        // REQUIRE(chemist::tensor::allclose(D, corr));
    }

    SECTION("Quadrupole") {
        // auto [Q]  = mm.at("EOctupole").run_as<q_type>(aos, r2, aos);
        // auto corr = tensors.at(mokup::property::quadrupole);
        // REQUIRE(chemist::tensor::allclose(Q, corr));
    }

    SECTION("Octupole") {
        // auto [O]  = mm.at("EOctupole").run_as<o_type>(aos, r3, aos);
        // auto corr = tensors.at(mokup::property::octopole);
        // REQUIRE(chemist::tensor::allclose(O, corr));
    }
}
