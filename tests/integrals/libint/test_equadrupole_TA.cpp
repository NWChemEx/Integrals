#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <chemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>
#include <simde/tensor_representation/tensor_representation.hpp>

using namespace integrals;
using namespace mokup;

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

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto mol        = get_molecule(name);
    auto aos        = get_bases(name, bs);
    std::vector bases{bs, bs};
    // auto tensors = mokup::get_ao_data(world).at(name).at(bases);
    d_op r;
    q_op r2;

    // SECTION("overlap matrix") {
    //     mm.change_input("EDipole", "Origin", origin);
    //     auto [S] = mm.at("EDipole").run_as<s_type>(aos, aos);
    //     REQUIRE(chemist::ta_helpers::allclose(S, X));
    // }

    SECTION("dipole matrix") {
        // auto [D]  = mm.at("EQuadrupole").run_as<d_type>(aos, r, aos);
        // auto corr = tensors.at(mokup::property::dipole);
        // REQUIRE(chemist::tensor::allclose(D, corr));
    }

    SECTION("Quadrupole") {
        // auto [Q]  = mm.at("EQuadrupole").run_as<q_type>(aos, r2, aos);
        // auto corr = tensors.at(mokup::property::quadrupole);
        // REQUIRE(chemist::tensor::allclose(Q, corr));
    }
}
