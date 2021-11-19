#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <chemist/tensor/allclose.hpp>
#include <mokup/mokup.hpp>

using namespace integrals;
using namespace mokup;

TEST_CASE("Dipole") {
    using i_op   = simde::type::el_identity;
    using d_op   = simde::type::el_dipole;
    using s_type = simde::AOTensorRepresentation<2, i_op>;
    using d_type = simde::AOTensorRepresentation<2, d_op>;

    auto& world = TA::get_default_world();
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto aos        = get_bases(name, bs);
    std::vector bases{bs, bs};
    auto corr = get_ao_data(name, bases, property::dipole, world);
    d_op r;

    // SECTION("overlap matrix") {
    //     mm.change_input("EDipole", "Origin", origin);
    //     auto [S] = mm.at("EDipole").run_as<s_type>(aos, aos);
    //     REQUIRE(chemist::ta_helpers::allclose(S, X));
    // }

    SECTION("dipole matrix") {
        auto [D] = mm.at("EDipole").run_as<d_type>(aos, r, aos);
        REQUIRE(chemist::tensor::allclose(D, corr));
    }
}
