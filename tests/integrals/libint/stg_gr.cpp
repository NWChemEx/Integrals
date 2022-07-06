#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>
#include <tensorwrapper/tensor/allclose.hpp>

using namespace mokup;

TEST_CASE("STG 4 Center GR") {
    using op_type       = simde::type::el_el_yukawa;
    using integral_type = simde::AOTensorRepresentation<4, op_type>;
    const auto key      = "Yukawa4";

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = mokup::molecule::h2;
    const auto prop = property::stg_gr;

    chemist::Electron e;
    op_type gr(chemist::operators::STG(), e, e);

    for(const auto& bs : {basis_set::sto3g, basis_set::ccpvdz}) {
        std::vector<basis_set> bs_key(4, bs);
        SECTION(as_string(name, bs)) {
            auto aos    = get_bases(name, bs);
            auto X_corr = get_ao_data(name, bs_key, prop);
            auto [X] = mm.at(key).run_as<integral_type>(aos, aos, gr, aos,
            aos); REQUIRE(tensorwrapper::tensor::allclose(X, X_corr));
        }
    }
}
