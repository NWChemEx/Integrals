#include "../../test_common_TA.hpp"
#include "integrals/integralsmm.hpp"
#include "integrals/property_types.hpp"
#include "testing/testing.hpp"
#include <libchemist/libchemist.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>

TEST_CASE("STG 4 Center GR") {
    using integral_type = integrals::pt::gr4c<double>;
    const auto key      = "STG 4 Center GR";

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto data            = testing::get_data(world);
    auto bases           = testing::get_bases();
    const auto mol_name  = "h2";
    const auto prop_name = "STG 4C GR";

    SECTION(mol_name) {
        for(auto bs : {"sto-3g", "cc-pvdz"}) {
            auto tensors = data.at(mol_name).at(bs);
            auto aos     = bases.at(mol_name).at(bs);
            auto X_corr  = tensors.at(prop_name);
            auto [X]     = mm.at(key).run_as<integral_type>(aos, aos, aos, aos);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
        }
    }
}