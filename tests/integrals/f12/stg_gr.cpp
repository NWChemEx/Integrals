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
    auto mols            = testing::get_molecules();
    auto data            = testing::get_data(world);
    const auto mol_name  = "h2";
    const auto prop_name = "STG 4C GR";

    SECTION(mol_name) {
        auto mol = mols.at(mol_name);
        for(auto bs : {"sto-3g", "cc-pvdz"}) {
            auto aos     = libchemist::apply_basis(bs, mol);
            auto tensors = data.at(mol_name).at(bs);
            auto X_corr  = tensors.at(prop_name);
            auto [X]     = mm.at(key).run_as<integral_type>(aos, aos, aos, aos);
            // double max_diff = (X("i,j,k,l") - X_corr("i,j,k,l")).abs_max();
            // std::cout << max_diff << std::endl;
            // testing::print_large_deviations(X, X_corr);
            // std::cout << std::fixed << std::setprecision(16) << X <<
            // std::endl;
            // REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
        }
    }
}