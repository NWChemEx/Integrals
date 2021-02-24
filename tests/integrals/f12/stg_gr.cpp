#include "../../test_common_TA.hpp"
#include "integrals/integralsmm.hpp"
#include "integrals/property_types.hpp"
#include "testing/testing.hpp"
#include <libchemist/libchemist.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>

// TEST_CASE("STG 2 Center Correlation Factor") {
//     using integral_type = integrals::pt::corre_2c<double>;
//     using overlap_type  = integrals::pt::overlap<double>;
//     const auto key      = "STG 2 Center Correlation Factor";

//     auto& world = TA::get_default_world();
//     sde::ModuleManager mm;
//     integrals::load_modules(mm);

//     SECTION("H2") {
//         auto mol      = nwx_testing::get_h2_mol();
//         const auto bs = "cc-pvdz";
//         auto aos      = libchemist::apply_basis(bs, mol);
//         auto tensors  = nwx_testing::h2_data(world).at(bs);
//         auto X_corr   = tensors.at("2 center correlation factor");

//         auto [X] = mm.at(key).run_as<integral_type>(aos, aos);
//         // REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
//     }
// }

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
        auto mol        = mols.at(mol_name);
        const auto bs   = "cc-pvdz";
        auto aos        = libchemist::apply_basis(bs, mol);
        auto tensors    = data.at(mol_name).at(bs);
        auto X_corr     = tensors.at(prop_name);
        auto [X]        = mm.at(key).run_as<integral_type>(aos, aos, aos, aos);
        double max_diff = (X("i,j,k,l") - X_corr("i,j,k,l")).abs_max();
        std::cout << max_diff << std::endl;
        REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
    }

    // SECTION("H2O") {
    //     auto mol      = nwx_testing::get_h2o_mol();
    //     const auto bs = "sto-3g";
    //     auto aos      = libchemist::apply_basis(bs, mol);
    //     auto tensors  = nwx_testing::h2o_data(world).at(bs);
    //     auto X_corr   = tensors.at("correlation factor");

    //     auto [X] = mm.at(key).run_as<integral_type>(aos, aos, aos, aos);
    //     // TA::TSpArrayD diff;
    //     // diff("i,j,k,l") = X("i,j,k,l") - X_corr("i,j,k,l");
    //     // auto t          = diff.find(0).get();
    //     // for(auto i = 0; i < 7; ++i)
    //     //     for(auto j = 0; j < 7; ++j)
    //     //         for(auto k = 0; k < 7; ++k)
    //     //             for(auto l = 0; l < 7; ++l) {
    //     //                 auto ijkl{i, j, k, l};
    //     //                 auto t_val = t(ijkl);
    //     //                 if(std::fabs(t_val) > 0.0001) {
    //     //                     std::cout << i << " " << j << " " << k << " "
    //     <<
    //     //                     l
    //     //                               << " ";
    //     //                     std::cout << t_val << std::endl;
    //     //                 }
    //     //             }
    //     // REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
    // }
}