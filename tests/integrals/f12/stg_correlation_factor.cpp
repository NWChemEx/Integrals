#include "../../test_common_TA.hpp"
#include "integrals/integralsmm.hpp"
#include "integrals/property_types.hpp"
#include "nwx_testing/nwx_testing.hpp"
#include <libchemist/libchemist.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <libint2.hpp>

TEST_CASE("STG 2 Center Correlation Factor") {
    using integral_type = integrals::pt::correlation_factor_2c<double>;
    using overlap_type  = integrals::pt::overlap<double>;
    const auto key      = "STG 2 Center Correlation Factor";

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto mols  = nwx_testing::get_mols();
    auto bases = nwx_testing::get_bases();
    auto data  = nwx_testing::get_data(world);

    SECTION("H2") {
        const auto name = "h2";
        auto mol        = mols.at(name);
        for(auto bs : {"sto-3g", "cc-pvdz"}) {
            SECTION(bs) {
                auto aos     = bases.at(name).at(bs);
                auto tensors = data.at(name).at(bs);
                auto X_corr  = tensors.at("2 center correlation factor");
                auto [X]     = mm.at(key).run_as<integral_type>(aos, aos);
                REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
            }
        }
    }
}

TEST_CASE("STG 4 Center Correlation Factor") {
    using integral_type = integrals::pt::correlation_factor_4c<double>;
    using overlap_type  = integrals::pt::overlap<double>;
    const auto key      = "STG 4 Center Correlation Factor";

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);

    SECTION("H2") {
        auto mol      = nwx_testing::get_h2_mol();
        const auto bs = "cc-pvdz";
        auto aos      = nwx_testing::get_h2_bases().at(bs);
        auto tensors  = nwx_testing::h2_data(world).at(bs);
        auto X_corr   = tensors.at("4 center correlation factor");
        auto [X]      = mm.at(key).run_as<integral_type>(aos, aos, aos, aos);
        nwx_testing::print_large_diffs(X, X_corr);
        REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
    }
    SECTION("H2O") {
        auto mol      = nwx_testing::get_h2o_mol();
        const auto bs = "sto-3g";
        auto aos      = libchemist::apply_basis(bs, mol);
        auto tensors  = nwx_testing::h2o_data(world).at(bs);
        auto X_corr   = tensors.at("correlation factor");

        auto [X] = mm.at(key).run_as<integral_type>(aos, aos, aos, aos);
        // TA::TSpArrayD diff;
        // diff("i,j,k,l") = X("i,j,k,l") - X_corr("i,j,k,l");
        // auto t          = diff.find(0).get();
        // for(auto i = 0; i < 7; ++i)
        //     for(auto j = 0; j < 7; ++j)
        //         for(auto k = 0; k < 7; ++k)
        //             for(auto l = 0; l < 7; ++l) {
        //                 auto ijkl{i, j, k, l};
        //                 auto t_val = t(ijkl);
        //                 if(std::fabs(t_val) > 0.0001) {
        //                     std::cout << i << " " << j << " " << k << " " <<
        //                     l
        //                               << " ";
        //                     std::cout << t_val << std::endl;
        //                 }
        //             }
        // REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
    }
}