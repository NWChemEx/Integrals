#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <integrals/nwx_libint/nwx_libint.hpp>
#include <integrals/nwx_libint/cauchy_schwarz_screener.hpp>
#include <property_types/ao_integrals/electron_repulsion.hpp>
#include <property_types/cauchy_schwarz_approximation.hpp>
#include <integrals/nwx_direct/eri_direct_type.hpp>
#include "H2O_STO3G_ERI_SCREENED.hpp"

std::vector<std::vector<double>> eri_check {
        {2.1874792352627042, 0.3699640325144978, 0.1564525878919850, 0.0606886000006391, 0.0606886000006391},
        {0.3699640325144978, 0.9039946468459082, 0.4248745604347619, 0.3222057218175323, 0.3222057218175323},
        {0.1564525878919850, 0.4248745604347619, 0.9381679451862797, 0.2740491521872362, 0.2740491521872362},
        {0.0606886000006391, 0.3222057218175323, 0.2740491521872362, 0.8801170058122373, 0.1336470421634229},
        {0.0606886000006391, 0.3222057218175323, 0.2740491521872362, 0.1336470421634229, 0.8801170058122373},
};

std::vector<std::vector<double>> stg_check {
        {0.8712382481019230, 0.1833574527284384, 0.0499799084050583, 0.0298548303674855, 0.0298548303674855},
        {0.1833574527284384, 0.5049372580585069, 0.2429078076543889, 0.1760787877695699, 0.1760787877695699},
        {0.0499799084050583, 0.2429078076543889, 0.5178471413265817, 0.1563617001398491, 0.1563617001398491},
        {0.0298548303674855, 0.1760787877695699, 0.1563617001398491, 0.4834953496227799, 0.0689823476694450},
        {0.0298548303674855, 0.1760787877695699, 0.1563617001398491, 0.0689823476694450, 0.4834953496227799},
};

std::vector<std::vector<double>> yuk_check {
        {1.9779965679341631, 0.3052904918032274, 0.1515386376453260, 0.0502970255598737, 0.0502970255598737},
        {0.3052904918032274, 0.5375014882330088, 0.3623472474047580, 0.1781992287690466, 0.1781992287690466},
        {0.1515386376453260, 0.3623472474047580, 0.5862346971968084, 0.1769470889613180, 0.1769470889613180},
        {0.0502970255598737, 0.1781992287690466, 0.1769470889613180, 0.5201331598096510, 0.0647323127080033},
        {0.0502970255598737, 0.1781992287690466, 0.1769470889613180, 0.0647323127080033, 0.5201331598096510},
};

const double eps  = 10000 * std::numeric_limits<double>::epsilon();
const double marg = 100 * std::numeric_limits<double>::epsilon();

TEST_CASE("Cauchy-Schwarz") {
    using integral_type = property_types::ERI4CIntegral<double>;
    using direct_type = property_types::ERI4CDirect<double>;
    using cs_approx = property_types::CauchySchwarzApprox<double>;

    auto [molecule, bs] = make_molecule();

    sde::ModuleManager mm;
    integrals::load_modules(mm);

    mm.at("CauchySchwarzSTG").change_input("STG Exponent", 1.0);
    mm.at("CauchySchwarzYukawa").change_input("STG Exponent", 1.0);

    // Check Cauchy-Schwarz Module
    auto [cs_mat_eri] = mm.at("CauchySchwarzERI").run_as<cs_approx>(bs, bs, std::size_t{0});
    auto [cs_mat_stg] = mm.at("CauchySchwarzSTG").run_as<cs_approx>(bs, bs, std::size_t{0});
    auto [cs_mat_yuk] = mm.at("CauchySchwarzYukawa").run_as<cs_approx>(bs, bs, std::size_t{0});
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            REQUIRE(cs_mat_eri[i][j] == Approx(eri_check[i][j]).epsilon(eps).margin(marg));
            REQUIRE(cs_mat_stg[i][j] == Approx(stg_check[i][j]).epsilon(eps).margin(marg));
            REQUIRE(cs_mat_yuk[i][j] == Approx(yuk_check[i][j]).epsilon(eps).margin(marg));
        }
    }

    // Check Screener Class
    auto li_bs = nwx_libint::make_basis(bs);
    std::vector<libint2::BasisSet> sets = {li_bs, li_bs, li_bs, li_bs};
    TA::TiledRange1 rng{0, 5, 6, 7};

    auto cs_screener_2c = nwx_libint::CauchySchwarzScreener<2>();
    auto cs_screener_3c = nwx_libint::CauchySchwarzScreener<3>();
    auto cs_screener_4c = nwx_libint::CauchySchwarzScreener<4>();
    cs_screener_3c.cs_mat2 = cs_mat_eri;
    cs_screener_4c.cs_mat1 = cs_mat_eri;
    cs_screener_4c.cs_mat2 = cs_mat_eri;
    REQUIRE(!cs_screener_2c.shellset({0, 0}, 0.0));
    REQUIRE(!cs_screener_3c.shellset({0, 0, 0}, 0.0));
    REQUIRE(!cs_screener_4c.shellset({0, 0, 0, 0}, 0.0));
    REQUIRE(cs_screener_3c.shellset({0, 0, 0}, 5.0));
    REQUIRE(cs_screener_4c.shellset({0, 0, 0, 0}, 5.0));

    // Check Screened Integrals
    mm.at("ERI4").change_input("Screening Threshold", 0.005);
    auto [X] = mm.at("ERI4").run_as<integral_type>(bs, bs, bs, bs, std::size_t{0});
    REQUIRE(libchemist::allclose(X, TensorType(X.world(), X.trange(), corr)));

    mm.at("ERI4Direct").change_input("Screening Threshold", 0.005);
    auto [X_direct] = mm.at("ERI4Direct").run_as<direct_type>(bs, bs, bs, bs, std::size_t{0});
    TensorType X_core(X.world(), X.trange());
    X_core("k, l, m, n") = X_direct("k, l, m, n");
    REQUIRE(libchemist::allclose(X_core, TensorType(X_core.world(), X_core.trange(), corr)));
}

