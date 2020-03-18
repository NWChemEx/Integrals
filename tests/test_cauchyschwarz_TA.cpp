#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <integrals/nwx_libint/nwx_libint.hpp>
#include <integrals/nwx_libint/cauchy_schwarz.hpp>
#include <integrals/nwx_libint/nwx_libint_factory.hpp>
#include <property_types/ao_integrals/electron_repulsion.hpp>
#include <integrals/nwx_direct/eri_direct_type.hpp>
#include "H2O_STO3G_ERI_SCREENED.hpp"


std::vector<double> two_index {1.0229550772140450,
                               5.1415197380182622,
                               2.4051606675132140,
                               6.3192223244527304,
                               6.3192223244527304};

std::vector<std::vector<double>> four_index {
        {2.1874792352627042, 0.3699640325144978, 0.1564525878919850, 0.0606886000006391, 0.0606886000006391},
        {0.3699640325144978, 0.9039946468459082, 0.4248745604347619, 0.3222057218175323, 0.3222057218175323},
        {0.1564525878919850, 0.4248745604347619, 0.9381679451862797, 0.2740491521872362, 0.2740491521872362},
        {0.0606886000006391, 0.3222057218175323, 0.2740491521872362, 0.8801170058122373, 0.1336470421634229},
        {0.0606886000006391, 0.3222057218175323, 0.2740491521872362, 0.1336470421634229, 0.8801170058122373},
};

const double eps  = 10000 * std::numeric_limits<double>::epsilon();
const double marg = 100 * std::numeric_limits<double>::epsilon();

TEST_CASE("Cauchy-Schwarz") {
    using integral_type = property_types::ERI4CIntegral<double>;
    using direct_type = property_types::ERI4CDirect<double>;

    auto [molecule, bs] = make_molecule();

    auto li_bs = nwx_libint::make_basis(bs);
    std::vector<libint2::BasisSet> sets = {li_bs, li_bs, li_bs};

    auto factory = nwx_libint::LibintFactory();
    factory.max_nprims = libint2::max_nprim(li_bs);
    factory.max_l = libint2::max_l(li_bs);
    factory.thresh = 1.0E-16;
    factory.deriv = 0;

    auto cs_screener = nwx_libint::CauchySchwarz<3, libint2::Operator::coulomb>();
    cs_screener.initialize(sets, factory);

    for (int i = 0; i < 5; ++i) {
        REQUIRE(cs_screener.cs_mat1[0][i] == Approx(two_index[i]).epsilon(eps).margin(marg));
        for (int j = 0; j < 5; ++j) {
            REQUIRE(cs_screener.cs_mat2[i][j] == Approx(four_index[i][j]).epsilon(eps).margin(marg));
        }
    }

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.at("ERI4").change_input("Screening Threshold", 0.005);
    auto [X] = mm.at("ERI4").run_as<integral_type>(bs, bs, bs, bs, std::size_t{0});

    compare_integrals(X, corr);

    mm.at("ERI4Direct").change_input("Screening Threshold", 0.005);
    auto [X_direct] = mm.at("ERI4Direct").run_as<direct_type>(bs, bs, bs, bs, std::size_t{0});

    TensorType X_core(X.world(), X.trange());
    X_core("k, l, m, n") = X_direct("k, l, m, n");

    compare_integrals(X_core, corr);
}

