#include "test_error.hpp"
#include <integrals/integrals.hpp>

using eri4_pt = simde::ERI4;
using pt      = integrals::property_types::Uncertainty<eri4_pt>;

using namespace integrals::libint::test;

TEST_CASE("AnalyticError") {
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    integrals::set_defaults(mm);

    auto& mod = mm.at("Analytic Error");

    simde::type::v_ee_type v_ee{};
    using float_type = double;

    SECTION("H2 Dimer/STO-3G (03|12)") {
        auto [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();

        simde::type::aos_squared bra(bra0, bra1);
        simde::type::aos_squared ket(ket0, ket1);
        chemist::braket::BraKet mnls(bra, v_ee, ket);
        float_type tol = 1.0e-10;
        auto error     = mod.run_as<pt>(mnls, tol);
        tensorwrapper::shape::Smooth shape({1, 1, 1, 1});
        auto benchmark = test::make_tw_buffer(0.0002263495626484, shape);
        auto screened  = test::make_tw_buffer(0.0002263440759384, shape);
        tensorwrapper::buffer::Contiguous corr(benchmark);
        corr("m,n,l,s") = screened("m,n,l,s") - benchmark("m,n,l,s");
        REQUIRE(error.buffer().approximately_equal(corr, 1.0e-12));
    }
}