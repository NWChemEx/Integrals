#include "test_error.hpp"
#include <integrals/integrals.hpp>

using eri4_pt = simde::ERI4;
using pt      = integrals::property_types::Uncertainty<eri4_pt>;

using namespace integrals::libint::test;

TEST_CASE("PrimitiveErrorModel") {
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    integrals::set_defaults(mm);

    auto& mod = mm.at("Primitive Error Model");

    simde::type::v_ee_type v_ee{};
    using float_type = double;
    SECTION("H2") {
        auto h2_aos = test::h2_sto3g_basis_set();
        simde::type::aos_squared h2_aos2(h2_aos, h2_aos);
        chemist::braket::BraKet mnls(h2_aos2, v_ee, h2_aos2);
        auto error = mod.run_as<pt>(mnls, 1.0e-10);

        tensorwrapper::shape::Smooth shape({2, 2, 2, 2});
        std::vector<float_type> data(shape.size());
        tensorwrapper::buffer::Contiguous buffer(std::move(data), shape);
        REQUIRE(buffer.approximately_equal(error.buffer(), 1.0e-12));
    }

    SECTION("H2 2") {
        float_type tol                = 1.0e-10;
        auto [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();

        simde::type::aos_squared bra(bra0, bra1);
        simde::type::aos_squared ket(ket0, ket1);
        chemist::braket::BraKet mnls(bra, v_ee, ket);
        auto error = mod.run_as<pt>(mnls, tol);

        auto corr = mm.at("Analytic Error").run_as<pt>(mnls, tol);
    }
}