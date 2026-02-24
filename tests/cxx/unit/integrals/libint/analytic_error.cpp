#include "../testing/testing.hpp"
#include <integrals/integrals.hpp>

using eri4_pt = simde::ERI4;
using pt      = integrals::property_types::Uncertainty<eri4_pt>;

using namespace integrals::testing;
using tensorwrapper::operations::approximately_equal;
namespace {

template<typename FloatType>
auto corr_error(FloatType diff) {
    tensorwrapper::shape::Smooth shape({1, 1, 1, 1});
    std::vector<FloatType> buffer(shape.size(), diff);
    tensorwrapper::buffer::Contiguous cont(std::move(buffer), shape);
    return simde::type::tensor(std::move(cont), shape);
}

} // namespace

TEST_CASE("AnalyticError") {
    auto mm   = initialize_integrals();
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
        auto corr      = corr_error<float_type>(-0.0000000054867100);
        REQUIRE(approximately_equal(error, corr, 1.0e-12));
    }

    SECTION("H2O/STO-3G (00|34)") {
        auto [bra0, bra1, ket0, ket1] = get_h2o_0034_bases();

        simde::type::aos_squared bra(bra0, bra1);
        simde::type::aos_squared ket(ket0, ket1);
        chemist::braket::BraKet mnls(bra, v_ee, ket);
        float_type tol = 1.0e-10;
        auto error     = mod.run_as<pt>(mnls, tol);
        auto corr      = corr_error<float_type>(-0.0000000000315610);
        REQUIRE(approximately_equal(error, corr, 1.0e-12));
    }
}