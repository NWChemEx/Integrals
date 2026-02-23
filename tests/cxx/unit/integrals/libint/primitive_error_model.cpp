#include "../testing/testing.hpp"
#include <integrals/integrals.hpp>

using eri4_pt = simde::ERI4;
using pt      = integrals::property_types::Uncertainty<eri4_pt>;

using namespace integrals::testing;

TEST_CASE("PrimitiveErrorModel") {
    auto mm = initialize_integrals();

    auto& mod                     = mm.at("Primitive Error Model");
    auto& anal_error              = mm.at("Analytic Error");
    auto [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();
    simde::type::aos_squared bra(bra0, bra1);
    simde::type::v_ee_type v_ee{};
    simde::type::aos_squared ket(ket0, ket1);
    chemist::braket::BraKet mnls(bra, v_ee, ket);
    double tol = 1E-10;
    auto error = mod.run_as<pt>(mnls, tol);
    std::cout << error << std::endl;
    std::cout << anal_error.run_as<pt>(mnls, tol) << std::endl;
}