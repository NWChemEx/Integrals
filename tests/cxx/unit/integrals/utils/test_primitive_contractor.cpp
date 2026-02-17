#include "../testing/testing.hpp"
using eri4_pt = simde::ERI4;

using namespace integrals::testing;
using simde::type::tensor;
using tensorwrapper::operations::approximately_equal;

namespace {

template<typename ModuleType, typename BraKetType>
void test_eri(ModuleType&& mod, ModuleType&& corr_mod,
              ModuleType&& corr_screen_mod, BraKetType&& mnls) {
    auto inputs = mod.inputs();

    eri4_pt::wrap_inputs(inputs, mnls);
    auto rv = mod.run(inputs);

    // Unpack the values we computed
    auto screen   = rv.at("Tensor Representation").template value<tensor>();
    auto eri      = rv.at("Corr ERI4").template value<tensor>();
    auto error    = rv.at("Actual Error").template value<tensor>();
    auto bb_error = rv.at("Black-Box Error Estimate").template value<tensor>();

    // Compute the correct values using the Libint-based modules
    auto corr_screen = corr_screen_mod.template run_as<eri4_pt>(mnls);
    auto corr_eri    = corr_mod.template run_as<eri4_pt>(mnls);

    tensor corr_error;
    corr_error("i,j,k,l") = corr_screen("i,j,k,l") - corr_eri("i,j,k,l");

    REQUIRE(approximately_equal(corr_eri, eri, 1.0E-15));
    // REQUIRE(approximately_equal(corr_screen, screen, 1.0E-10));

    std::cout << "Correct error vs our error: " << corr_error << " " << error
              << std::endl;
    std::cout << "Error vs BB error: " << error << " " << bb_error << std::endl;
    std::cout << "Correct screened vs our screened: " << corr_screen << " "
              << screen << std::endl;
    // std::cout << "Correct ERI vs our ERI: " << corr_eri << " " << eri
    //           << std::endl;
}
} // namespace

TEST_CASE("PrimitiveContractor") {
    auto mm = initialize_integrals();

    double tol = 1.0E-10;

    mm.copy_module("ERI4", "Benchmark ERI4");
    mm.change_submod("Prototype: ERI4", "Primitive ERI4", "Benchmark ERI4");
    mm.change_input("ERI4", "threshold", tol);
    mm.change_input("Prototype: ERI4", "threshold", tol);

    // Module that we wrote to screen the integrals and estimate the error
    auto mod = mm.at("Prototype: ERI4");

    // Module that uses Libint to compute the correct value
    auto corr_mod = mm.at("Benchmark ERI4");

    // Module that uses Libint to compute the screened value
    auto corr_screen_mod = mm.at("ERI4");

    // Make the bases and the operator we are going to test on
    auto [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();
    simde::type::aos_squared bra(bra0, bra1);
    simde::type::aos_squared ket(ket0, ket1);
    simde::type::v_ee_type v_ee;

    // SECTION("one shell") {
    //     simde::type::aos_squared bra(bra0, bra0);
    //     chemist::braket::BraKet mnls(bra, v_ee, bra);
    //     auto I      = mod.run_as<eri4_pt>(mnls);
    //     auto corr_I = corr_mod.run_as<eri4_pt>(mnls);
    //     // REQUIRE(approximately_equal(I, corr_I, 1.0E-12));
    //     std::cout << I << std::endl;
    //     std::cout << corr_I << std::endl;
    // }

    SECTION("(ss|ss)") {
        chemist::braket::BraKet mnls(bra, v_ee, ket);
        test_eri(mod, corr_mod, corr_screen_mod, mnls);
    }

    SECTION("(ps|ss)") {
        bra.first.ao_basis_set().shell(0).l() = 1;
        chemist::braket::BraKet mnls(bra, v_ee, ket);
        test_eri(mod, corr_mod, corr_screen_mod, mnls);
    }
}