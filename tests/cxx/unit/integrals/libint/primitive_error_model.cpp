// #include "../testing/testing.hpp"
// #include <integrals/integrals.hpp>

// using eri4_pt = simde::ERI4;
// using pt      = integrals::property_types::Uncertainty<eri4_pt>;

// using namespace integrals::testing;

// TEST_CASE("PrimitiveErrorModel") {
//     auto mm = initialize_integrals();

//     auto& mod = mm.at("Primitive Error Model");

//     simde::type::v_ee_type v_ee{};
//     using float_type = double;
//     SECTION("H2") {
//         auto h2_aos = test::h2_sto3g_basis_set();
//         simde::type::aos_squared h2_aos2(h2_aos, h2_aos);
//         chemist::braket::BraKet mnls(h2_aos2, v_ee, h2_aos2);
//         auto error = mod.run_as<pt>(mnls, 1.0e-10);

//         tensorwrapper::shape::Smooth shape({2, 2, 2, 2});
//         std::vector<float_type> data(shape.size());
//         tensorwrapper::buffer::Contiguous buffer(std::move(data), shape);
//         REQUIRE(buffer.approximately_equal(error.buffer(), 1.0e-12));
//     }

//     SECTION("H2 2") {
//         float_type tol                = 1.0e-10;
//         auto [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();
//         auto value = manual_contract_shell(bra0, bra1, ket0, ket1, mm);
//         std::cout << "Screened value: " << value << std::endl;

//         simde::type::aos_squared bra(bra0, bra1);
//         simde::type::aos_squared ket(ket0, ket1);
//         chemist::braket::BraKet mnls(bra, v_ee, ket);

//         auto error = mod.run_as<pt>(mnls, tol);
//         auto corr  = mm.at("Analytic Error").run_as<pt>(mnls, tol);

//         auto eri4_mod = mm.at("ERI4").unlocked_copy();
//         eri4_mod.change_input("Threshold", float_type(1E-16));
//         auto i_value = eri4_mod.run_as<eri4_pt>(mnls);
//         std::cout << "Corr Screened Value: " << i_value << std::endl;

//         // Our value:      0.0002263495591894
//         // Libint's value: 0.0002263440759384
//         // No screening:   0.0002263495626484
//     }
// }