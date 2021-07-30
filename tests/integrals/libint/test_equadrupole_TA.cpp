// #include "integrals/integralsmm.hpp"
// #include "integrals/property_types.hpp"
// #include <catch2/catch.hpp>
// #include <libchemist/ta_helpers/ta_helpers.hpp>
// #include <property_types/ao_integrals/emultipole.hpp>
// #include <property_types/ao_integrals/overlap.hpp>
// #include <testing/testing.hpp>

// using namespace integrals;
// using namespace testing;

// TEST_CASE("Quadrupole") {
//     using s_type = pt::overlap<double>;
//     using d_type = pt::edipole<double>;
//     using q_type = pt::equadrupole<double>;

//     auto& world = TA::get_default_world();
//     sde::ModuleManager mm;
//     integrals::load_modules(mm);

//     const auto name = molecule::h2o;
//     const auto bs   = basis_set::sto3g;
//     auto mol        = testing::get_molecules().at(name);
//     auto aos        = testing::get_bases().at(name).at(bs);
//     std::vector bases{bs, bs};
//     auto tensors = testing::get_ao_data(world).at(name).at(bases);
//     auto origin  = std::array<double, 3>{0, 0, 0};
//     SECTION("Overlap") {
//         mm.change_input("EQuadrupole", "Origin", origin);
//         auto [S] = mm.at("EQuadrupole").run_as<s_type>(aos, aos);
//         auto X   = TA::retile(tensors.at(property::overlap), S.trange());
//         REQUIRE(libchemist::ta_helpers::allclose(S, X));
//     }
//     SECTION("Dipole") {
//         auto [D] = mm.at("EQuadrupole").run_as<d_type>(origin, aos, aos);
//         auto X   = TA::retile(tensors.at(property::dipole), D.trange());
//         REQUIRE(libchemist::ta_helpers::allclose(D, X));
//     }
//     SECTION("Quadrupole") {
//         auto [Q] = mm.at("EQuadrupole").run_as<q_type>(origin, aos, aos);
//         auto X   = TA::retile(tensors.at(property::quadrupole), Q.trange());
//         REQUIRE(libchemist::ta_helpers::allclose(Q, X));
//     }
// }
