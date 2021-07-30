// #include "integrals/integrals.hpp"
// #include <catch2/catch.hpp>
// #include <libchemist/ta_helpers/ta_helpers.hpp>
// #include <testing/testing.hpp>

// using namespace testing;

// TEST_CASE("Yukawa2C") {
//     using integral_type = integrals::pt::yukawa2c<double>;
//     using size_vector   = integrals::type::size_vector;

//     auto& world = TA::get_default_world();
//     sde::ModuleManager mm;
//     integrals::load_modules(mm);
//     auto name = molecule::h2o;
//     auto bs   = basis_set::sto3g;
//     auto aos  = testing::get_bases().at(name).at(bs);
//     std::vector bases{bs, bs};
//     auto corr         = testing::get_ao_data(world).at(name).at(bases);
//     auto stg_exponent = 1.0;
//     mm.at("Yukawa2").change_input("Tile size", size_vector{6});
//     auto [X] = mm.at("Yukawa2").run_as<integral_type>(stg_exponent, aos,
//     aos); auto corr_R = TA::retile(corr.at(property::yukawa), X.trange());
//     REQUIRE(libchemist::ta_helpers::allclose(X, corr_R));
// }
