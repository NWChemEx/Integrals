// #include "integrals/integrals.hpp"
// #include <catch2/catch.hpp>
// #include <mokup/mokup.hpp>
// #include <tensorwrapper/tensor/allclose.hpp>

// using namespace mokup;

// TEST_CASE("Overlap") {
//     using op_type       = simde::type::el_identity;
//     using integral_type = simde::AOTensorRepresentation<2, op_type>;

//     // TODO: Better names which describe what's being tested
//     pluginplay::ModuleManager mm;
//     integrals::load_modules(mm);

//     const auto name = molecule::h2o;
//     const auto bs   = basis_set::sto3g;
//     auto mol        = get_molecule(name);
//     auto aos        = get_bases(name, bs);
//     std::vector bases{bs, bs};
//     auto corr_S = get_ao_data(name, bases, property::overlap);

//     op_type I;
//     auto [S] = mm.at("Overlap").run_as<integral_type>(aos, I, aos);
//     REQUIRE(tensorwrapper::tensor::allclose(S, corr_S));
// }
