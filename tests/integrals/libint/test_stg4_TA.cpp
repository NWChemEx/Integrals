// #include "integrals/integrals.hpp"
// #include <catch2/catch.hpp>
// #include <mokup/mokup.hpp>
// #include <tensorwrapper/tensor/allclose.hpp>

// using namespace mokup;

// TEST_CASE("STG4C") {
//     using op_type       = simde::type::el_el_stg;
//     using integral_type = simde::AOTensorRepresentation<4, op_type>;

//     pluginplay::ModuleManager mm;
//     integrals::load_modules(mm);

//     auto name = molecule::h2o;
//     auto bs   = basis_set::sto3g;
//     auto aos  = get_bases(name, bs);
//     std::vector bases{bs, bs, bs, bs};
//     auto corr = get_ao_data(name, bases, property::stg);

//     chemist::Electron e;
//     op_type stg(chemist::operators::STG(1.0, 1.0));
//     auto [X] = mm.at("STG4").run_as<integral_type>(aos, aos, stg, aos, aos);
//     REQUIRE(tensorwrapper::tensor::allclose(X, corr));
// }
