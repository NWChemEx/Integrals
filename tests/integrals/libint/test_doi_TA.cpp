// #include "integrals/integrals.hpp"
// #include <catch2/catch.hpp>
// #include <mokup/mokup.hpp>
// #include <tensorwrapper/tensor/allclose.hpp>

// using namespace mokup;

// TEST_CASE("DOI") {
//     using op_type       = simde::type::el_el_delta;
//     using integral_type = simde::EDOI;

//     pluginplay::ModuleManager mm;
//     integrals::load_modules(mm);

//     const auto name = molecule::h2o;
//     const auto bs   = basis_set::sto3g;
//     auto mol        = get_molecule(name);
//     auto aos        = get_bases(name, bs);
//     std::vector bases{bs, bs};
//     auto corr = get_ao_data(name, bases, property::dois);
//     op_type d;
//     auto [X] = mm.at("DOI").run_as<integral_type>(aos, d, aos);
//     REQUIRE(tensorwrapper::tensor::allclose(X, corr));
// }
