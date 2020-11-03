#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/electron_repulsion.hpp>

using matrix_t = TA::detail::matrix_il<double>;

static matrix_t corr{
  {
    1.0464370899978459,
    3.4291996305312606,
    0.0000000000000000,
    0.0000000000000000,
    0.0000000000000000,
    2.6052624057150817,
    2.6052624057150817,
  },
  {
    3.4291996305312606,
    26.435225216427671,
    0.0000000000000000,
    0.0000000000000000,
    0.0000000000000000,
    25.3420821293274088,
    25.3420821293274088,
  },
  {
    0.0000000000000000,
    0.0000000000000000,
    5.7847978365504318,
    0.0000000000000000,
    0.0000000000000000,
    3.2924421173969143,
    3.2924421173969143,
  },
  {
    0.0000000000000000,
    0.0000000000000000,
    0.0000000000000000,
    5.7847978365504300,
    0.0000000000000000,
    0.0000000000000000,
    0.0000000000000000,
  },
  {
    0.0000000000000000,
    0.0000000000000000,
    0.0000000000000000,
    0.0000000000000000,
    5.7847978365504318,
    4.2141100538676941,
    -4.2141100538676941,
  },
  {
    2.6052624057150817,
    25.3420821293274088,
    3.2924421173969143,
    0.0000000000000000,
    4.2141100538676941,
    39.9325707858561643,
    26.6712894368540034,
  },
  {
    2.6052624057150817,
    25.3420821293274088,
    3.2924421173969143,
    0.0000000000000000,
    -4.2141100538676941,
    26.6712894368540034,
    39.9325707858561643,
  },
};

TEST_CASE("ERI2C") {
    using integral_type = property_types::ERI2CIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    mm.at("ERI2").change_input("Tile size", std::vector<std::size_t>{1});
    auto [X] = mm.at("ERI2").run_as<integral_type>(bs, bs, std::size_t{0});

    REQUIRE(libchemist::ta_helpers::allclose(X, TensorType(X.world(), X.trange(), corr)));
}
