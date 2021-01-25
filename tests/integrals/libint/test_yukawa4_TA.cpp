#include "../../tensors/H2O_STO3G_Yukawa[1].hpp"
#include "../../test_common_TA.hpp"
#include "integrals/integrals.hpp"

TEST_CASE("Yukawa4C") {
    using integral_type = integrals::pt::yukawa4c<double>;
    using size_vector   = integrals::type::size_vector;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto stg_exponent   = 1.0;
    mm.at("Yukawa4").change_input("Screening Threshold", 0.000001);
    auto [X] =
      mm.at("Yukawa4").run_as<integral_type>(stg_exponent, bs, bs, bs, bs);
    TensorType corr_R(X.world(), X.trange(), corr);
    REQUIRE(libchemist::ta_helpers::allclose(X, corr_R));
}
