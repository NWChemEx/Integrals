#include "../../tensors/H2O_STO3G_Yukawa[1]_3C.hpp"
#include "../../test_common_TA.hpp"
#include "integrals/integrals.hpp"

TEST_CASE("Yukawa3C") {
    using integral_type = integrals::pt::yukawa3c<double>;
    using size_vector   = integrals::type::size_vector;
    const auto key      = "Yukawa3";

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto stg_exponent   = 1.0;
    mm.change_input(key, "Tile size", size_vector{7});
    // mm.change_input(key, "Screening Threshold", 0.000001);
    auto [X] = mm.run_as<integral_type>(key, stg_exponent, bs, bs, bs);
    TensorType corr_R(X.world(), X.trange(), corr);
    REQUIRE(libchemist::ta_helpers::allclose(X, corr_R));
}