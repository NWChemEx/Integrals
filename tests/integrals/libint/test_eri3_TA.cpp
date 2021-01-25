#include "../../tensors/H2O_STO3G_DF.hpp"
#include "../../test_common_TA.hpp"
#include "integrals/integrals.hpp"

TEST_CASE("ERI3C") {
    using integral_type = integrals::pt::eri3c<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    //     mm.at("ERI3").change_input("Screening Threshold", 0.000001);
    //     auto [X] = mm.at("ERI3").run_as<integral_type>(bs, bs, bs);

    //     REQUIRE(libchemist::ta_helpers::allclose(
    //       X, TensorType(X.world(), X.trange(), corr)));
}
