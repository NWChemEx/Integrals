#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <integrals/nwx_direct/eri_direct_type.hpp>
#include "H2O_STO3G_ERI.hpp"

TEST_CASE("ERI4CDirect") {
    using integral_type = property_types::ERI4CDirect<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto [X] = mm.at("ERI4Direct").run_as<integral_type>(bs, bs, bs, bs, std::size_t{0});

    TensorType real_X(X.world(), X.trange());
    real_X("k, l, m, n") = X("k, l, m, n");

    REQUIRE(libchemist::ta_helpers::allclose(real_X, TensorType(real_X.world(), real_X.trange(), corr)));
}
