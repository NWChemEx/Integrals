#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <integrals/nwx_direct/stg_direct_type.hpp>
#include "H2O_STO3G_STG[1].hpp"

TEST_CASE("STG4CDirect") {
    using integral_type = property_types::STG4CDirect<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto [X] = mm.at("STG4Direct").run_as<integral_type>(bs, bs, bs, bs, std::size_t{0});

    TensorType real_X(X.world(), X.trange());
    real_X("k, l, m, n") = X("k, l, m, n");

    compare_integrals(real_X, stg1ref, 0.0, 1e-12);
}