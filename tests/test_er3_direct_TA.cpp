#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <integrals/nwx_direct/eri_direct_type.hpp>
#include "H2O_STO3G_DF.hpp"

TEST_CASE("ERI3CDirect") {
    using integral_type = property_types::ERI3CDirect<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto [X] = mm.at("ERI3Direct").run_as<integral_type>(bs, bs, bs, std::size_t{0});

    TensorType real_X(X.world(), X.trange());
    real_X("k, l, m") = X("k, l, m");

    compare_integrals(real_X, corr);
}
