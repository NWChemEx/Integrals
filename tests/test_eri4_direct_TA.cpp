#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <integrals/nwx_direct/eri_direct_type.hpp>
#include "H2O_STO3G_ERI.hpp"

TEST_CASE("ERI4CDirect") {
    using integral_type = property_types::ERI4CDirect<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto [I] = mm.at("ERI4Direct").run_as<integral_type>(bs, bs, bs, bs, std::size_t{0});

    TensorType real_I(I.world(), I.trange());
    real_I("k, l, m, n") = I("k, l, m, n");

    compare_integrals(real_I, corr);
}
