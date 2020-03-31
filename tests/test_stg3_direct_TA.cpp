#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <integrals/nwx_direct/stg_direct_type.hpp>
#include "H2O_STO3G_STG[1]_3C.hpp"

TEST_CASE("STG3CDirect") {
    using integral_type = property_types::STG3CDirect<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto stg_exponent = 1.0;
    mm.at("STG3Direct").change_input("Tile size", std::vector<std::size_t>{7});
    auto [X] = mm.at("STG3Direct").run_as<integral_type>(bs, bs, bs, std::size_t{0}, stg_exponent);

    TensorType real_X(X.world(), X.trange());
    real_X("k, l, m") = X("k, l, m");

    compare_integrals(real_X, corr);
}