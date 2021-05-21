#include "../../test_common_TA.hpp"
#include "integrals/integralsmm.hpp"
#include "integrals/property_types.hpp"
#include "testing/testing.hpp"
#include <libchemist/libchemist.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>

using namespace testing;

TEST_CASE("STG 4 Center GR") {
    using integral_type = integrals::pt::gr4c<double>;
    const auto key      = "STG 4 Center GR";

    auto& world = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2;
    const auto bs   = basis_set::sto3g;

    for(const auto& bs : {basis_set::sto3g, basis_set::ccpvdz}) {
        std::vector<basis_set> bs_key(4, bs);
        SECTION(as_string(name) + "/" + as_string(bs)) {
            auto aos     = get_bases().at(name).at(bs);
            auto tensors = get_ao_data(world).at(name).at(bs_key);
            auto X_corr  = tensors.at(property::stg_gr);
            auto [X]     = mm.at(key).run_as<integral_type>(aos, aos, aos, aos);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_corr));
        }
    }
}