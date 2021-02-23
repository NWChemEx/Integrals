#include "../../test_common_TA.hpp"
#include "integrals/integralsmm.hpp"
#include "integrals/property_types.hpp"
#include "nwx_tesing/H2O_STO3G_OVLP.hpp"
#include "nwx_testing/H2O_STO3G_Multipole.hpp"
#include <libchemist/ta_helpers/ta_helpers.hpp>

using namespace integrals;

TEST_CASE("Dipole") {
    using s_type = pt::overlap<double>;
    using d_type = pt::edipole<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    std::array<double, 3> origin{0, 0, 0};

    SECTION("overlap matrix") {
        mm.change_input("EDipole", "Origin", origin);
        auto [S] = mm.at("EDipole").run_as<s_type>(bs, bs);
        TensorType corr_s(S.world(), S.trange(), corr);

        REQUIRE(libchemist::ta_helpers::allclose(
          S, TensorType(S.world(), S.trange(), corr)));
    }

    SECTION("dipole matrix") {
        auto [D] = mm.at("EDipole").run_as<d_type>(origin, bs, bs);
        REQUIRE(libchemist::ta_helpers::allclose(
          D, TensorType(D.world(), D.trange(), corr_dipole)));
    }
}
