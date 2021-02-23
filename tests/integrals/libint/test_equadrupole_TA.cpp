#include "../../test_common_TA.hpp"
#include "integrals/integralsmm.hpp"
#include "integrals/property_types.hpp"
#include "nwx_testing/H2O_STO3G_Multipole.hpp"
#include "nwx_testing/H2O_STO3G_OVLP.hpp"
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <property_types/ao_integrals/emultipole.hpp>
#include <property_types/ao_integrals/overlap.hpp>

using namespace integrals;

TEST_CASE("Quadrupole") {
    using s_type = pt::overlap<double>;
    using d_type = pt::edipole<double>;
    using q_type = pt::equadrupole<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto origin         = std::array<double, 3>{0, 0, 0};
    SECTION("Overlap") {
        mm.change_input("EQuadrupole", "Origin", origin);
        auto [S] = mm.at("EQuadrupole").run_as<s_type>(bs, bs);
        REQUIRE(libchemist::ta_helpers::allclose(
          S, TensorType(S.world(), S.trange(), corr)));
    }
    SECTION("Dipole") {
        auto [D] = mm.at("EQuadrupole").run_as<d_type>(origin, bs, bs);
        REQUIRE(libchemist::ta_helpers::allclose(
          D, TensorType(D.world(), D.trange(), corr_dipole)));
    }
    SECTION("Quadrupole") {
        auto [Q] = mm.at("EQuadrupole").run_as<q_type>(origin, bs, bs);
        REQUIRE(libchemist::ta_helpers::allclose(
          Q, TensorType(Q.world(), Q.trange(), corr_quad)));
    }
}
