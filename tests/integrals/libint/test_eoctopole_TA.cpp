#include "../../test_common_TA.hpp"
#include "integrals/integralsmm.hpp"
#include "integrals/property_types.hpp"
#include "nwx_tesing/H2O_STO3G_OVLP.hpp"
#include "nwx_testing/H2O_STO3G_Multipole.hpp"
#include <libchemist/ta_helpers/ta_helpers.hpp>

using namespace integrals;

TEST_CASE("Octopole") {
    using s_type = pt::overlap<double>;
    using d_type = pt::edipole<double>;
    using q_type = pt::equadrupole<double>;
    using o_type = pt::eoctopole<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto origin         = std::array<double, 3>{0, 0, 0};

    SECTION("Overlap") {
        mm.change_input("EOctopole", "Origin", origin);
        auto [S] = mm.at("EOctopole").run_as<s_type>(bs, bs);
        REQUIRE(libchemist::ta_helpers::allclose(
          S, TensorType(S.world(), S.trange(), corr)));
    }
    SECTION("Dipole") {
        auto [D] = mm.at("EOctopole").run_as<d_type>(origin, bs, bs);
        REQUIRE(libchemist::ta_helpers::allclose(
          D, TensorType(D.world(), D.trange(), corr_dipole)));
    }
    SECTION("Quadrupole") {
        auto [Q] = mm.at("EOctopole").run_as<q_type>(origin, bs, bs);
        REQUIRE(libchemist::ta_helpers::allclose(
          Q, TensorType(Q.world(), Q.trange(), corr_quad)));
    }
    SECTION("Octopole") {
        auto [O] = mm.at("EOctopole").run_as<o_type>(origin, bs, bs);
        REQUIRE(libchemist::ta_helpers::allclose(
          O, TensorType(O.world(), O.trange(), corr_octo)));
    }
}