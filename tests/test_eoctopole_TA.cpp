#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/overlap.hpp>
#include <property_types/ao_integrals/emultipole.hpp>
#include "H2O_STO3G_Multipole.hpp"
#include "H2O_STO3G_OVLP.hpp"

TEST_CASE("Octopole") {
    using s_type = property_types::OverlapIntegral<double>;
    using d_type = property_types::EDipoleIntegral<double>;
    using q_type = property_types::EQuadrupoleIntegral<double>;
    using o_type = property_types::EOctopoleIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto origin = std::array<double, 3>{0,0,0};
    auto [S] = mm.at("EDipole").run_as<s_type>(bs, bs, std::size_t{0});
    auto [D] = mm.at("EDipole").run_as<d_type>(bs, bs, std::size_t{0}, origin);
    auto [Q] = mm.at("EQuadrupole").run_as<q_type>(bs, bs, std::size_t{0}, origin);
    auto [O] = mm.at("EOctopole").run_as<o_type>(bs, bs, std::size_t{0}, origin);

    REQUIRE(libchemist::allclose(S, TensorType(S.world(), S.trange(), corr)));
    REQUIRE(libchemist::allclose(D, TensorType(D.world(), D.trange(), corr_dipole)));
    REQUIRE(libchemist::allclose(Q, TensorType(Q.world(), Q.trange(), corr_quad)));
    REQUIRE(libchemist::allclose(O, TensorType(O.world(), O.trange(), corr_octo)));
}