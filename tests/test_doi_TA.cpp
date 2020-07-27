#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/doi.hpp>
#include "H2O_STO3G_DOI.hpp"

TEST_CASE("DOI") {
    using integral_type = property_types::DOI<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto [X] = mm.at("DOI").run_as<integral_type>(bs, bs, std::size_t{0});

    REQUIRE(libchemist::allclose(X, TensorType(X.world(), X.trange(), corr)));
}