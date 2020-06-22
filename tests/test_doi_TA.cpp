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

    compare_integrals(X, corr);

    mm.copy_module("DOI", "DOI_1");
    std::vector<std::vector<std::size_t>> atom_groups{{0, 1, 2}};
    mm.at("DOI_1").change_input("Atom Tile Groups", atom_groups);
    auto [X_1] = mm.at("DOI_1").run_as<integral_type>(bs, bs, std::size_t{0});

    compare_integrals(X_1, corr);
}