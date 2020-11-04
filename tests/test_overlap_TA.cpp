#include "test_common_TA.hpp"
#include "H2O_STO3G_OVLP.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/overlap.hpp>

TEST_CASE("Overlap") {
    using integral_type = property_types::OverlapIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.at("Overlap").change_input("Tile size",
                                  std::vector<std::size_t>{3, 1, 1});
    auto [molecule, bs] = make_molecule();
    auto [S] = mm.at("Overlap").run_as<integral_type>(bs, bs, std::size_t{0});

    REQUIRE(libchemist::ta_helpers::allclose(S, TensorType(S.world(), S.trange(), corr)));

    mm.copy_module("Overlap", "Overlap_1");
    std::vector<std::pair<std::size_t, std::size_t>> atom_groups{
      {0, 1}, {1, 2}, {2, 3}};
    mm.at("Overlap_1").change_input("Tile size", std::vector<std::size_t>{100});
    mm.at("Overlap_1").change_input("Atom Tile Groups", atom_groups);
    auto [S_1] =
      mm.at("Overlap_1").run_as<integral_type>(bs, bs, std::size_t{0});

    REQUIRE(
      libchemist::ta_helpers::allclose(S_1, TensorType(S_1.world(), S_1.trange(), corr)));
}
