#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/overlap.hpp>

static BlockTensor corr{
  {{0, 0,},
   {
     1.0000000000000004, 0.2367039365108476,  0.0000000000000000,
     0.0000000000000000, 0.0000000000000000,  0.2367039365108476,
     1.0000000000000002, -0.0000000000000000, 0.0000000000000000,
     0.0000000000000000, 0.0000000000000000,  -0.0000000000000000,
     1.0000000000000007, 0.0000000000000000,  0.0000000000000000,
     0.0000000000000000, 0.0000000000000000,  0.0000000000000000,
     1.0000000000000002, 0.0000000000000000,  0.0000000000000000,
     0.0000000000000000, 0.0000000000000000,  0.0000000000000000,
     1.0000000000000007,
   }},
  {{0, 1,},
   {
     0.0384055905135491,
     0.3861387813310929,
     0.2097276494226498,
     0.0000000000000000,
     0.2684376412681763,
   }},
  {{0, 2,},
   {
     0.0384055905135491,
     0.3861387813310929,
     0.2097276494226498,
     0.0000000000000000,
     -0.2684376412681763,
   }},
  {{1, 0,},
   {
     0.0384055905135491,
     0.3861387813310928,
     0.2097276494226498,
     0.0000000000000000,
     0.2684376412681763,
   }},
  {{1, 1,},
   {
     1.0000000000000002,
   }},
  {{1, 2,},
   {
     0.1817608668218930,
   }},
  {{2, 0,},
   {
     0.0384055905135491,
     0.3861387813310928,
     0.2097276494226498,
     0.0000000000000000,
     -0.2684376412681763,
   }},
  {{2, 1,},
   {
     0.1817608668218930,
   }},
  {{2, 2,},
   {
     1.0000000000000002,
   }}
};

TEST_CASE("Overlap") {
    using integral_type = property_types::OverlapIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    mm.at("Overlap").change_input("Tile size", std::vector<std::size_t>{3, 1, 1});
    auto [molecule, bs] = make_molecule();
    auto [S] = mm.at("Overlap").run_as<integral_type>(bs, bs, std::size_t{0});

    compare_integrals(S, corr);

    mm.copy_module("Overlap", "Overlap_1");
    std::vector<std::vector<std::size_t>> atom_groups{{0}, {1}, {2}};
    mm.at("Overlap_1").change_input("Atom Tile Groups", atom_groups);
    auto [S_1] = mm.at("Overlap_1").run_as<integral_type>(bs, bs, std::size_t{0});

    compare_integrals(S_1, corr);
}
