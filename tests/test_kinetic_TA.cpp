#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/kinetic.hpp>

static BlockTensor corr{
        {{0, 0,}, {29.00319994553958, -0.1680109393164922, 0.0000000000000000, 0.0000000000000000,  0.0000000000000000,
                   -0.1680109393164923, 0.8081279549303477, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
                   0.0000000000000000, 0.0000000000000000, 2.528731198194765, 0.0000000000000000, 0.0000000000000000,
                   0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 2.528731198194765, 0.0000000000000000,
                   0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 2.528731198194765, } },
        {{0, 1,}, {-0.008416385185447427, -0.008416385185447427, 0.07051733851899882, 0.07051733851899882, 0.1149203802569082,
                   0.1149203802569082, 0.0000000000000000, 0.0000000000000000, 0.1470905524127557, -0.1470905524127557, } },
        {{1, 0,}, {-0.008416385185447427, 0.07051733851899882, 0.1149203802569082,  0.0000000000000000, 0.1470905524127557,
                   -0.008416385185447427, 0.07051733851899882, 0.1149203802569082, 0.0000000000000000,  -0.1470905524127557, } },
        {{1, 1,}, {0.760031883566609, -0.003979736727037247, -0.003979736727037247, 0.760031883566609, } },
};

TEST_CASE("Kinetic") {
    using integral_type = property_types::KineticIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    mm.at("Kinetic").change_input("Tile size", std::vector<std::size_t>{1, 2});
    auto [T] = mm.at("Kinetic").run_as<integral_type>(bs, bs, std::size_t{0});

    compare_integrals(T, corr);

    mm.copy_module("Kinetic", "Kinetic_1");
    std::vector<std::vector<std::size_t>> atom_groups{{0}, {1, 2}};
    mm.at("Kinetic_1").change_input("Atom Tile Groups", atom_groups);
    auto [T_1] = mm.at("Kinetic_1").run_as<integral_type>(bs, bs, std::size_t{0});

    compare_integrals(T_1, corr);
}
