#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/ao_integrals/nuclear.hpp>

static BlockTensor corr{
        {{0, 0,}, {-61.5805952694322, -7.410821856331163, -0.0144738837457361, 0.0000000000000000, 0.0000000000000000, -1.231685572142488,
                   -7.410821856331163, -10.00907114207003, -0.1768908347336431, 0.0000000000000000, 0.0000000000000000, -2.977226853578134,
                   -0.01447388374573611, -0.1768908347336431, -9.944043341698766, 0.0000000000000000, 0.0000000000000000, -1.471793338712961,
                   0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -9.875875995090944, 0.0000000000000000, 0.0000000000000000,
                   0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -9.987549935088563, -1.822236913476131,
                   -1.231685572142488, -2.977226853578134, -1.471793338712961, 0.0000000000000000, -1.822236913476131, -5.300203252295022, } },
        {{0, 1,}, {-1.231685572142488, -2.977226853578134, -1.471793338712961, 0.0000000000000000, 1.822236913476131, -1.067171080472437, } },
        {{1, 0,}, {-1.231685572142488, -2.977226853578134, -1.471793338712961, 0.0000000000000000,  1.82223691347613, -1.067171080472437, } },
        {{1, 1,}, {-5.30020325229502, } },
};

TEST_CASE("Nuclear") {
    using integral_type = property_types::NuclearIntegral<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    mm.at("Nuclear").change_input("Tile size", std::vector<std::size_t>{6, 1});
    auto [V] = mm.at("Nuclear").run_as<integral_type>(bs, bs, molecule, std::size_t{0});

    compare_integrals(V, corr);

    mm.copy_module("Nuclear", "Nuclear_1");
    std::vector<std::vector<std::size_t>> atom_groups{{0, 1}, {2}};
    mm.at("Nuclear_1").change_input("Atom Tile Groups", atom_groups);
    auto [V_1] = mm.at("Nuclear_1").run_as<integral_type>(bs, bs, molecule, std::size_t{0});

    compare_integrals(V_1, corr);
}
