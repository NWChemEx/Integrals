#include "integrals/integralsmm.hpp"
#include "integrals/libint_integral.hpp"

namespace integrals {
namespace libint {

auto core = integrals::libint::detail_::implementation_type::core;

void load_modules(sde::ModuleManager& mm, std::size_t tile_size) {
    mm.add_module("Overlap", std::make_shared<Overlap>(core));
    mm.at("Overlap").change_input("Tile Size",tile_size);
    mm.add_module("Kinetic", std::make_shared<Kinetic>(core));
    mm.at("Kinetic").change_input("Tile Size",tile_size);
    mm.add_module("Nuclear", std::make_shared<Nuclear>(core));
    mm.at("Nuclear").change_input("Tile Size",tile_size);
    mm.add_module("ERI2", std::make_shared<Metric>(core));
    mm.at("ERI2").change_input("Tile Size",tile_size);
    mm.add_module("ERI3", std::make_shared<DFERI>());
    mm.at("ERI3").change_input("Tile Size",tile_size);
    mm.add_module("ERI4", std::make_shared<ERI>());
    mm.at("ERI4").change_input("Tile Size",tile_size);
    mm.add_module("EDipole", std::make_shared<EDipole>());
    mm.at("EDipole").change_input("Tile Size",tile_size);
    mm.add_module("EQuadrupole", std::make_shared<EQuadrupole>());
    mm.at("EQuadrupole").change_input("Tile Size",tile_size);
    mm.add_module("EOctopole", std::make_shared<EOctopole>());
    mm.at("EOctopole").change_input("Tile Size",tile_size);
}

} // namespace libint
} // namespace integrals
