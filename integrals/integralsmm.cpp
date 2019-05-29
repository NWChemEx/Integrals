#include "integrals/integralsmm.hpp"
#include "integrals/libint_integral.hpp"

namespace integrals {
namespace libint {


void load_modules(sde::ModuleManager& mm,
                  const std::size_t tile_size,
                  const detail_::implementation_type& impl,
                  const double screen) {
    mm.add_module("Overlap", std::make_shared<Overlap>(impl));
    mm.at("Overlap").change_input("Tile Size",tile_size);

    mm.add_module("Kinetic", std::make_shared<Kinetic>(impl));
    mm.at("Kinetic").change_input("Tile Size",tile_size);

    mm.add_module("Nuclear", std::make_shared<Nuclear>(impl));
    mm.at("Nuclear").change_input("Tile Size",tile_size);

    mm.add_module("ERI2", std::make_shared<Metric>(impl));
    mm.at("ERI2").change_input("Tile Size",tile_size);

    mm.add_module("ERI3", std::make_shared<DFERI>(impl));
    mm.at("ERI3").change_input("Tile Size",tile_size);
    mm.at("ERI3").change_input("Screening Threshold",screen);

    mm.add_module("ERI4", std::make_shared<ERI>(impl));
    mm.at("ERI4").change_input("Tile Size",tile_size);
    mm.at("ERI4").change_input("Screening Threshold",screen);

    mm.add_module("EDipole", std::make_shared<EDipole>(impl));
    mm.at("EDipole").change_input("Tile Size",tile_size);

    mm.add_module("EQuadrupole", std::make_shared<EQuadrupole>(impl));
    mm.at("EQuadrupole").change_input("Tile Size",tile_size);

    mm.add_module("EOctopole", std::make_shared<EOctopole>(impl));
    mm.at("EOctopole").change_input("Tile Size",tile_size);
}

} // namespace libint
} // namespace integrals
