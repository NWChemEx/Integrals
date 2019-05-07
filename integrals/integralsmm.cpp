#include "integrals/integralsmm.hpp"
#include "integrals/libint_integral.hpp"

namespace integrals {
namespace libint {

void load_modules(sde::ModuleManager& mm) {
    mm.add_module("Overlap", std::make_shared<Overlap>());
    mm.add_module("Kinetic", std::make_shared<Kinetic>());
    mm.add_module("Nuclear", std::make_shared<Nuclear>());
    mm.add_module("ERI2", std::make_shared<Metric>());
    mm.add_module("ERI3", std::make_shared<DFERI>());
    mm.add_module("ERI4", std::make_shared<ERI>());
    mm.add_module("EDipole", std::make_shared<EDipole>());
    mm.add_module("EQuadrupole", std::make_shared<EQuadrupole>());
    mm.add_module("EOctopole", std::make_shared<EOctopole>());
}

} // namespace libint
} // namespace integrals
