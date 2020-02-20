#include "integrals/integralsmm.hpp"
#include "integrals/overlap_integral.hpp"
#include "integrals/kinetic_integral.hpp"
#include "integrals/nuclear_integral.hpp"
#include "integrals/eri_integrals.hpp"
#include "integrals/stg_integrals.hpp"
#include "integrals/yukawa_integrals.hpp"
#include "integrals/emultipole_integrals.hpp"
#include "integrals/doi.hpp"

namespace integrals {

void load_modules(sde::ModuleManager& mm) {
    mm.add_module("Overlap", std::make_shared<Overlap>());

    mm.add_module("Kinetic", std::make_shared<Kinetic>());

    mm.add_module("Nuclear", std::make_shared<Nuclear>());

    mm.add_module("ERI2", std::make_shared<ERI2>());

    mm.add_module("ERI3", std::make_shared<ERI3>());

    mm.add_module("ERI4", std::make_shared<ERI4>());

    mm.add_module("STG2", std::make_shared<STG2>());

    mm.add_module("STG3", std::make_shared<STG3>());

    mm.add_module("STG4", std::make_shared<STG4>());

    mm.add_module("Yukawa2", std::make_shared<Yukawa2>());

    mm.add_module("Yukawa3", std::make_shared<Yukawa3>());

    mm.add_module("Yukawa4", std::make_shared<Yukawa4>());

    mm.add_module("EDipole", std::make_shared<EDipole>());

    mm.add_module("EQuadrupole", std::make_shared<EQuadrupole>());

    mm.add_module("EOctopole", std::make_shared<EOctopole>());

    mm.add_module("DOI", std::make_shared<DOI>());
}

} // namespace integrals
