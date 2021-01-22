#include "integrals/libint/libint.hpp"

namespace integrals {

void load_modules(sde::ModuleManager& mm) {
    mm.add_module<LibintDOI<double>>("DOI");
    mm.add_module<LibintEDipole<double>>("EDipole");
    mm.add_module<LibintEQuadrupole<double>>("EQuadrupole");
    mm.add_module<LibintEOctopole<double>>("EOctopole");
}

// void load_modules(sde::ModuleManager& mm) {
//     mm.add_module("Overlap", std::make_shared<Overlap>());

//     mm.add_module("Kinetic", std::make_shared<Kinetic>());

//     mm.add_module("Nuclear", std::make_shared<Nuclear>());

//     mm.add_module("ERI2", std::make_shared<ERI2>());

//     mm.add_module("ERI3", std::make_shared<ERI3>());

//     mm.add_module("ERI4", std::make_shared<ERI4>());

//     mm.add_module("STG2", std::make_shared<STG2>());

//     mm.add_module("STG3", std::make_shared<STG3>());

//     mm.add_module("STG4", std::make_shared<STG4>());

//     mm.add_module("Yukawa2", std::make_shared<Yukawa2>());

//     mm.add_module("Yukawa3", std::make_shared<Yukawa3>());

//     mm.add_module("Yukawa4", std::make_shared<Yukawa4>());

//     mm.add_module("ERI3Direct", std::make_shared<ERI3Direct>());

//     mm.add_module("ERI4Direct", std::make_shared<ERI4Direct>());

//     mm.add_module("STG3Direct", std::make_shared<STG3Direct>());

//     mm.add_module("STG4Direct", std::make_shared<STG4Direct>());

//     mm.add_module("Yukawa3Direct", std::make_shared<Yukawa3Direct>());

//     mm.add_module("Yukawa4Direct", std::make_shared<Yukawa4Direct>());

//     mm.add_module("CauchySchwarzERI", std::make_shared<CS_ERI>());

//     mm.add_module("CauchySchwarzSTG", std::make_shared<CS_STG>());

//     mm.add_module("CauchySchwarzYukawa", std::make_shared<CS_Yukawa>());

//     mm.at("ERI3").change_submod("Cauchy-Schwarz",
//     std::make_shared<sde::Module>(
//                                                     mm.at("CauchySchwarzERI")));
//     mm.at("ERI4").change_submod("Cauchy-Schwarz",
//     std::make_shared<sde::Module>(
//                                                     mm.at("CauchySchwarzERI")));
//     mm.at("ERI3Direct")
//       .change_submod("Cauchy-Schwarz",
//                      std::make_shared<sde::Module>(mm.at("CauchySchwarzERI")));
//     mm.at("ERI4Direct")
//       .change_submod("Cauchy-Schwarz",
//                      std::make_shared<sde::Module>(mm.at("CauchySchwarzERI")));

//     mm.at("STG3").change_submod("Cauchy-Schwarz",
//     std::make_shared<sde::Module>(
//                                                     mm.at("CauchySchwarzSTG")));
//     mm.at("STG4").change_submod("Cauchy-Schwarz",
//     std::make_shared<sde::Module>(
//                                                     mm.at("CauchySchwarzSTG")));
//     mm.at("STG3Direct")
//       .change_submod("Cauchy-Schwarz",
//                      std::make_shared<sde::Module>(mm.at("CauchySchwarzSTG")));
//     mm.at("STG4Direct")
//       .change_submod("Cauchy-Schwarz",
//                      std::make_shared<sde::Module>(mm.at("CauchySchwarzSTG")));

//     mm.at("Yukawa3").change_submod(
//       "Cauchy-Schwarz",
//       std::make_shared<sde::Module>(mm.at("CauchySchwarzYukawa")));
//     mm.at("Yukawa4").change_submod(
//       "Cauchy-Schwarz",
//       std::make_shared<sde::Module>(mm.at("CauchySchwarzYukawa")));
//     mm.at("Yukawa3Direct")
//       .change_submod("Cauchy-Schwarz", std::make_shared<sde::Module>(
//                                          mm.at("CauchySchwarzYukawa")));
//     mm.at("Yukawa4Direct")
//       .change_submod("Cauchy-Schwarz", std::make_shared<sde::Module>(
//                                          mm.at("CauchySchwarzYukawa")));
// }

} // namespace integrals
