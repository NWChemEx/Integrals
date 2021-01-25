#include "integrals/libint/libint.hpp"

namespace integrals {

void load_modules(sde::ModuleManager& mm) {
    mm.add_module<LibintDOI<double>>("DOI");
    mm.add_module<LibintEDipole<double>>("EDipole");
    mm.add_module<LibintEQuadrupole<double>>("EQuadrupole");
    mm.add_module<LibintEOctopole<double>>("EOctopole");
    mm.add_module<LibintERI2C<double>>("ERI2");
    mm.add_module<LibintERI3C<double>>("ERI3");
    mm.add_module<LibintERI4C<double>>("ERI4");
    mm.add_module<LibintKinetic<double>>("Kinetic");
    mm.add_module<LibintNuclear<double>>("Nuclear");
    mm.add_module<LibintOverlap<double>>("Overlap");
    mm.add_module<LibintSTG2C<double>>("STG2");
    mm.add_module<LibintSTG3C<double>>("STG3");
    mm.add_module<LibintSTG4C<double>>("STG4");
    mm.add_module<LibintYukawa2C<double>>("Yukawa2");
    mm.add_module<LibintYukawa3C<double>>("Yukawa3");
    mm.add_module<LibintYukawa4C<double>>("Yukawa4");
}

// void load_modules(sde::ModuleManager& mm) {

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
