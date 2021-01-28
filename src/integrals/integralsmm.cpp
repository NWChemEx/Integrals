#include "libint/libint.hpp"
#include "transformed.hpp"

namespace integrals {

template<typename ElementType>
void load_libint_integrals(sde::ModuleManager& mm) {
    mm.add_module<LibintDOI<ElementType>>("DOI");
    mm.add_module<LibintEDipole<ElementType>>("EDipole");
    mm.add_module<LibintEQuadrupole<ElementType>>("EQuadrupole");
    mm.add_module<LibintEOctopole<ElementType>>("EOctopole");
    mm.add_module<LibintERI2C<ElementType>>("ERI2");
    mm.add_module<LibintERI3C<ElementType>>("ERI3");
    mm.add_module<LibintERI4C<ElementType>>("ERI4");
    mm.add_module<LibintKinetic<ElementType>>("Kinetic");
    mm.add_module<LibintNuclear<ElementType>>("Nuclear");
    mm.add_module<LibintOverlap<ElementType>>("Overlap");
    mm.add_module<LibintSTG2C<ElementType>>("STG2");
    mm.add_module<LibintSTG3C<ElementType>>("STG3");
    mm.add_module<LibintSTG4C<ElementType>>("STG4");
    mm.add_module<LibintYukawa2C<ElementType>>("Yukawa2");
    mm.add_module<LibintYukawa3C<ElementType>>("Yukawa3");
    mm.add_module<LibintYukawa4C<ElementType>>("Yukawa4");
}

void load_modules(sde::ModuleManager& mm) { load_libint_integrals<double>(mm); }

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
