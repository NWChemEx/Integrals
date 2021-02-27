#include "f12/f12.hpp"
#include "integrals/transformed.hpp"
#include "libint/libint.hpp"

namespace integrals {

// TODO: These keys are going to clobber each other if double and float are
//       loaded.

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

template<typename T>
void load_transformed_libint_integrals(sde::ModuleManager& mm) {
    register_transformed_integral<pt::edipole<T>>(mm, "EDipole");
    register_transformed_integral<pt::equadrupole<T>>(mm, "EQuadrupole");
    register_transformed_integral<pt::eoctopole<T>>(mm, "EOctopole");
    register_transformed_integral<pt::eri2c<T>>(mm, "ERI2");
    register_transformed_integral<pt::eri3c<T>>(mm, "ERI3");
    register_transformed_integral<pt::eri4c<T>>(mm, "ERI4");
    register_transformed_integral<pt::kinetic<T>>(mm, "Kinetic");
    register_transformed_integral<pt::nuclear<T>>(mm, "Nuclear");
    register_transformed_integral<pt::overlap<T>>(mm, "Overlap");
    register_transformed_integral<pt::stg2c<T>>(mm, "STG2");
    register_transformed_integral<pt::stg3c<T>>(mm, "STG3");
    register_transformed_integral<pt::stg4c<T>>(mm, "STG4");
    register_transformed_integral<pt::yukawa2c<T>>(mm, "Yukawa2");
    register_transformed_integral<pt::yukawa3c<T>>(mm, "Yukawa3");
    register_transformed_integral<pt::yukawa4c<T>>(mm, "Yukawa4");
}

template<typename T>
void load_f12_integrals(sde::ModuleManager& mm) {
    mm.add_module<f12::stg_correlation_factor_2c<T>>(
      "STG 2 Center Correlation Factor");
    mm.add_module<f12::stg_correlation_factor_4c<T>>(
      "STG 4 Center Correlation Factor");
    mm.add_module<f12::stg_correlation_factor_squared_4c<T>>(
      "STG 4 Center Correlation Factor Squared");
    mm.add_module<f12::stg_dfdr_squared_4c<T>>("STG 4 Center dfdr Squared");
    mm.add_module<f12::stg_gr2c<T>>("STG 2 Center GR");
    mm.add_module<f12::stg_gr4c<T>>("STG 4 Center GR");
}

template<typename T>
void set_f12_integral_defaults(sde::ModuleManager& mm) {
    mm.change_submod("STG 2 Center Correlation Factor", "STG kernel", "STG2");
    mm.change_submod("STG 4 Center Correlation Factor", "STG kernel", "STG4");
    mm.change_submod("STG 4 Center Correlation Factor Squared", "STG kernel",
                     "STG4");
    mm.change_submod("STG 4 Center dfdr Squared", "STG Kernel", "STG4");
    mm.change_submod("STG 2 Center GR", "Yukawa kernel", "Yukawa2");
    mm.change_submod("STG 4 Center GR", "Yukawa kernel", "Yukawa4");
}

void load_modules(sde::ModuleManager& mm) {
    load_libint_integrals<double>(mm);
    load_transformed_libint_integrals<double>(mm);
    load_f12_integrals<double>(mm);

    // See TODO at top of file before enabling
    // load_libint_integrals<float>(mm);
    // load_transformed_integrals<float>(mm);
    // load_f12_integrals<float>(mm);

    set_f12_integral_defaults<double>(mm);
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
