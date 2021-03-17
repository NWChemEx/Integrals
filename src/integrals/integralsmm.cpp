#include "integrals/transformed.hpp"
#include "libint/cs_screened_integrals.hpp"
#include "libint/libint.hpp"
#include "libint/shellnorms.hpp"

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

    mm.add_module<ScreenedERI3C<ElementType>>("ERI3 CS");
    mm.add_module<ScreenedERI4C<ElementType>>("ERI4 CS");
    mm.add_module<ScreenedSTG3C<ElementType>>("STG3 CS");
    mm.add_module<ScreenedSTG4C<ElementType>>("STG4 CS");
    mm.add_module<ScreenedYukawa3C<ElementType>>("Yukawa3 CS");
    mm.add_module<ScreenedYukawa4C<ElementType>>("Yukawa4 CS");
    mm.add_module<ShellNormCoulomb<ElementType>>("Shell Norms Coulomb");
    mm.add_module<ShellNormSTG<ElementType>>("Shell Norms STG");
    mm.add_module<ShellNormYukawa<ElementType>>("Shell Norms Yukawa");
    mm.change_submod("ERI3 CS", "Shell Norms", "Shell Norms Coulomb");
    mm.change_submod("ERI4 CS", "Shell Norms", "Shell Norms Coulomb");
    mm.change_submod("STG3 CS", "Shell Norms", "Shell Norms STG");
    mm.change_submod("STG4 CS", "Shell Norms", "Shell Norms STG");
    mm.change_submod("Yukawa3 CS", "Shell Norms", "Shell Norms Yukawa");
    mm.change_submod("Yukawa4 CS", "Shell Norms", "Shell Norms Yukawa");
}

template<typename T>
void load_transformed_integrals(sde::ModuleManager& mm) {
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

#undef REGISTER_TRANSFORMED_INT

void load_modules(sde::ModuleManager& mm) {
    load_libint_integrals<double>(mm);
    load_transformed_integrals<double>(mm);
}

} // namespace integrals
