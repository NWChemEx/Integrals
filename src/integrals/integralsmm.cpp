#include "f12/f12.hpp"
#include "integrals/integralsmm.hpp"
#include "libint/cs_screened_integrals.hpp"
#include "libint/emultipole.hpp"
#include "libint/libint.hpp"
#include "libint/shellnorms.hpp"

namespace integrals {

void load_libint_integrals(pluginplay::ModuleManager& mm) {
    // mm.add_module<LibintDOI<ElementType>>("DOI");
    mm.add_module<LibintDipole>("EDipole");
    mm.add_module<LibintQuadrupole>("EQuadrupole");
    mm.add_module<LibintOctupole>("EOctupole");
    mm.add_module<Libint<2, simde::type::el_el_coulomb>>("ERI2");
    mm.add_module<Libint<3, simde::type::el_el_coulomb>>("ERI3");
    mm.add_module<Libint<4, simde::type::el_el_coulomb>>("ERI4");
    mm.add_module<Libint<2, simde::type::el_kinetic>>("Kinetic");
    mm.add_module<Libint<2, simde::type::el_nuc_coulomb>>("Nuclear");
    mm.add_module<Libint<2, simde::type::el_identity>>("Overlap");
    mm.add_module<Libint<2, simde::type::el_el_stg>>("STG2");
    mm.add_module<Libint<3, simde::type::el_el_stg>>("STG3");
    mm.add_module<Libint<4, simde::type::el_el_stg>>("STG4");
    mm.add_module<Libint<2, simde::type::el_el_yukawa>>("Yukawa2");
    mm.add_module<Libint<3, simde::type::el_el_yukawa>>("Yukawa3");
    mm.add_module<Libint<4, simde::type::el_el_yukawa>>("Yukawa4");

    // mm.add_module<ScreenedERI3C<ElementType>>("ERI3 CS");
    // mm.add_module<ScreenedERI4C<ElementType>>("ERI4 CS");
    // mm.add_module<ScreenedSTG3C<ElementType>>("STG3 CS");
    // mm.add_module<ScreenedSTG4C<ElementType>>("STG4 CS");
    // mm.add_module<ScreenedYukawa3C<ElementType>>("Yukawa3 CS");
    // mm.add_module<ScreenedYukawa4C<ElementType>>("Yukawa4 CS");
    // mm.add_module<ShellNormCoulomb<ElementType>>("Shell Norms Coulomb");
    // mm.add_module<ShellNormSTG<ElementType>>("Shell Norms STG");
    // mm.add_module<ShellNormYukawa<ElementType>>("Shell Norms Yukawa");
    // mm.change_submod("ERI3 CS", "Shell Norms", "Shell Norms Coulomb");
    // mm.change_submod("ERI4 CS", "Shell Norms", "Shell Norms Coulomb");
    // mm.change_submod("STG3 CS", "Shell Norms", "Shell Norms STG");
    // mm.change_submod("STG4 CS", "Shell Norms", "Shell Norms STG");
    // mm.change_submod("Yukawa3 CS", "Shell Norms", "Shell Norms Yukawa");
    // mm.change_submod("Yukawa4 CS", "Shell Norms", "Shell Norms Yukawa");
}

void load_transformed_libint_integrals(pluginplay::ModuleManager& mm) {
    // register_transformed_integral<pt::edipole<T>>(mm, "EDipole");
    // register_transformed_integral<pt::equadrupole<T>>(mm, "EQuadrupole");
    // register_transformed_integral<pt::eoctopole<T>>(mm, "EOctopole");
    // register_transformed_integral<pt::eri2c<T>>(mm, "ERI2");
    // register_transformed_integral<pt::eri3c<T>>(mm, "ERI3");
    // register_transformed_integral<pt::eri4c<T>>(mm, "ERI4");
    // register_transformed_integral<pt::kinetic<T>>(mm, "Kinetic");
    // register_transformed_integral<pt::nuclear<T>>(mm, "Nuclear");
    // register_transformed_integral<pt::overlap<T>>(mm, "Overlap");
    // register_transformed_integral<pt::stg2c<T>>(mm, "STG2");
    // register_transformed_integral<pt::stg3c<T>>(mm, "STG3");
    // register_transformed_integral<pt::stg4c<T>>(mm, "STG4");
    // register_transformed_integral<pt::yukawa2c<T>>(mm, "Yukawa2");
    // register_transformed_integral<pt::yukawa3c<T>>(mm, "Yukawa3");
    // register_transformed_integral<pt::yukawa4c<T>>(mm, "Yukawa4");
}

void load_f12_integrals(pluginplay::ModuleManager& mm) {
    // mm.add_module<f12::stg_dfdr_squared_4c<T>>("STG 4 Center dfdr Squared");
}

void load_transformed_f12_integrals(pluginplay::ModuleManager& mm) {
    // register_transformed_integral<pt::dfdr_squared_4c<T>>(
    //   mm, "STG 4 Center dfdr Squared");
}

void set_f12_integral_defaults(pluginplay::ModuleManager& mm) {
    // mm.change_submod("STG 4 Center dfdr Squared", "STG Kernel", "STG4");
}

void load_modules(pluginplay::ModuleManager& mm) {
    load_libint_integrals(mm);
    load_transformed_libint_integrals(mm);
    load_f12_integrals(mm);
    load_transformed_f12_integrals(mm);
    set_f12_integral_defaults(mm);
}

} // namespace integrals
