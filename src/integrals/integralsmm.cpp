#include "f12/f12.hpp"
#include "integrals/integralsmm.hpp"
#include "libint/cs_screened_integrals.hpp"
#include "libint/emultipole.hpp"
#include "libint/libint.hpp"
#include "libint/shellnorms.hpp"
#include "transforms/transforms.hpp"

namespace integrals {

template<std::size_t N, typename OpType>
void register_transformed_integral(pluginplay::ModuleManager& mm,
                                   std::string key) {
    using module_t = StandardTransform<N, OpType>;
    auto new_key   = "Transformed " + key;
    mm.add_module<module_t>(new_key);
    mm.change_submod(new_key, "integral kernel", key);
}

void load_libint_integrals(pluginplay::ModuleManager& mm) {
    mm.add_module<LibintDOI>("DOI");
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
    mm.add_module<CauchySchwarzScreened<simde::ERI3>>("ERI3 CS");
    mm.add_module<CauchySchwarzScreened<simde::ERI4>>("ERI4 CS");
    mm.add_module<CauchySchwarzScreened<simde::STG3>>("STG3 CS");
    mm.add_module<CauchySchwarzScreened<simde::STG4>>("STG4 CS");
    mm.add_module<CauchySchwarzScreened<simde::Yukawa3>>("Yukawa3 CS");
    mm.add_module<CauchySchwarzScreened<simde::Yukawa4>>("Yukawa4 CS");
    mm.add_module<ShellNormCoulomb<double>>("Shell Norms Coulomb");
    mm.add_module<ShellNormSTG<double>>("Shell Norms STG");
    mm.add_module<ShellNormYukawa<double>>("Shell Norms Yukawa");
    mm.change_submod("ERI3 CS", "Shell Norms", "Shell Norms Coulomb");
    mm.change_submod("ERI4 CS", "Shell Norms", "Shell Norms Coulomb");
    mm.change_submod("STG3 CS", "Shell Norms", "Shell Norms STG");
    mm.change_submod("STG4 CS", "Shell Norms", "Shell Norms STG");
    mm.change_submod("Yukawa3 CS", "Shell Norms", "Shell Norms Yukawa");
    mm.change_submod("Yukawa4 CS", "Shell Norms", "Shell Norms Yukawa");
}

void load_transformed_libint_integrals(pluginplay::ModuleManager& mm) {
    // register_transformed_integral<pt::edipole<T>>(mm, "EDipole");
    // register_transformed_integral<pt::equadrupole<T>>(mm, "EQuadrupole");
    // register_transformed_integral<pt::eoctopole<T>>(mm, "EOctopole");
    // register_transformed_integral<pt::eri2c<T>>(mm, "ERI2");
    register_transformed_integral<3, simde::type::el_el_coulomb>(mm, "ERI3");
    register_transformed_integral<4, simde::type::el_el_coulomb>(mm, "ERI4");
    // register_transformed_integral<pt::kinetic<T>>(mm, "Kinetic");
    // register_transformed_integral<pt::nuclear<T>>(mm, "Nuclear");
    // register_transformed_integral<pt::overlap<T>>(mm, "Overlap");
    register_transformed_integral<2, simde::type::el_kinetic>(mm, "Kinetic");
    register_transformed_integral<2, simde::type::el_nuc_coulomb>(mm, "Nuclear");
}

void load_f12_integrals(pluginplay::ModuleManager& mm) {
    mm.add_module<Libint<4, simde::type::el_el_f12_commutator>>(
      "STG 4 Center dfdr Squared");
}

void load_transformed_f12_integrals(pluginplay::ModuleManager& mm) {
    register_transformed_integral<4, simde::type::el_el_f12_commutator>(
      mm, "STG 4 Center dfdr Squared");
    register_transformed_integral<4, simde::type::el_el_stg>(mm, "STG4");
    register_transformed_integral<4, simde::type::el_el_yukawa>(mm, "Yukawa4");
}

void load_misc_transforms(pluginplay::ModuleManager& mm) {
    mm.add_module<StandardTransform<2, simde::type::el_scf_k>>("Transformed K");
    mm.add_module<StandardTransform<2, simde::type::fock>>("Transformed Fock");
}

void load_modules(pluginplay::ModuleManager& mm) {
    load_libint_integrals(mm);
    load_transformed_libint_integrals(mm);
    load_f12_integrals(mm);
    load_transformed_f12_integrals(mm);
    load_misc_transforms(mm);
}

} // namespace integrals
