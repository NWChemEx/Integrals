#include "transforms.hpp"

namespace integrals::transforms {

// -----------------------------------------------------------------------------
// -- Define Module Load Functions
// -----------------------------------------------------------------------------

template<std::size_t N, typename OpType>
void register_transformed_integral(pluginplay::ModuleManager& mm,
                                   std::string key) {
    using module_t = StandardTransform<N, OpType>;
    auto new_key   = "Transformed " + key;
    mm.add_module<module_t>(new_key);
    mm.change_submod(new_key, "integral kernel", key);
}

void load_transformed_integrals(pluginplay::ModuleManager& mm) {
    using namespace simde::type;

    // register_transformed_integral<pt::edipole<T>>(mm, "EDipole");
    // register_transformed_integral<pt::equadrupole<T>>(mm, "EQuadrupole");
    // register_transformed_integral<pt::eoctopole<T>>(mm, "EOctopole");
    // register_transformed_integral<pt::eri2c<T>>(mm, "ERI2");
    register_transformed_integral<3, el_el_coulomb>(mm, "ERI3");
    register_transformed_integral<4, el_el_coulomb>(mm, "ERI4");
    // register_transformed_integral<pt::kinetic<T>>(mm, "Kinetic");
    // register_transformed_integral<pt::nuclear<T>>(mm, "Nuclear");
    // register_transformed_integral<pt::overlap<T>>(mm, "Overlap");
    register_transformed_integral<2, el_kinetic>(mm, "Kinetic");
    register_transformed_integral<2, el_nuc_coulomb>(mm, "Nuclear");
    register_transformed_integral<4, el_el_stg>(mm, "STG4");
    register_transformed_integral<4, el_el_yukawa>(mm, "Yukawa4");
    register_transformed_integral<4, el_el_f12_commutator>(
      mm, "STG 4 Center dfdr Squared");

    mm.add_module<StandardTransform<2, el_scf_k>>("Transformed K");
    mm.add_module<StandardTransform<2, fock>>("Transformed Fock");
}

} // namespace integrals::transforms
