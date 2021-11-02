#include "detail_/run_fundamental.hpp"
#include "standard.hpp"
#include <simde/simde.hpp>
#include <utilities/strings/string_tools.hpp>

namespace integrals {

template<std::size_t N, typename OpType>
TEMPLATED_MODULE_CTOR(StandardTransform, N, OpType) {
    using pt     = simde::TransformedTensorRepresentation<N, OpType>;
    using sub_pt = simde::AOTensorRepresentation<N, OpType>;

    satisfies_property_type<pt>();

    add_submodule<sub_pt>("integral kernel");
}

template<std::size_t N, typename OpType>
TEMPLATED_MODULE_RUN(StandardTransform, N, OpType) {
    using pt     = simde::TransformedTensorRepresentation<N, OpType>;
    using sub_pt = simde::AOTensorRepresentation<N, OpType>;

    const auto& [ao_spaces, derived_spaces, op] = pt::unwrap_inputs(inputs);

    if(ao_spaces.size() + derived_spaces.size() != N)
        throw std::runtime_error("Did not get enough (or got too many) bases.");

    // Grab the AO basis sets for each mode
    std::decay_t<decltype(ao_spaces)> all_aos;
    for(const auto& [mode, aos] : ao_spaces) all_aos.emplace(mode, aos);
    for(const auto& [mode, mos] : derived_spaces)
        all_aos.emplace(mode, std::cref(mos.get().from_space()));

    // Get the integral in the AO basis
    auto& submod = submods.at("integral kernel");
    auto t       = detail_::run_fundamental<N>(submod.value(), all_aos, op);

    // Sort transforms
    std::multimap<std::size_t, std::size_t> size2mode;
    for(const auto& [mode, mos] : derived_spaces) {
        size2mode.emplace(mos.get().size(), mode);
    }

    // Do the transforms
    simde::type::tensor temp;
    const auto out_idx = t.make_annotation();
    for(const auto& [size, mode] : size2mode) {
        const auto& C          = derived_spaces.at(mode).get().C();
        std::string mode2trans = "i" + std::to_string(mode);
        auto in_idx   = utilities::strings::replace(mode2trans, "mu", out_idx);
        temp(out_idx) = C("mu," + mode2trans) * t(in_idx);
        t             = temp;
    }

    auto rv = results();
    return pt::wrap_results(rv, t);
}

template class StandardTransform<2, simde::type::el_scf_k>;
template class StandardTransform<2, simde::type::fock>;
template class StandardTransform<3, simde::type::el_el_coulomb>;
template class StandardTransform<4, simde::type::el_el_coulomb>;
template class StandardTransform<4, simde::type::el_el_f12_commutator>;
template class StandardTransform<4, simde::type::el_el_stg>;
template class StandardTransform<4, simde::type::el_el_yukawa>;

} // namespace integrals
