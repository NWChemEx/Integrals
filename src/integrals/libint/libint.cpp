#include "detail_/fill_ND_functor.hpp"
#include "detail_/nwx_TA_utils.hpp"
#include "detail_/nwx_libint.hpp"
#include "detail_/special_setup.hpp"
#include "detail_/type_traits.hpp"
#include "libint.hpp"
#include <simde/tensor_representation/ao_tensor_representation.hpp>

using ao_space_t  = simde::type::ao_space;
using ao_ref_wrap = std::reference_wrapper<const ao_space_t>;
using size_type   = std::size_t;
using size_vector = std::vector<size_type>;
using pair_vector = std::vector<std::pair<size_type, size_type>>;

namespace integrals {
namespace {

template<std::size_t N, typename ModuleInputs>
auto unpack_bases(const ModuleInputs& inputs) {
    using ao_space_ref = const ao_space_t&;
    std::vector<ao_space_t> aos(N);
    if constexpr(N == 2) {
        aos[0] = inputs.at("bra").template value<ao_space_t>();
        aos[1] = inputs.at("ket").template value<ao_space_t>();
    } else if constexpr(N == 3) {
        aos[0] = inputs.at("bra").template value<ao_space_t>();
        aos[1] = inputs.at("ket 1").template value<ao_space_t>();
        aos[2] = inputs.at("ket 2").template value<ao_space_t>();
    } else if constexpr(N == 4) {
        aos[0] = inputs.at("bra 1").template value<ao_space_t>();
        aos[1] = inputs.at("bra 2").template value<ao_space_t>();
        aos[2] = inputs.at("ket 1").template value<ao_space_t>();
        aos[3] = inputs.at("ket 2").template value<ao_space_t>();
    }
    std::vector<libchemist::AOBasisSet<double>> rv;
    for(auto i = 0u; i < N; ++i) rv.emplace_back(aos[i].basis_set());
    return rv;
}

} // namespace

template<std::size_t N, typename OperatorType>
TEMPLATED_MODULE_CTOR(Libint, N, OperatorType) {
    using my_pt = simde::AOTensorRepresentation<N, OperatorType>;

    satisfies_property_type<my_pt>();

    add_input<double>("Threshold").set_default(1.0E-16);
    add_input<size_vector>("Tile Size").set_default(size_vector{180});
    add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});
}

template<std::size_t N, typename OperatorType>
TEMPLATED_MODULE_RUN(Libint, N, OperatorType) {
    using my_pt = simde::AOTensorRepresentation<N, OperatorType>;

    auto& world      = TA::get_default_world(); // TODO: Get from runtime
    auto thresh      = inputs.at("Threshold").value<double>();
    auto tile_size   = inputs.at("Tile Size").value<size_vector>();
    auto atom_ranges = inputs.at("Atom Tile Groups").value<pair_vector>();

    auto aos                 = unpack_bases<N>(inputs);
    auto op                  = inputs.at("op").value<const OperatorType&>();
    constexpr auto libint_op = integrals::op_v<OperatorType>;
    auto trange            = nwx_TA::select_tiling(aos, tile_size, atom_ranges);
    using tensor_type      = TA::TSpArrayD;
    using value_type       = typename tensor_type::value_type;
    auto fill              = nwx_TA::FillNDFunctor<value_type, libint_op, N>();
    const double cs_thresh = 0.0; // Just to satisfy initialize
    fill.initialize(nwx_libint::make_basis_sets(aos), 0, thresh, cs_thresh);

    SpecialSetup<OperatorType>::setup(fill, op);

    // else if constexpr(property_types::ao_integrals::is_stg_v<PropType> ||
    //                   property_types::ao_integrals::is_yukawa_v<PropType>) {
    //     auto gamma = inputs.at("STG Exponent").value<element_type>();
    //     fill.factory.stg_exponent = gamma;
    // }

    auto I  = TiledArray::make_array<tensor_type>(world, trange, fill);
    auto rv = results();
    return my_pt::wrap_results(rv, simde::type::tensor(I));
}

template class Libint<2, simde::type::el_el_coulomb>;
template class Libint<3, simde::type::el_el_coulomb>;
template class Libint<4, simde::type::el_el_coulomb>;
template class Libint<2, simde::type::el_kinetic>;
template class Libint<2, simde::type::el_nuc_coulomb>;
} // namespace integrals
