
#include "detail_/fill_ND_functor.hpp"
#include "detail_/nwx_TA_utils.hpp"
#include "detail_/nwx_libint.hpp"
#include "detail_/special_setup.hpp"
#include "detail_/type_traits.hpp"
#include <simde/tensor_representations/ao_tensor_representation.hpp>

using ao_space_t  = simde::type::ao_space;
using ao_ref_wrap = std::reference_wrapper<const ao_space_t>;
using size_type   = std::size_t;
using size_vector = std::vector<size_vector>;
using pair_vector = std::vector<std::pair<size_type, size_type>>;

namespace integrals {
namespace {

template<std::size_t N>
auto unpack_bases(const ModuleInputs& inputs) {
    using ao_space_ref = const ao_space_t&;
    std::array<ao_ref_wrap, N> rv;
    if constexpr(N == 2) {
        rv[0] = std::cref(inputs.at("bra").value<ao_space_ref>());
        rv[1] = std::cref(inputs.at("ket").value<ao_space_ref>());
    } else if constexpr(N == 3) {
        rv[0] = std::cref(inputs.at("bra").value<ao_space_ref>());
        rv[1] = std::cref(inputs.at("ket 1").value<ao_space_ref>());
        rv[2] = std::cref(inputs.at("ket 2").value<ao_space_ref>());
    } else if constexpr(N == 4) {
        rv[0] = std::cref(inputs.at("bra 1").value<ao_space_ref>());
        rv[1] = std::cref(inptus.at("bra 2").value<ao_space_ref>());
        rv[2] = std::cref(inputs.at("ket 1").value<ao_space_ref>());
        rv[3] = std::cref(inputs.at("ket 2").value<ao_space_ref>());
    }
    return rv;
}

} // namespace

template<std::size_t N, typename OperatorType>
TEMPLATED_MODULE_CTOR(Libint, N, OperatorType) {
    using my_pt = simde::AOTensorRepresentation<N, OperatorType>;

    add_input<double>("Threshold").set_default(1.0E-16);
    add_input<size_vector>("Tile Size").set_default(size_vector{180});
    add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});
}

template<std::size_t N, typename OperatorType>
TEMPLATED_MODULE_RUN(Libint, N, OperatorType) {
    using my_pt = simde::AOTensorRepresentation<N, OperatorType>;

    auto& world      = TA::get_default_world(); // TODO: Get from runtime
    auto thresh      = inputs.at("Threshold").value<element_type>();
    auto tile_size   = inputs.at("Tile Size").value<size_vector>();
    auto atom_ranges = inputs.at("Atom Tile Groups").value<pair_vector>();

    auto aos                 = unpack_bases<T>(inputs);
    auto op                  = inputs.at("op").value<const OperatorType&>();
    constexpr auto libint_op = detail_::op_v<OperatorType>;
    auto trange = nwx_TA::select_tiling(bs, tile_size, atom_ranges);
    auto fill   = nwx_TA::FillNDFunctor<value_type, op, n_centers>();
    const double cs_cthresh = 0.0; // Just to satisfy initialize
    fill.initialize(nwx_libint::make_basis_sets(bs), deriv, thresh, cs_thresh);

    SpecialSetup<OperatorType>::setup(fill, op);

    // else if constexpr(property_types::ao_integrals::is_stg_v<PropType> ||
    //                   property_types::ao_integrals::is_yukawa_v<PropType>) {
    //     auto gamma = inputs.at("STG Exponent").value<element_type>();
    //     fill.factory.stg_exponent = gamma;
    // }

    using tensor_type = TA::SparseArrayD;
    auto I  = TiledArray::make_array<tensor_type>(world, trange, fill);
    auto rv = results();
    return PropType::wrap_results(rv, simde::tensor(I));
}

template class Libint<2, simde::type::el_el_coulomb>;
template class Libint<3, simde::type::el_el_coulomb>;
template class Libint<4, simde::type::el_el_coulomb>;
template class Libint<2, simde::type::el_kinetic>;

} // namespace integrals
