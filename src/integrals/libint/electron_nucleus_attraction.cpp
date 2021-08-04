#include "detail_/fill_ND_functor.hpp"
#include "detail_/nwx_TA_utils.hpp"
#include "detail_/nwx_libint.hpp"
#include "detail_/type_traits.hpp"
#include "libint.hpp"
#include <simde/tensor_representation/tensor_representation.hpp>

using ao_space_t  = simde::type::ao_space;
using ao_ref_wrap = std::reference_wrapper<const ao_space_t>;
using size_type   = std::size_t;
using size_vector = std::vector<size_type>;
using pair_vector = std::vector<std::pair<size_type, size_type>>;
using tensor_type = TA::TSpArrayD;
using value_type  = typename tensor_type::value_type;

namespace integrals {

MODULE_CTOR(LibintDOI) {
    using my_pt = simde::AOTensorRepresentation<2, simde::type::el_el_delta>;

    satisfies_property_type<my_pt>();

    add_input<double>("Threshold").set_default(1.0E-16);
    add_input<size_vector>("Tile Size").set_default(size_vector{180});
    add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});
}

MODULE_RUN(LibintDOI) {
    using my_pt = simde::AOTensorRepresentation<2, simde::type::el_el_delta>;

    auto& world      = TA::get_default_world(); // TODO: Get from runtime
    auto thresh      = inputs.at("Threshold").value<double>();
    auto tile_size   = inputs.at("Tile Size").value<size_vector>();
    auto atom_ranges = inputs.at("Atom Tile Groups").value<pair_vector>();

    auto [bra, op, ket]      = my_pt::unwrap_inputs(inputs);
    constexpr auto libint_op = libint2::Operator::delta;
    const auto& bra_bs       = bra.basis_set();
    const auto& ket_bs       = ket.basis_set();
    std::vector<simde::type::ao_basis_set> bs{bra_bs, bra_bs, ket_bs, ket_bs};
    auto trange            = nwx_TA::select_tiling(bs, tile_size, atom_ranges);
    auto fill              = nwx_TA::FillNDFunctor<value_type, libint_op, 4>();
    const double cs_thresh = 0.0; // Just to satisfy initialize
    fill.initialize(nwx_libint::make_basis_sets(bs), 0, thresh, cs_thresh);

    auto I  = TiledArray::make_array<tensor_type>(world, trange, fill);
    auto rv = results();
    return my_pt::wrap_results(rv, simde::type::tensor(I));
}

} // namespace integrals
