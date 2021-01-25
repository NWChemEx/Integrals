#include "fill_ND_functor.hpp"
#include "integrals/libint/libint.hpp"
#include "integrals/types.hpp"
#include "nwx_TA_utils.hpp"
#include "nwx_libint.hpp"
#include "traits.hpp"
#include <property_types/ao_integrals/type_traits.hpp>

namespace integrals {

template<typename PropType>
TEMPLATED_MODULE_CTOR(Libint, PropType) {
    using element_type = double; // TODO: Get from PropType
    using type::pair_vector;
    using type::size_vector;

    description("Computes an in-core integral with libint");
    satisfies_property_type<PropType>();

    add_input<element_type>("Threshold").set_default(1.0E-16);
    add_input<size_vector>("Tile Size").set_default(size_vector{180});
    add_input<element_type>("Screening Threshold").set_default(0.0);
    add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});
}

template<typename PropType>
TEMPLATED_MODULE_RUN(Libint, PropType) {
    using element_type = double; // TODO: Get from PropType
    using tensor_type  = type::tensor<element_type>;
    using value_type   = typename tensor_type::value_type;
    using basis_set    = type::basis_set<element_type>;
    using basis_vector = std::vector<basis_set>;

    using type::pair_vector;
    using type::size_vector;

    auto& world      = TA::get_default_world(); // TODO: Get from runtime
    auto thresh      = inputs.at("Threshold").value<element_type>();
    auto tile_size   = inputs.at("Tile Size").value<size_vector>();
    auto cs_thresh   = inputs.at("Screening Threshold").value<element_type>();
    auto atom_ranges = inputs.at("Atom Tile Groups").value<pair_vector>();

    constexpr auto n_centers =
      property_types::ao_integrals::n_centers_v<PropType>;
    constexpr auto op = op_v<PropType>;

    basis_vector bs;

    // TODO: This logic should be encapsulated in the FillNDFunctor, which
    //       should be specialized on the property type
    if constexpr(is_doi_v<PropType>) {
        // DOI is a bit special in that it's 4 center with 2 basis sets
        auto [bra_space, ket_space] = PropType::unwrap_inputs(inputs);

        auto& bra = bra_space.basis_set();
        auto& ket = ket_space.basis_set();
        bs        = basis_vector{bra, bra, ket, ket};
    } else if constexpr(n_centers == 2) {
        auto [bra_space, ket_space] = PropType::unwrap_inputs(inputs);

        auto& bra = bra_space.basis_set();
        auto& ket = ket_space.basis_set();
        bs        = basis_vector{bra, ket};
    } else if constexpr(n_centers == 3) {
        auto [bra_space, ket1_space, ket2_space] =
          PropType::unwrap_inputs(inputs);

        auto& bra  = bra_space.basis_set();
        auto& ket1 = ket1_space.basis_set();
        auto& ket2 = ket2_space.basis_set();
        bs         = basis_vector{bra, ket1, ket2};
    } else if constexpr(n_centers == 4) {
        auto [bra1_space, bra2_space, ket1_space, ket2_space] =
          PropType::unwrap_inputs(inputs);

        auto& bra1 = bra1_space.basis_set();
        auto& bra2 = bra2_space.basis_set();
        auto& ket1 = ket1_space.basis_set();
        auto& ket2 = ket2_space.basis_set();
        bs         = basis_vector{bra1, bra2, ket1, ket2};
    }

    auto trange = nwx_TA::select_tiling(bs, tile_size, atom_ranges);

    auto fill = nwx_TA::FillNDFunctor<value_type, op, n_centers>();
    fill.initialize(nwx_libint::make_basis_sets(bs), 0, thresh, cs_thresh);

    auto I  = TiledArray::make_array<tensor_type>(world, trange, fill);
    auto rv = results();
    return PropType::wrap_results(rv, I);
}

template class Libint<pt::doi<double>>;
template class Libint<pt::edipole<double>>;
template class Libint<pt::equadrupole<double>>;
template class Libint<pt::eoctopole<double>>;
template class Libint<pt::eri2c<double>>;
template class Libint<pt::eri3c<double>>;
template class Libint<pt::eri4c<double>>;

} // namespace integrals