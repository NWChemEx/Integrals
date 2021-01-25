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
    using ao_space_t   = const type::ao_space_t<element_type>&;

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
    // TODO: satisfying a derived property type should automatically satisfy the
    //       the base types. Once that happens use TWoCenter<T>::unwrap etc.
    //       instead of manually grabbing the bras and kets

    // DOI is a bit special in that it's 4 center with 2 basis sets
    if constexpr(is_doi_v<PropType>) {
        auto bra_space = inputs.at("bra").value<ao_space_t>();
        auto ket_space = inputs.at("ket").value<ao_space_t>();

        auto& bra = bra_space.basis_set();
        auto& ket = ket_space.basis_set();
        bs        = basis_vector{bra, bra, ket, ket};
    } else if constexpr(n_centers == 2) {
        auto bra_space = inputs.at("bra").value<ao_space_t>();
        auto ket_space = inputs.at("ket").value<ao_space_t>();

        auto& bra = bra_space.basis_set();
        auto& ket = ket_space.basis_set();
        bs        = basis_vector{bra, ket};
    } else if constexpr(n_centers == 3) {
        auto bra_space  = inputs.at("bra").value<ao_space_t>();
        auto ket1_space = inputs.at("ket 1").value<ao_space_t>();
        auto ket2_space = inputs.at("ket 2").value<ao_space_t>();

        auto& bra  = bra_space.basis_set();
        auto& ket1 = ket1_space.basis_set();
        auto& ket2 = ket2_space.basis_set();
        bs         = basis_vector{bra, ket1, ket2};
    } else if constexpr(n_centers == 4) {
        auto bra1_space = inputs.at("bra 1").value<ao_space_t>();
        auto bra2_space = inputs.at("bra 2").value<ao_space_t>();
        auto ket1_space = inputs.at("ket 1").value<ao_space_t>();
        auto ket2_space = inputs.at("ket 2").value<ao_space_t>();

        auto& bra1 = bra1_space.basis_set();
        auto& bra2 = bra2_space.basis_set();
        auto& ket1 = ket1_space.basis_set();
        auto& ket2 = ket2_space.basis_set();
        bs         = basis_vector{bra1, bra2, ket1, ket2};
    }

    auto trange = nwx_TA::select_tiling(bs, tile_size, atom_ranges);

    auto fill = nwx_TA::FillNDFunctor<value_type, op, n_centers>();
    fill.initialize(nwx_libint::make_basis_sets(bs), 0, thresh, cs_thresh);

    // Take care of any special parameters the fill function needs

    if constexpr(is_nuclear_v<PropType>) {
        using mol_type  = const libchemist::Molecule&;
        const auto& mol = inputs.at("Molecule").value<mol_type>();
        std::vector<std::pair<double, std::array<double, 3>>> qs;
        for(const auto& ai : mol)
            qs.emplace_back(static_cast<const double&>(ai.Z()), ai.coords());
        fill.factory.qs = qs;
    } else if constexpr(is_stg_v<PropType> || is_yukawa_v<PropType>) {
        auto gamma = inputs.at("STG Exponent").value<element_type>();
        fill.factory.stg_exponent = gamma;
    }

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
template class Libint<pt::kinetic<double>>;
template class Libint<pt::nuclear<double>>;
template class Libint<pt::overlap<double>>;
template class Libint<pt::stg2c<double>>;
template class Libint<pt::stg3c<double>>;
template class Libint<pt::stg4c<double>>;
template class Libint<pt::yukawa2c<double>>;
template class Libint<pt::yukawa3c<double>>;
template class Libint<pt::yukawa4c<double>>;

} // namespace integrals