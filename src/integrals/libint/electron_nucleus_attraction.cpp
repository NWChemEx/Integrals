
// #include "detail_/fill_ND_functor.hpp"
// #include "detail_/nwx_TA_utils.hpp"
// #include "detail_/nwx_libint.hpp"
// #include "detail_/type_traits.hpp"

// using ao_space_t  = simde::type::ao_space;
// using ao_ref_wrap = std::reference_wrapper<const ao_space_t>;
// using size_type   = std::size_t;
// using size_vector = std::vector<size_vector>;
// using pair_vector = std::vector<std::pair<size_type, size_type>>;

// namespace integrals {

// template<>
// TEMPLATED_MODULE_CTOR(Libint, 2, simde::type::el_nuc_coulomb) {
//     using my_pt = simde::AOTensorRepresentation<2,
//     simde::type::el_nuc_coulomb>;

//     add_input<double>("Threshold").set_default(1.0E-16);
//     add_input<size_vector>("Tile Size").set_default(size_vector{180});
//     add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});
// }

// template<>
// TEMPLATED_MODULE_RUN(Libint, 2, simde::type::el_nuc_coulomb) {
//     using my_pt = simde::AOTensorRepresentation<2,
//     simde::type::el_nuc_coulomb>;

//     auto& world      = TA::get_default_world(); // TODO: Get from runtime
//     auto thresh      = inputs.at("Threshold").value<element_type>();
//     auto tile_size   = inputs.at("Tile Size").value<size_vector>();
//     auto atom_ranges = inputs.at("Atom Tile Groups").value<pair_vector>();

//     auto [bra, op, ket]      = my_pt::unpack_inputs(inputs);
//     constexpr auto libint_op = detail_::op_v<simde::type::el_nuc_coulomb>;
//     auto trange = nwx_TA::select_tiling(bs, tile_size, atom_ranges);
//     auto fill   = nwx_TA::FillNDFunctor<value_type, op, n_centers>();
//     const double cs_cthresh = 0.0; // Just to satisfy initialize
//     fill.initialize(nwx_libint::make_basis_sets(bs), deriv, thresh,
//     cs_thresh);

//     const auto& pot = op.

//                       if constexpr(is_el_nuc) {
//         using mol_type  = const libchemist::Molecule&;
//         const auto& mol = inputs.at("Molecule").value<mol_type>();
//         std::vector<std::pair<double, std::array<double, 3>>> qs;
//         for(const auto& ai : mol)
//             qs.emplace_back(static_cast<const double&>(ai.Z()), ai.coords());
//         fill.factory.qs = qs;
//     }
//     // else if constexpr(property_types::ao_integrals::is_stg_v<PropType> ||
//     //                   property_types::ao_integrals::is_yukawa_v<PropType>)
//     {
//     //     auto gamma = inputs.at("STG Exponent").value<element_type>();
//     //     fill.factory.stg_exponent = gamma;
//     // }

//     using tensor_type = TA::SparseArrayD;
//     auto I  = TiledArray::make_array<tensor_type>(world, trange, fill);
//     auto rv = results();
//     return PropType::wrap_results(rv, simde::tensor(I));
// }
// template class Libint<2, simde::type::el_nuc_coulomb>;

// } // namespace integrals
