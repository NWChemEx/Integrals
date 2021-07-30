// #include "cs_screened_integrals.hpp"
// #include "fill_ND_functor.hpp"
// #include "integrals/types.hpp"
// #include "nwx_TA_utils.hpp"
// #include "nwx_libint.hpp"
// #include "traits.hpp"
// #include <libchemist/ta_helpers/ta_hashers.hpp>
// #include <property_types/ao_integrals/type_traits.hpp>
// #include <property_types/cauchy_schwarz_approximation.hpp>
// #include "../unpack_basis_sets.hpp"

// namespace integrals {

// template<typename PropType>
// TEMPLATED_MODULE_CTOR(CauchySchwarzScreened, PropType) {
//     using element_type = double; // TODO: Get from PropType
//     using type::pair_vector;
//     using type::size_vector;
//     using cs_approx_type = property_types::ShellNorms<element_type>;

//     description("Computes an in-core integral with libint");
//     satisfies_property_type<PropType>();

//     add_input<element_type>("Threshold").set_default(1.0E-16);
//     add_input<size_vector>("Tile Size").set_default(size_vector{180});
//     add_input<element_type>("Screening Threshold").set_default(0.0);
//     add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});

//     add_submodule<cs_approx_type>("Shell Norms")
//       .set_description(
//         "Computes the Cauchy-Schwarz Matrix for a pair of basis sets");
// }

// template<typename PropType>
// TEMPLATED_MODULE_RUN(CauchySchwarzScreened, PropType) {
//     using element_type = double; // TODO: Get from PropType
//     using tensor_type  = type::tensor<element_type>;
//     using value_type   = typename tensor_type::value_type;
//     using type::pair_vector;
//     using type::size_vector;
//     using cs_approx_type = property_types::ShellNorms<element_type>;
//     using aospace = property_types::type::ao_space_t<element_type>;

//     auto& world      = TA::get_default_world(); // TODO: Get from runtime
//     auto thresh      = inputs.at("Threshold").value<element_type>();
//     auto tile_size   = inputs.at("Tile Size").value<size_vector>();
//     auto cs_thresh   = inputs.at("Screening
//     Threshold").value<element_type>(); auto atom_ranges = inputs.at("Atom
//     Tile Groups").value<pair_vector>();

//     constexpr auto n_centers =
//       property_types::ao_integrals::n_centers_v<PropType>;
//     constexpr auto op = op_v<PropType>;

//     auto bs     = unpack_basis_sets<PropType>(inputs);
//     auto trange = nwx_TA::select_tiling(bs, tile_size, atom_ranges);
//     auto fill   = nwx_TA::FillNDFunctor<value_type, op, n_centers>();
//     const std::size_t deriv = 0; // TODO: Template on derivative order
//     fill.initialize(nwx_libint::make_basis_sets(bs), deriv, thresh,
//     cs_thresh);

//     // Take care of any special parameters the fill function needs
//     // (should probably we done inside FillNDFunctor to encapsulate the
//     setup) if constexpr(property_types::ao_integrals::is_stg_v<PropType> ||
//                         property_types::ao_integrals::is_yukawa_v<PropType>)
//                         {
//         auto gamma = inputs.at("STG Exponent").value<element_type>();
//         fill.factory.stg_exponent = gamma;
//     }

//     if(cs_thresh > 0.0) {
//         if constexpr(n_centers == 4) {
//             auto bra1 = inputs.at("bra 1").value<aospace>();
//             auto bra2 = inputs.at("bra 1").value<aospace>();
//             auto [cs_mat1] = submods.at("Shell Norms")
//               .run_as<cs_approx_type>(bra1, bra2, deriv);
//             fill.screen.cs_mat1 = cs_mat1;
//         }

//         auto ket1 = inputs.at("ket 1").value<aospace>();
//         auto ket2 = inputs.at("ket 1").value<aospace>();
//         auto [cs_mat2] = submods.at("Shell Norms")
//           .run_as<cs_approx_type>(ket1, ket2, deriv);
//         fill.screen.cs_mat2 = cs_mat2;
//     }

//     auto I  = TiledArray::make_array<tensor_type>(world, trange, fill);
//     auto rv = results();
//     return PropType::wrap_results(rv, I);
// }

// template class CauchySchwarzScreened<pt::eri3c<double>>;
// template class CauchySchwarzScreened<pt::eri4c<double>>;
// template class CauchySchwarzScreened<pt::stg3c<double>>;
// template class CauchySchwarzScreened<pt::stg4c<double>>;
// template class CauchySchwarzScreened<pt::yukawa3c<double>>;
// template class CauchySchwarzScreened<pt::yukawa4c<double>>;

// } // namespace integrals
