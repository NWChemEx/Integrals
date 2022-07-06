// #include "detail_/aos2shells.hpp"
// #include "detail_/bases_helper.hpp"
// #include "detail_/make_engine.hpp"
// #include "detail_/make_shape.hpp"
// #include "detail_/shells2ord.hpp"
// #include "libint.hpp"
// #include <simde/tensor_representation/ao_tensor_representation.hpp>

// /// TODO: Unify implementations. Maybe with recursion?

// namespace integrals {

// using identity_op   = simde::type::el_identity;
// using dipole_op     = simde::type::el_dipole;
// using quadrupole_op = simde::type::el_quadrupole;
// using octupole_op   = simde::type::el_octupole;

// using overlap_pt    = simde::AOTensorRepresentation<2, identity_op>;
// using dipole_pt     = simde::AOTensorRepresentation<2, dipole_op>;
// using quadrupole_pt = simde::AOTensorRepresentation<2, quadrupole_op>;
// using octupole_pt   = simde::AOTensorRepresentation<2, octupole_op>;

// /// Grab the various detail_ functions
// using namespace detail_;

// MODULE_CTOR(LibintDipole) {
//     description("Computes in-core dipole integrals with libint");

//     satisfies_property_type<overlap_pt>();
//     identity_op I;
//     change_input(I.as_string()).change(std::move(I));

//     satisfies_property_type<dipole_pt>();
//     dipole_op r;
//     change_input(r.as_string()).change(std::move(r));

//     add_input<double>("Threshold")
//       .set_default(1.0E-16)
//       .set_description(
//         "The target precision with which the integrals will be computed");
// }

// MODULE_RUN(LibintDipole) {
//     /// Typedefs
//     using my_pt         = dipole_pt;
//     using size_vector_t = std::vector<std::size_t>;
//     using tensor_t      = simde::type::tensor;
//     using field_t       = typename tensor_t::field_type;

//     /// Grab input information
//     auto bases  = unpack_bases<2>(inputs);
//     auto op_str = dipole_op().as_string();
//     auto op     = inputs.at(op_str).template value<const dipole_op&>();
//     auto thresh = inputs.at("Threshold").value<double>();

//     /// Lambda to calculate values
//     auto l = [&](const auto& lo, const auto& up, auto* data) {
//         /// Convert index values from AOs to shells
//         size_vector_t lo_shells, up_shells;
//         for(auto i = 0; i < N; ++i) {
//             auto shells_in_tile = aos2shells(bases[i], lo[i], up[i]);
//             lo_shells.push_back(shells_in_tile.front());
//             up_shells.push_back(shells_in_tile.back());
//         }

//         /// Make the libint engine to calculate integrals
//         auto engine     = make_engine(bases, op, thresh);
//         const auto& buf = engine.results();

//         /// Loop through shell combinations
//         size_vector_t curr_shells = lo_shells;
//         while(curr_shells[0] <= up_shells[0]) {
//             /// Determine which values will be computed this time
//             auto ord_pos = shells2ord(bases, curr_shells);

//             /// Compute values
//             run_engine_(engine, bases, curr_shells,
//                         std::make_index_sequence<N>());
//             auto vals = buf[0];

//             /// Copy libint values into tile data;
//             for(auto i = 0; i < ord_pos.size(); ++i) {
//                 data[ord_pos[i]] = vals[i];
//             }

//             /// Increment curr_shells
//             curr_shells[N - 1] += 1;
//             for(auto i = 1; i < N; ++i) {
//                 if(curr_shells[N - i] > up_shells[N - i]) {
//                     /// Reset this dimension and increment the next one
//                     /// curr_shells[0] accumulates until we reach the end
//                     curr_shells[N - i] = lo_shells[N - i];
//                     curr_shells[N - i - 1] += 1;
//                 }
//             }
//         }
//     };
//     tensor_t I(l, make_shape(bases),
//                tensorwrapper::tensor::default_allocator<field_t>());


//     auto rv = results();
//     // rv      = overlap_pt::wrap_results(rv, simde::type::tensor(S));
//     rv = dipole_pt::wrap_results(rv, simde::type::tensor(D));
//     return rv;
// }

// MODULE_CTOR(LibintQuadrupole) {
//     description("Computes an in-core integral with libint");
    
//     satisfies_property_type<overlap_pt>();
//     identity_op I;
//     change_input(I.as_string()).change(std::move(I));

//     satisfies_property_type<dipole_pt>();
//     dipole_op r;
//     change_input(r.as_string()).change(std::move(r));

//     satisfies_property_type<quadrupole_pt>();
//     quadrupole_op r2;
//     change_input(r2.as_string()).change(std::move(r2));

//     add_input<double>("Threshold")
//       .set_default(1.0E-16)
//       .set_description(
//         "The target precision with which the integrals will be computed");
// }

// MODULE_RUN(LibintQuadrupole) {
//     auto [bra_space, r2, ket_space] = quadrupole_pt::unwrap_inputs(inputs);
//     auto thresh                     = inputs.at("Threshold").value<double>();
//     auto tile_size   = inputs.at("Tile Size").value<size_vector>();
//     auto cs_thresh   = inputs.at("Screening Threshold").value<double>();
//     auto atom_ranges = inputs.at("Atom Tile Groups").value<pair_vector>();

//     auto& bra   = bra_space.basis_set();
//     auto& ket   = ket_space.basis_set();
//     auto& world = TA::get_default_world();

//     constexpr auto libint_op = libint2::Operator::emultipole2;
//     auto fill = nwx_TA::FillNDFunctor<value_type, libint_op, 2>();
//     auto bs   = nwx_libint::make_basis_sets({bra, ket});
//     fill.initialize(bs, 0, thresh, cs_thresh);

//     for(size_type i = 0; i < 3; ++i)
//         fill.factory.origin[i] = r2.gauge_origin().coord(i);

//     auto nopers          = libint2::operator_traits<libint_op>::nopers;
//     auto component_range = nwx_TA::make_tiled_range(nopers, nopers);
//     auto trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
//                                         {component_range});

//     auto X = TiledArray::make_array<tensor_type>(world, trange, fill);
//     auto separate_comps = nwx_TA::make_tiled_range(nopers, {1, 3, 6});
//     trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
//                                    {separate_comps});
//     X      = TA::retile(X, trange);

//     // Separate out components
//     tensor_type D, Q;
//     auto upper      = trange.tiles_range().upbound();
//     using size_type = long;
//     using il_type   = std::initializer_list<size_type>;
//     il_type lo_D{1, 0, 0}, hi_D{2, upper[1], upper[2]};
//     il_type lo_Q{2, 0, 0}, hi_Q{3, upper[1], upper[2]};

//     D("i,j,k") = X("i,j,k").block(lo_D, hi_D);
//     Q("i,j,k") = X("i,j,k").block(lo_Q, hi_Q);

//     // Make overlap 2D
//     // tensor_type S;
//     // il_type lo_S{0, 0, 0}, hi_S{1, upper[1], upper[2]};
//     // S("i,j,k") = X("i,j,k").block(lo_S, hi_S);
//     // auto I = TA::diagonal_array<tensor_type, element_type>(
//     //   world, TA::TiledRange{S.trange().dim(0)});
//     // S = chemist::ta_helpers::einsum::einsum("j,k", "i,j,k", "i", S, I);

//     auto rv = results();
//     // rv      = overlap_pt::wrap_results(rv, S);
//     rv = dipole_pt::wrap_results(rv, simde::type::tensor(D));
//     rv = quadrupole_pt::wrap_results(rv, simde::type::tensor(Q));
//     return rv;
// }

// MODULE_CTOR(LibintOctupole) {
//     description("Computes an in-core integral with libint");

//     satisfies_property_type<overlap_pt>();
//     identity_op I;
//     change_input(I.as_string()).change(std::move(I));

//     satisfies_property_type<dipole_pt>();
//     dipole_op r;
//     change_input(r.as_string()).change(std::move(r));

//     satisfies_property_type<quadrupole_pt>();
//     quadrupole_op r2;
//     change_input(r2.as_string()).change(std::move(r2));

//     satisfies_property_type<octupole_pt>();
//     octupole_op r3;
//     change_input(r3.as_string()).change(std::move(r3));

//     add_input<double>("Threshold")
//       .set_default(1.0E-16)
//       .set_description(
//         "The target precision with which the integrals will be computed");
// }

// MODULE_RUN(LibintOctupole) {
//     auto [bra_space, r3, ket_space] = octupole_pt::unwrap_inputs(inputs);
//     auto thresh                     = inputs.at("Threshold").value<double>();
//     auto tile_size   = inputs.at("Tile Size").value<size_vector>();
//     auto cs_thresh   = inputs.at("Screening Threshold").value<double>();
//     auto atom_ranges = inputs.at("Atom Tile Groups").value<pair_vector>();

//     auto& bra   = bra_space.basis_set();
//     auto& ket   = ket_space.basis_set();
//     auto& world = TA::get_default_world();

//     constexpr auto libint_op = libint2::Operator::emultipole3;
//     auto fill = nwx_TA::FillNDFunctor<value_type, libint_op, 2>();
//     auto bs   = nwx_libint::make_basis_sets({bra, ket});
//     fill.initialize(bs, 0, thresh, cs_thresh);

//     for(size_type i = 0; i < 3; ++i)
//         fill.factory.origin[i] = r3.gauge_origin().coord(i);

//     auto nopers          = libint2::operator_traits<libint_op>::nopers;
//     auto component_range = nwx_TA::make_tiled_range(nopers, nopers);
//     auto trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
//                                         {component_range});

//     auto X = TiledArray::make_array<tensor_type>(world, trange, fill);
//     auto separate_comps = nwx_TA::make_tiled_range(nopers, {1, 3, 6, 10});
//     trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
//                                    {separate_comps});
//     X      = TA::retile(X, trange);

//     // Separate out components
//     tensor_type D, Q, O;
//     auto upper      = trange.tiles_range().upbound();
//     using size_type = long;
//     using il_type   = std::initializer_list<size_type>;

//     il_type lo_D{1, 0, 0}, hi_D{2, upper[1], upper[2]};
//     il_type lo_Q{2, 0, 0}, hi_Q{3, upper[1], upper[2]};
//     il_type lo_O{3, 0, 0}, hi_O{4, upper[1], upper[2]};

//     D("i,j,k") = X("i,j,k").block(lo_D, hi_D);
//     Q("i,j,k") = X("i,j,k").block(lo_Q, hi_Q);
//     O("i,j,k") = X("i,j,k").block(lo_O, hi_O);

//     // Make overlap 2D
//     // tensor_type S;
//     // il_type lo_S{0, 0, 0}, hi_S{1, upper[1], upper[2]};
//     // S("i,j,k") = X("i,j,k").block(lo_S, hi_S);
//     // auto I = TA::diagonal_array<tensor_type, element_type>(
//     //   world, TA::TiledRange{S.trange().dim(0)});
//     // S = chemist::ta_helpers::einsum::einsum("j,k", "i,j,k", "i", S, I);

//     auto rv = results();
//     // rv      = overlap_pt::wrap_results(rv, S);
//     rv = dipole_pt::wrap_results(rv, simde::type::tensor(D));
//     rv = quadrupole_pt::wrap_results(rv, simde::type::tensor(Q));
//     rv = octupole_pt::wrap_results(rv, simde::type::tensor(O));
//     return rv;
// }

// } // namespace integrals
