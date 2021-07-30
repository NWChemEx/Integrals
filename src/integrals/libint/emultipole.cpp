// #include "fill_ND_functor.hpp"
// #include "integrals/types.hpp"
// #include "libint.hpp"
// #include "nwx_TA_utils.hpp"
// #include "nwx_libint.hpp"
// #include "traits.hpp"
// #include <libchemist/ta_helpers/einsum/einsum.hpp>

// // TODO: Unify implementations. Maybe with recursion?

// namespace integrals {

// template<typename T>
// TEMPLATED_MODULE_CTOR(Libint, pt::edipole<T>) {
//     using element_type = T;
//     using type::pair_vector;
//     using type::size_vector;

//     description("Computes an in-core integral with libint");
//     satisfies_property_type<pt::overlap<T>>();
//     satisfies_property_type<pt::edipole<T>>();

//     add_input<element_type>("Threshold").set_default(1.0E-16);
//     add_input<size_vector>("Tile Size").set_default(size_vector{180});
//     add_input<element_type>("Screening Threshold").set_default(0.0);
//     add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});
// }

// template<typename T>
// TEMPLATED_MODULE_RUN(Libint, pt::edipole<T>) {
//     using element_type = T;
//     using edipole_pt   = pt::edipole<T>;
//     using overlap_pt   = pt::overlap<T>;
//     using tensor_type  = type::tensor<T>;
//     using value_type   = typename tensor_type::value_type;
//     using type::pair_vector;
//     using type::size_vector;

//     auto [origin, bra_space, ket_space] = edipole_pt::unwrap_inputs(inputs);

//     auto& bra         = bra_space.basis_set();
//     auto& ket         = ket_space.basis_set();
//     std::size_t deriv = 0;

//     auto thresh      = inputs.at("Threshold").value<element_type>();
//     auto tile_size   = inputs.at("Tile Size").value<size_vector>();
//     auto cs_thresh   = inputs.at("Screening
//     Threshold").value<element_type>(); auto atom_ranges = inputs.at("Atom
//     Tile Groups").value<pair_vector>(); auto& world      =
//     TA::get_default_world();

//     auto fill =
//       nwx_TA::FillNDFunctor<value_type, libint2::Operator::emultipole1, 2>();
//     fill.initialize(nwx_libint::make_basis_sets({bra, ket}), deriv, thresh,
//                     cs_thresh);
//     fill.factory.origin = origin;

//     auto nopers =
//       libint2::operator_traits<libint2::Operator::emultipole1>::nopers;
//     auto component_range = nwx_TA::make_tiled_range(nopers, nopers);
//     auto trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
//                                         {component_range});

//     auto X = TiledArray::make_array<tensor_type>(world, trange, fill);
//     auto separate_comps = nwx_TA::make_tiled_range(nopers, {1, 3});
//     trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
//                                    {separate_comps});
//     X      = TA::retile(X, trange);

//     // Separate out components
//     tensor_type S, D;
//     auto upper      = trange.tiles_range().upbound();
//     using size_type = long;
//     using il_type   = std::initializer_list<size_type>;
//     il_type lo_S{0, 0, 0}, hi_S{1, upper[1], upper[2]};
//     il_type lo_D{1, 0, 0}, hi_D{2, upper[1], upper[2]};

//     S("i,j,k") = X("i,j,k").block(lo_S, hi_S);
//     D("i,j,k") = X("i,j,k").block(lo_D, hi_D);

//     // Make overlap 2D
//     auto I = TA::diagonal_array<tensor_type, element_type>(
//       world, TA::TiledRange{S.trange().dim(0)});
//     S = libchemist::ta_helpers::einsum::einsum("j,k", "i,j,k", "i", S, I);

//     auto rv = results();
//     rv      = overlap_pt::wrap_results(rv, S);
//     rv      = edipole_pt::wrap_results(rv, D);
//     return rv;
// }

// template<typename T>
// TEMPLATED_MODULE_CTOR(Libint, pt::equadrupole<T>) {
//     using element_type = T;
//     using type::pair_vector;
//     using type::size_vector;

//     description("Computes an in-core integral with libint");
//     satisfies_property_type<pt::overlap<T>>();
//     satisfies_property_type<pt::edipole<T>>();
//     satisfies_property_type<pt::equadrupole<T>>();

//     add_input<element_type>("Threshold").set_default(1.0E-16);
//     add_input<size_vector>("Tile Size").set_default(size_vector{180});
//     add_input<element_type>("Screening Threshold").set_default(0.0);
//     add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});
// }

// template<typename T>
// TEMPLATED_MODULE_RUN(Libint, pt::equadrupole<T>) {
//     using element_type   = T;
//     using overlap_pt     = pt::overlap<T>;
//     using edipole_pt     = pt::edipole<T>;
//     using equadrupole_pt = pt::equadrupole<T>;
//     using eoctopole_pt   = pt::eoctopole<T>;
//     using tensor_type    = type::tensor<T>;
//     using value_type     = typename tensor_type::value_type;
//     using type::pair_vector;
//     using type::size_vector;

//     auto [origin, bra_space, ket_space] =
//     equadrupole_pt::unwrap_inputs(inputs);

//     auto& bra         = bra_space.basis_set();
//     auto& ket         = ket_space.basis_set();
//     std::size_t deriv = 0;

//     auto thresh      = inputs.at("Threshold").value<element_type>();
//     auto tile_size   = inputs.at("Tile Size").value<size_vector>();
//     auto cs_thresh   = inputs.at("Screening
//     Threshold").value<element_type>(); auto atom_ranges = inputs.at("Atom
//     Tile Groups").value<pair_vector>(); auto& world      =
//     TA::get_default_world();

//     auto fill =
//       nwx_TA::FillNDFunctor<value_type, libint2::Operator::emultipole2, 2>();
//     fill.initialize(nwx_libint::make_basis_sets({bra, ket}), deriv, thresh,
//                     cs_thresh);
//     fill.factory.origin = origin;

//     auto nopers =
//       libint2::operator_traits<libint2::Operator::emultipole2>::nopers;
//     auto component_range = nwx_TA::make_tiled_range(nopers, nopers);
//     auto trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
//                                         {component_range});

//     auto X = TiledArray::make_array<tensor_type>(world, trange, fill);
//     auto separate_comps = nwx_TA::make_tiled_range(nopers, {1, 3, 6});
//     trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
//                                    {separate_comps});
//     X      = TA::retile(X, trange);

//     // Separate out components
//     tensor_type S, D, Q;
//     auto upper      = trange.tiles_range().upbound();
//     using size_type = long;
//     using il_type   = std::initializer_list<size_type>;
//     il_type lo_S{0, 0, 0}, hi_S{1, upper[1], upper[2]};
//     il_type lo_D{1, 0, 0}, hi_D{2, upper[1], upper[2]};
//     il_type lo_Q{2, 0, 0}, hi_Q{3, upper[1], upper[2]};

//     S("i,j,k") = X("i,j,k").block(lo_S, hi_S);
//     D("i,j,k") = X("i,j,k").block(lo_D, hi_D);
//     Q("i,j,k") = X("i,j,k").block(lo_Q, hi_Q);

//     // Make overlap 2D
//     auto I = TA::diagonal_array<tensor_type, element_type>(
//       world, TA::TiledRange{S.trange().dim(0)});
//     S = libchemist::ta_helpers::einsum::einsum("j,k", "i,j,k", "i", S, I);

//     auto rv = results();
//     rv      = overlap_pt::wrap_results(rv, S);
//     rv      = edipole_pt::wrap_results(rv, D);
//     rv      = equadrupole_pt::wrap_results(rv, Q);
//     return rv;
// }

// template<typename T>
// TEMPLATED_MODULE_CTOR(Libint, pt::eoctopole<T>) {
//     using element_type = T;
//     using type::pair_vector;
//     using type::size_vector;

//     description("Computes an in-core integral with libint");
//     satisfies_property_type<pt::overlap<T>>();
//     satisfies_property_type<pt::edipole<T>>();
//     satisfies_property_type<pt::equadrupole<T>>();
//     satisfies_property_type<pt::eoctopole<T>>();

//     add_input<element_type>("Threshold").set_default(1.0E-16);
//     add_input<size_vector>("Tile Size").set_default(size_vector{180});
//     add_input<element_type>("Screening Threshold").set_default(0.0);
//     add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});
// }

// template<typename T>
// TEMPLATED_MODULE_RUN(Libint, pt::eoctopole<T>) {
//     using element_type   = T;
//     using overlap_pt     = pt::overlap<T>;
//     using edipole_pt     = pt::edipole<T>;
//     using equadrupole_pt = pt::equadrupole<T>;
//     using eoctopole_pt   = pt::eoctopole<T>;
//     using tensor_type    = type::tensor<T>;
//     using value_type     = typename tensor_type::value_type;
//     using type::pair_vector;
//     using type::size_vector;

//     auto [origin, bra_space, ket_space] =
//     eoctopole_pt::unwrap_inputs(inputs);

//     auto& bra         = bra_space.basis_set();
//     auto& ket         = ket_space.basis_set();
//     std::size_t deriv = 0;

//     auto thresh      = inputs.at("Threshold").value<element_type>();
//     auto tile_size   = inputs.at("Tile Size").value<size_vector>();
//     auto cs_thresh   = inputs.at("Screening
//     Threshold").value<element_type>(); auto atom_ranges = inputs.at("Atom
//     Tile Groups").value<pair_vector>(); auto& world      =
//     TA::get_default_world();

//     auto fill =
//       nwx_TA::FillNDFunctor<value_type, libint2::Operator::emultipole3, 2>();
//     fill.initialize(nwx_libint::make_basis_sets({bra, ket}), deriv, thresh,
//                     cs_thresh);
//     fill.factory.origin = origin;

//     auto nopers =
//       libint2::operator_traits<libint2::Operator::emultipole3>::nopers;
//     auto component_range = nwx_TA::make_tiled_range(nopers, nopers);
//     auto trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
//                                         {component_range});

//     auto X = TiledArray::make_array<tensor_type>(world, trange, fill);
//     auto separate_comps = nwx_TA::make_tiled_range(nopers, {1, 3, 6, 10});
//     trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
//                                    {separate_comps});
//     X      = TA::retile(X, trange);

//     // Separate out components
//     tensor_type S, D, Q, O;
//     auto upper      = trange.tiles_range().upbound();
//     using size_type = long;
//     using il_type   = std::initializer_list<size_type>;
//     il_type lo_S{0, 0, 0}, hi_S{1, upper[1], upper[2]};
//     il_type lo_D{1, 0, 0}, hi_D{2, upper[1], upper[2]};
//     il_type lo_Q{2, 0, 0}, hi_Q{3, upper[1], upper[2]};
//     il_type lo_O{3, 0, 0}, hi_O{4, upper[1], upper[2]};

//     S("i,j,k") = X("i,j,k").block(lo_S, hi_S);
//     D("i,j,k") = X("i,j,k").block(lo_D, hi_D);
//     Q("i,j,k") = X("i,j,k").block(lo_Q, hi_Q);
//     O("i,j,k") = X("i,j,k").block(lo_O, hi_O);

//     // Make overlap 2D
//     auto I = TA::diagonal_array<tensor_type, element_type>(
//       world, TA::TiledRange{S.trange().dim(0)});
//     S = libchemist::ta_helpers::einsum::einsum("j,k", "i,j,k", "i", S, I);

//     auto rv = results();
//     rv      = overlap_pt::wrap_results(rv, S);
//     rv      = edipole_pt::wrap_results(rv, D);
//     rv      = equadrupole_pt::wrap_results(rv, Q);
//     rv      = eoctopole_pt::wrap_results(rv, O);
//     return rv;
// }

// template class Libint<pt::edipole<double>>;
// template class Libint<pt::equadrupole<double>>;
// template class Libint<pt::eoctopole<double>>;

// } // namespace integrals
