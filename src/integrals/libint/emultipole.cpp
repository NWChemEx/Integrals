#include "detail_/fill_ND_functor.hpp"
#include "detail_/nwx_TA_utils.hpp"
#include "detail_/nwx_libint.hpp"
#include "emultipole.hpp"
#include <simde/tensor_representation/tensor_representation.hpp>

// TODO: Unify implementations. Maybe with recursion?

namespace integrals {

using size_type   = std::size_t;
using size_vector = std::vector<size_type>;
using pair_vector = std::vector<std::pair<size_type, size_type>>;
using tensor_type = TA::TSpArrayD;
using value_type  = typename tensor_type::value_type;

using identity_op   = simde::type::el_identity;
using dipole_op     = simde::type::el_dipole;
using quadrupole_op = simde::type::el_quadrupole;
using octupole_op   = simde::type::el_octupole;

using overlap_pt    = simde::AOTensorRepresentation<2, identity_op>;
using dipole_pt     = simde::AOTensorRepresentation<2, dipole_op>;
using quadrupole_pt = simde::AOTensorRepresentation<2, quadrupole_op>;
using octupole_pt   = simde::AOTensorRepresentation<2, octupole_op>;

MODULE_CTOR(LibintDipole) {
    description("Computes in-core dipole integrals with libint");
    satisfies_property_type<overlap_pt>();
    satisfies_property_type<dipole_pt>();

    change_input("[I_1]").change(identity_op{});

    add_input<double>("Threshold").set_default(1.0E-16);
    add_input<size_vector>("Tile Size").set_default(size_vector{180});
    add_input<double>("Screening Threshold").set_default(0.0);
    add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});
}

MODULE_RUN(LibintDipole) {
    auto [bra_space, r, ket_space] = dipole_pt::unwrap_inputs(inputs);
    auto thresh                    = inputs.at("Threshold").value<double>();
    auto tile_size   = inputs.at("Tile Size").value<size_vector>();
    auto cs_thresh   = inputs.at("Screening Threshold").value<double>();
    auto atom_ranges = inputs.at("Atom Tile Groups").value<pair_vector>();

    auto& bra   = bra_space.basis_set();
    auto& ket   = ket_space.basis_set();
    auto& world = TA::get_default_world();

    constexpr auto libint_op = libint2::Operator::emultipole1;

    auto fill = nwx_TA::FillNDFunctor<value_type, libint_op, 2>();
    auto bs   = nwx_libint::make_basis_sets({bra, ket});
    fill.initialize(bs, 0, thresh, cs_thresh);

    for(size_type i = 0; i < 3; ++i)
        fill.factory.origin[i] = r.gauge_origin().coord(i);

    auto nopers          = libint2::operator_traits<libint_op>::nopers;
    auto component_range = nwx_TA::make_tiled_range(nopers, nopers);
    auto trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
                                        {component_range});

    auto X = TiledArray::make_array<tensor_type>(world, trange, fill);
    auto separate_comps = nwx_TA::make_tiled_range(nopers, {1, 3});
    trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
                                   {separate_comps});
    X      = TA::retile(X, trange);

    // Separate out components
    auto upper      = trange.tiles_range().upbound();
    using size_type = long;
    using il_type   = std::initializer_list<size_type>;
    tensor_type D;
    il_type lo_D{1, 0, 0}, hi_D{2, upper[1], upper[2]};
    D("i,j,k") = X("i,j,k").block(lo_D, hi_D);

    // tensor_type S;
    // il_type lo_S{0, 0, 0}, hi_S{1, upper[1], upper[2]};
    // S("i,j,k") = X("i,j,k").block(lo_S, hi_S);

    // Make overlap 2D
    // auto I_trange = TA::TiledRange{S.trange().dim(0)};
    // auto I        = TA::diagonal_array<tensor_type, double>(world, I_trange);
    // S("j,k")      = S("i,j,k") * I("i");

    auto rv = results();
    // rv      = overlap_pt::wrap_results(rv, simde::type::tensor(S));
    rv = dipole_pt::wrap_results(rv, simde::type::tensor(D));
    return rv;
}

MODULE_CTOR(LibintQuadrupole) {
    description("Computes an in-core integral with libint");
    satisfies_property_type<overlap_pt>();
    satisfies_property_type<dipole_pt>();
    satisfies_property_type<quadrupole_pt>();

    change_input("[I_1]").change(identity_op{});
    change_input("[r_1]^2").change(quadrupole_op{});
    change_input("[r_1]").change(dipole_op{});

    add_input<double>("Threshold").set_default(1.0E-16);
    add_input<size_vector>("Tile Size").set_default(size_vector{180});
    add_input<double>("Screening Threshold").set_default(0.0);
    add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});
}

MODULE_RUN(LibintQuadrupole) {
    auto [bra_space, r2, ket_space] = quadrupole_pt::unwrap_inputs(inputs);
    auto thresh                     = inputs.at("Threshold").value<double>();
    auto tile_size   = inputs.at("Tile Size").value<size_vector>();
    auto cs_thresh   = inputs.at("Screening Threshold").value<double>();
    auto atom_ranges = inputs.at("Atom Tile Groups").value<pair_vector>();

    auto& bra   = bra_space.basis_set();
    auto& ket   = ket_space.basis_set();
    auto& world = TA::get_default_world();

    constexpr auto libint_op = libint2::Operator::emultipole2;
    auto fill = nwx_TA::FillNDFunctor<value_type, libint_op, 2>();
    auto bs   = nwx_libint::make_basis_sets({bra, ket});
    fill.initialize(bs, 0, thresh, cs_thresh);

    for(size_type i = 0; i < 3; ++i)
        fill.factory.origin[i] = r2.gauge_origin().coord(i);

    auto nopers          = libint2::operator_traits<libint_op>::nopers;
    auto component_range = nwx_TA::make_tiled_range(nopers, nopers);
    auto trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
                                        {component_range});

    auto X = TiledArray::make_array<tensor_type>(world, trange, fill);
    auto separate_comps = nwx_TA::make_tiled_range(nopers, {1, 3, 6});
    trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
                                   {separate_comps});
    X      = TA::retile(X, trange);

    // Separate out components
    tensor_type D, Q;
    auto upper      = trange.tiles_range().upbound();
    using size_type = long;
    using il_type   = std::initializer_list<size_type>;
    il_type lo_D{1, 0, 0}, hi_D{2, upper[1], upper[2]};
    il_type lo_Q{2, 0, 0}, hi_Q{3, upper[1], upper[2]};

    D("i,j,k") = X("i,j,k").block(lo_D, hi_D);
    Q("i,j,k") = X("i,j,k").block(lo_Q, hi_Q);

    // Make overlap 2D
    // tensor_type S;
    // il_type lo_S{0, 0, 0}, hi_S{1, upper[1], upper[2]};
    // S("i,j,k") = X("i,j,k").block(lo_S, hi_S);
    // auto I = TA::diagonal_array<tensor_type, element_type>(
    //   world, TA::TiledRange{S.trange().dim(0)});
    // S = libchemist::ta_helpers::einsum::einsum("j,k", "i,j,k", "i", S, I);

    auto rv = results();
    // rv      = overlap_pt::wrap_results(rv, S);
    rv = dipole_pt::wrap_results(rv, simde::type::tensor(D));
    rv = quadrupole_pt::wrap_results(rv, simde::type::tensor(Q));
    return rv;
}

MODULE_CTOR(LibintOctupole) {
    description("Computes an in-core integral with libint");
    satisfies_property_type<overlap_pt>();
    satisfies_property_type<dipole_pt>();
    satisfies_property_type<quadrupole_pt>();
    satisfies_property_type<octupole_pt>();

    change_input("[I_1]").change(identity_op{});
    change_input("[r_1]^2").change(quadrupole_op{});
    change_input("[r_1]").change(dipole_op{});
    change_input("[r_1]^3").change(octupole_op{});

    add_input<double>("Threshold").set_default(1.0E-16);
    add_input<size_vector>("Tile Size").set_default(size_vector{180});
    add_input<double>("Screening Threshold").set_default(0.0);
    add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});
}

MODULE_RUN(LibintOctupole) {
    auto [bra_space, r3, ket_space] = octupole_pt::unwrap_inputs(inputs);
    auto thresh                     = inputs.at("Threshold").value<double>();
    auto tile_size   = inputs.at("Tile Size").value<size_vector>();
    auto cs_thresh   = inputs.at("Screening Threshold").value<double>();
    auto atom_ranges = inputs.at("Atom Tile Groups").value<pair_vector>();

    auto& bra   = bra_space.basis_set();
    auto& ket   = ket_space.basis_set();
    auto& world = TA::get_default_world();

    constexpr auto libint_op = libint2::Operator::emultipole3;
    auto fill = nwx_TA::FillNDFunctor<value_type, libint_op, 2>();
    auto bs   = nwx_libint::make_basis_sets({bra, ket});
    fill.initialize(bs, 0, thresh, cs_thresh);

    for(size_type i = 0; i < 3; ++i)
        fill.factory.origin[i] = r3.gauge_origin().coord(i);

    auto nopers          = libint2::operator_traits<libint_op>::nopers;
    auto component_range = nwx_TA::make_tiled_range(nopers, nopers);
    auto trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
                                        {component_range});

    auto X = TiledArray::make_array<tensor_type>(world, trange, fill);
    auto separate_comps = nwx_TA::make_tiled_range(nopers, {1, 3, 6, 10});
    trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
                                   {separate_comps});
    X      = TA::retile(X, trange);

    // Separate out components
    tensor_type D, Q, O;
    auto upper      = trange.tiles_range().upbound();
    using size_type = long;
    using il_type   = std::initializer_list<size_type>;

    il_type lo_D{1, 0, 0}, hi_D{2, upper[1], upper[2]};
    il_type lo_Q{2, 0, 0}, hi_Q{3, upper[1], upper[2]};
    il_type lo_O{3, 0, 0}, hi_O{4, upper[1], upper[2]};

    D("i,j,k") = X("i,j,k").block(lo_D, hi_D);
    Q("i,j,k") = X("i,j,k").block(lo_Q, hi_Q);
    O("i,j,k") = X("i,j,k").block(lo_O, hi_O);

    // Make overlap 2D
    // tensor_type S;
    // il_type lo_S{0, 0, 0}, hi_S{1, upper[1], upper[2]};
    // S("i,j,k") = X("i,j,k").block(lo_S, hi_S);
    // auto I = TA::diagonal_array<tensor_type, element_type>(
    //   world, TA::TiledRange{S.trange().dim(0)});
    // S = libchemist::ta_helpers::einsum::einsum("j,k", "i,j,k", "i", S, I);

    auto rv = results();
    // rv      = overlap_pt::wrap_results(rv, S);
    rv = dipole_pt::wrap_results(rv, simde::type::tensor(D));
    rv = quadrupole_pt::wrap_results(rv, simde::type::tensor(Q));
    rv = octupole_pt::wrap_results(rv, simde::type::tensor(O));
    return rv;
}

} // namespace integrals
