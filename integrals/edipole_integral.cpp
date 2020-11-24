#include "integrals/emultipole_integrals.hpp"
#include "integrals/libint_integral.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include <libchemist/ta_helpers/einsum/einsum.hpp>
#include <property_types/ao_integrals/emultipole.hpp>
#include <property_types/ao_integrals/overlap.hpp>

namespace integrals {

template<typename element_type>
using overlap_type = property_types::OverlapIntegral<element_type>;
template<typename element_type>
using eDipole_type = property_types::EDipoleIntegral<element_type>;
template<typename element_type>
using libint_type = property_types::LibIntIntegral<element_type>;
template<typename element_type>
using tensor = typename type::tensor<element_type>;
template<typename element_type>
using value_type = typename tensor<element_type>::value_type;

template<typename element_type>
EDipoleInt<element_type>::EDipoleInt() : sde::ModuleBase(this) {
    description("Computes dipole integrals with Libint");
    satisfies_property_type<overlap_type<element_type>>();
    satisfies_property_type<eDipole_type<element_type>>();
    satisfies_property_type<libint_type<element_type>>();
}

template<typename element_type>
sde::type::result_map EDipoleInt<element_type>::run_(
  sde::type::input_map inputs, sde::type::submodule_map submods) const {
    auto [bra, ket, deriv, origin] =
      eDipole_type<element_type>::unwrap_inputs(inputs);
    auto [thresh, tile_size, cs_thresh, atom_ranges] =
      libint_type<element_type>::unwrap_inputs(inputs);
    auto& world = TA::get_default_world();

    auto fill = nwx_TA::FillNDFunctor<value_type<element_type>,
                                      libint2::Operator::emultipole1, 2>();
    fill.initialize(nwx_libint::make_basis_sets({bra, ket}), deriv, thresh,
                    cs_thresh);
    fill.factory.origin = origin;

    auto nopers =
      libint2::operator_traits<libint2::Operator::emultipole1>::nopers;
    auto component_range = nwx_TA::make_tiled_range(nopers, nopers);
    auto trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
                                        {component_range});

    auto X = TiledArray::make_array<tensor<element_type>>(world, trange, fill);
    auto separate_comps = nwx_TA::make_tiled_range(nopers, {1, 3});
    trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges,
                                   {separate_comps});
    X      = TA::retile(X, trange);

    // Separate out components
    tensor<element_type> S, D;
    auto upper = trange.tiles_range().upbound();
    using size_type = long;
    using il_type = std::initializer_list<size_type>;
    il_type lo_S{0, 0, 0}, hi_S{1, upper[1], upper[2]};
    il_type lo_D{1, 0, 0}, hi_D{2, upper[1], upper[2]};

    S("i,j,k") = X("i,j,k").block(lo_S, hi_S);
    D("i,j,k") = X("i,j,k").block(lo_D, hi_D);

    // Make overlap 2D
    auto I = TA::diagonal_array<tensor<element_type>, element_type>(
      world, TA::TiledRange{S.trange().dim(0)});
    S = libchemist::ta_helpers::einsum::einsum("j,k", "i,j,k", "i", S, I);

    auto rv = results();
    rv      = overlap_type<element_type>::wrap_results(rv, S);
    rv      = eDipole_type<element_type>::wrap_results(rv, D);
    return rv;
}

template class EDipoleInt<double>;

} // namespace integrals
