#include "integrals/libint_integral.hpp"
#include "kinetic_integral.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include <property_types/ao_integrals/kinetic.hpp>

namespace integrals {

template<typename element_type>
using kinetic_type = property_types::ao_integrals::Kinetic<element_type>;
template<typename element_type>
using libint_type = property_types::LibIntIntegral<element_type>;
template<typename element_type>
using tensor = typename type::tensor<element_type>;
template<typename element_type>
using value_type = typename tensor<element_type>::value_type;

template<typename element_type>
KineticInt<element_type>::KineticInt() : sde::ModuleBase(this) {
    description("Computes kinetic integrals with Libint");
    satisfies_property_type<kinetic_type<element_type>>();
    satisfies_property_type<libint_type<element_type>>();
}

template<typename element_type>
sde::type::result_map KineticInt<element_type>::run_(
  sde::type::input_map inputs, sde::type::submodule_map submods) const {
    auto [bra_space, ket_space] =
      kinetic_type<element_type>::unwrap_inputs(inputs);

    auto& bra         = bra_space.basis_set();
    auto& ket         = ket_space.basis_set();
    std::size_t deriv = 0;

    auto [thresh, tile_size, cs_thresh, atom_ranges] =
      libint_type<element_type>::unwrap_inputs(inputs);
    auto& world = TA::get_default_world();

    auto fill = nwx_TA::FillNDFunctor<value_type<element_type>,
                                      libint2::Operator::kinetic, 2>();
    fill.initialize(nwx_libint::make_basis_sets({bra, ket}), deriv, thresh,
                    cs_thresh);

    auto trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges);

    auto T = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

    auto rv = results();
    return kinetic_type<element_type>::wrap_results(rv, T);
}

template class KineticInt<double>;

} // namespace integrals