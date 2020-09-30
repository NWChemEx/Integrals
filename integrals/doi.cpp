#include "integrals/doi.hpp"
#include "integrals/libint_integral.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include <property_types/ao_integrals/doi.hpp>

namespace integrals {

template<typename element_type>
using doi_type = property_types::DOI<element_type>;
template<typename element_type>
using libint_type = property_types::LibIntIntegral<element_type>;
template<typename element_type>
using tensor = typename type::tensor<element_type>;
template<typename element_type>
using value_type = typename tensor<element_type>::value_type;

template<typename element_type>
DOInt<element_type>::DOInt() : sde::ModuleBase(this) {
    description("Computes differential overlap integrals with Libint");
    satisfies_property_type<doi_type<element_type>>();
    satisfies_property_type<libint_type<element_type>>();
}

template<typename element_type>
sde::type::result_map DOInt<element_type>::run_(
  sde::type::input_map inputs, sde::type::submodule_map submods) const {
    auto [bra, ket, deriv] = doi_type<element_type>::unwrap_inputs(inputs);
    auto [thresh, tile_size, cs_thresh, atom_ranges] =
      libint_type<element_type>::unwrap_inputs(inputs);
    auto& world = TA::get_default_world();

    auto fill = nwx_TA::FillNDFunctor<value_type<element_type>,
                                      libint2::Operator::delta, 4>();
    fill.initialize(nwx_libint::make_basis_sets({bra, bra, ket, ket}), deriv,
                    thresh, cs_thresh);

    auto trange =
      nwx_TA::select_tiling({bra, bra, ket, ket}, tile_size, atom_ranges);

    auto I = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

    auto rv = results();
    return doi_type<element_type>::wrap_results(rv, I);
}

template class DOInt<double>;

} // namespace integrals