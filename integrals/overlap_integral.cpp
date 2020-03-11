#include "integrals/overlap_integral.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "integrals/libint_integral.hpp"
#include <property_types/ao_integrals/overlap.hpp>

namespace integrals {

    template<typename element_type>
    using overlap_type = property_types::OverlapIntegral<element_type>;
    template<typename element_type>
    using libint_type = property_types::LibIntIntegral<element_type>;
    template<typename element_type>
    using tensor = typename type::tensor<element_type>;
    template<typename element_type>
    using value_type = typename tensor<element_type>::value_type;

    template<typename element_type>
    OverlapInt<element_type>::OverlapInt() : sde::ModuleBase(this) {
        description("Computes overlap integrals with Libint");
        satisfies_property_type<overlap_type<element_type>>();
        satisfies_property_type<libint_type<element_type>>();
    }

    template<typename element_type>
    sde::type::result_map OverlapInt<element_type>::run_(sde::type::input_map inputs,
                                                         sde::type::submodule_map submods) const {
        auto [bra, ket, deriv] = overlap_type<element_type>::unwrap_inputs(inputs);
        auto [thresh, tile_size, cs_thresh] = libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        auto fill = nwx_TA::FillNDFunctor<value_type<element_type>, libint2::Operator::overlap, 2>();
        fill.initialize(nwx_libint::make_basis_sets({bra, ket}), deriv, thresh, cs_thresh);

        auto trange = nwx_TA::make_trange(fill.LIBasis_sets, tile_size);

        auto S = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

        auto rv = results();
        return overlap_type<element_type>::wrap_results(rv, S);
    }

    template class OverlapInt<double>;

} // namespace integrals