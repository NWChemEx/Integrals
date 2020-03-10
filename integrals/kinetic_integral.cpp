#include "kinetic_integral.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_libint/nwx_libint_factory.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "integrals/libint_integral.hpp"
#include <property_types/ao_integrals/kinetic.hpp>

namespace integrals {

    template<typename element_type>
    using kinetic_type = property_types::KineticIntegral<element_type>;
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
    sde::type::result_map KineticInt<element_type>::run_(sde::type::input_map inputs,
                                                         sde::type::submodule_map submods) const {
        auto [bra, ket, deriv] = kinetic_type<element_type>::unwrap_inputs(inputs);
        auto [thresh, tile_size, cs_thresh] = libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        auto fill = nwx_TA::FillNDFunctor<value_type<element_type>, libint2::Operator::kinetic, 2>();

        fill.LIBasis_sets = nwx_libint::make_basis_sets({bra, ket});

        fill.factory = nwx_libint::LibintFactory<2, libint2::Operator::kinetic>();
        fill.factory.max_nprims = nwx_libint::sets_max_nprims(fill.LIBasis_sets);
        fill.factory.max_l = nwx_libint::sets_max_l(fill.LIBasis_sets);
        fill.factory.thresh = thresh;
        fill.factory.deriv = deriv;

        if (cs_thresh != 0.0) {
            fill.cs_thresh = cs_thresh;
            fill.screen.initialize(fill.LIBasis_sets, fill.factory);
        }

        auto trange = nwx_TA::make_trange(fill.LIBasis_sets, tile_size);

        auto T = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

        auto rv = results();
        return kinetic_type<element_type>::wrap_results(rv, T);
    }

    template class KineticInt<double>;

} // namespace integrals
