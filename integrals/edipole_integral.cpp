#include "integrals/emultipole_integrals.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "integrals/libint_integral.hpp"
#include <property_types/ao_integrals/emultipole.hpp>

namespace integrals {

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
        satisfies_property_type<eDipole_type<element_type>>();
        satisfies_property_type<libint_type<element_type>>();
    }

    template<typename element_type>
    sde::type::result_map EDipoleInt<element_type>::run_(sde::type::input_map inputs,
                                                         sde::type::submodule_map submods) const {
        auto [bra, ket, deriv, origin] = eDipole_type<element_type>::unwrap_inputs(inputs);
        auto [thresh, tile_size, cs_thresh, atom_ranges] = libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        auto fill = nwx_TA::FillNDFunctor<value_type<element_type>, libint2::Operator::emultipole1, 2>();
        fill.initialize(nwx_libint::make_basis_sets({bra, ket}), deriv, thresh, cs_thresh);
        fill.factory.origin = origin;

        auto nopers = libint2::operator_traits<libint2::Operator::emultipole1>::nopers;
        auto component_range = nwx_TA::make_tiled_range(nopers, nopers);
        TA::TiledRange trange;
        if (atom_ranges.empty()) {
            trange = nwx_TA::make_trange({bra, ket}, tile_size, {component_range});
        } else {
            trange = nwx_TA::make_trange({bra, ket}, atom_ranges, {component_range});
        }

        auto S = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

        auto rv = results();
        return eDipole_type<element_type>::wrap_results(rv, S);
    }

    template class EDipoleInt<double>;

}
