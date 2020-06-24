#include "integrals/stg_integrals.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "integrals/libint_integral.hpp"
#include <property_types/ao_integrals/stg.hpp>

namespace integrals {

    template<typename element_type>
    using stg2c_type = property_types::STG2CIntegral<element_type>;
    template<typename element_type>
    using libint_type = property_types::LibIntIntegral<element_type>;
    template<typename element_type>
    using tensor = typename type::tensor<element_type>;
    template<typename element_type>
    using value_type = typename tensor<element_type>::value_type;

    template<typename element_type>
    STG2CInt<element_type>::STG2CInt() : sde::ModuleBase(this) {
        description("Computes 2-center Slater geminal integrals with Libint");
        satisfies_property_type<stg2c_type<element_type>>();
        satisfies_property_type<libint_type<element_type>>();
    }

    template<typename element_type>
    sde::type::result_map STG2CInt<element_type>::run_(sde::type::input_map inputs,
                                                       sde::type::submodule_map submods) const {
        auto [bra, ket, deriv, stg_exponent] = stg2c_type<element_type>::unwrap_inputs(inputs);
        auto [thresh, tile_size, cs_thresh, atom_ranges] = libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        auto fill = nwx_TA::FillNDFunctor<value_type<element_type>, libint2::Operator::stg, 2>();
        fill.initialize(nwx_libint::make_basis_sets({bra, ket}), deriv, thresh, cs_thresh);
        fill.factory.stg_exponent = stg_exponent;

        auto trange = nwx_TA::select_tiling({bra, ket}, tile_size, atom_ranges);

        auto I = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

        auto rv = results();
        return stg2c_type<element_type>::wrap_results(rv, I);
    }

    template class STG2CInt<double>;

} // namespace integrals
