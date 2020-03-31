#include "integrals/yukawa_integrals.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "integrals/libint_integral.hpp"
#include <property_types/ao_integrals/yukawa.hpp>

namespace integrals {

    template<typename element_type>
    using yukawa3c_type = property_types::Yukawa3CIntegral<element_type>;
    template<typename element_type>
    using libint_type = property_types::LibIntIntegral<element_type>;
    template<typename element_type>
    using tensor = typename type::tensor<element_type>;
    template<typename element_type>
    using value_type = typename tensor<element_type>::value_type;

    template<typename element_type>
    Yukawa3CInt<element_type>::Yukawa3CInt() : sde::ModuleBase(this) {
        description("Computes 2-center Slater geminal integrals with Libint");
        satisfies_property_type<yukawa3c_type<element_type>>();
        satisfies_property_type<libint_type<element_type>>();
    }

    template<typename element_type>
    sde::type::result_map Yukawa3CInt<element_type>::run_(sde::type::input_map inputs,
                                                          sde::type::submodule_map submods) const {
        auto [bra, ket1, ket2, deriv, stg_exponent] = yukawa3c_type<element_type>::unwrap_inputs(inputs);
        auto [thresh, tile_size, cs_thresh] = libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        auto fill = nwx_TA::FillNDFunctor<value_type<element_type>, libint2::Operator::yukawa, 3>();
        fill.initialize(nwx_libint::make_basis_sets({bra, ket1, ket2}), deriv, thresh, cs_thresh);
        fill.factory.stg_exponent = stg_exponent;

        auto trange = nwx_TA::make_trange({bra, ket1, ket2}, tile_size);

        auto I = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

        auto rv = results();
        return yukawa3c_type<element_type>::wrap_results(rv, I);
    }

    template class Yukawa3CInt<double>;

} // namespace integrals
