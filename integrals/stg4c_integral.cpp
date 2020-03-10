#include "integrals/stg_integrals.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_libint/nwx_libint_factory.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "integrals/libint_integral.hpp"
#include <property_types/ao_integrals/stg.hpp>

namespace integrals {

    template<typename element_type>
    using stg4c_type = property_types::STG4CIntegral<element_type>;
    template<typename element_type>
    using libint_type = property_types::LibIntIntegral<element_type>;
    template<typename element_type>
    using tensor = typename type::tensor<element_type>;
    template<typename element_type>
    using value_type = typename tensor<element_type>::value_type;

    template<typename element_type>
    STG4CInt<element_type>::STG4CInt() : sde::ModuleBase(this) {
        description("Computes 2-center Slater geminal integrals with Libint");
        satisfies_property_type<stg4c_type<element_type>>();
        satisfies_property_type<libint_type<element_type>>();
    }

    template<typename element_type>
    sde::type::result_map STG4CInt<element_type>::run_(sde::type::input_map inputs,
                                                       sde::type::submodule_map submods) const {
        auto [bra1, bra2, ket1, ket2, deriv, stg_exponent] = stg4c_type<element_type>::unwrap_inputs(inputs);
        auto [thresh, tile_size, cs_thresh] = libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        auto fill = nwx_TA::FillNDFunctor<value_type<element_type>, libint2::Operator::stg, 4>();

        fill.LIBasis_sets = nwx_libint::make_basis_sets({bra1, bra2, ket1, ket2});

        fill.factory = nwx_libint::LibintFactory<4, libint2::Operator::stg>();
        fill.factory.max_nprims = nwx_libint::sets_max_nprims(fill.LIBasis_sets);
        fill.factory.max_l = nwx_libint::sets_max_l(fill.LIBasis_sets);
        fill.factory.thresh = thresh;
        fill.factory.deriv = deriv;
        fill.factory.stg_exponent = stg_exponent;

        if (cs_thresh != 0.0) {
            fill.cs_thresh = cs_thresh;
            fill.screen.initialize(fill.LIBasis_sets, fill.factory);
        }

        auto trange = nwx_TA::make_trange(fill.LIBasis_sets, tile_size);

        auto I = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

        auto rv = results();
        return stg4c_type<element_type>::wrap_results(rv, I);
    }

    template class STG4CInt<double>;

} // namespace integrals
