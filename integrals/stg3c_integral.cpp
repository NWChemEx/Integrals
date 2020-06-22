#include "integrals/stg_integrals.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "integrals/libint_integral.hpp"
#include <property_types/ao_integrals/stg.hpp>
#include <property_types/cauchy_schwarz_approximation.hpp>

namespace integrals {

    template<typename element_type>
    using stg3c_type = property_types::STG3CIntegral<element_type>;
    template<typename element_type>
    using libint_type = property_types::LibIntIntegral<element_type>;
    template<typename element_type>
    using tensor = typename type::tensor<element_type>;
    template<typename element_type>
    using value_type = typename tensor<element_type>::value_type;
    template<typename element_type>
    using cs_approx_type = property_types::CauchySchwarzApprox<element_type>;

    template<typename element_type>
    STG3CInt<element_type>::STG3CInt() : sde::ModuleBase(this) {
        description("Computes 2-center Slater geminal integrals with Libint");
        satisfies_property_type<stg3c_type<element_type>>();
        satisfies_property_type<libint_type<element_type>>();

        add_submodule<cs_approx_type<element_type>>("Cauchy-Schwarz")
                .set_description("Computes the Cauchy-Schwarz Matrix for a pair of basis sets");

    }

    template<typename element_type>
    sde::type::result_map STG3CInt<element_type>::run_(sde::type::input_map inputs,
                                                       sde::type::submodule_map submods) const {
        auto [bra, ket1, ket2, deriv, stg_exponent] = stg3c_type<element_type>::unwrap_inputs(inputs);
        auto [thresh, tile_size, cs_thresh, atom_ranges] = libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        auto fill = nwx_TA::FillNDFunctor<value_type<element_type>, libint2::Operator::stg, 3>();
        fill.initialize(nwx_libint::make_basis_sets({bra, ket1, ket2}), deriv, thresh, cs_thresh);
        fill.factory.stg_exponent = stg_exponent;

        if (cs_thresh > 0.0) {
            auto [cs_mat] = submods.at("Cauchy-Schwarz").run_as<cs_approx_type<element_type>>(ket1, ket2, deriv);
            fill.screen.cs_mat2 = cs_mat;
        }

        TA::TiledRange trange;
        if (atom_ranges.empty()) {
            trange = nwx_TA::make_trange({bra, ket1, ket2}, tile_size);
        } else {
            trange = nwx_TA::make_trange({bra, ket1, ket2}, atom_ranges);
        }

        auto I = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

        auto rv = results();
        return stg3c_type<element_type>::wrap_results(rv, I);
    }

    template class STG3CInt<double>;

} // namespace integrals
