#include "nuclear_integral.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_libint/nwx_libint_factory.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include "integrals/libint_integral.hpp"
#include <property_types/ao_integrals/nuclear.hpp>

namespace integrals {

    template<typename element_type>
    using nuclear_type = property_types::NuclearIntegral<element_type>;
    template<typename element_type>
    using libint_type = property_types::LibIntIntegral<element_type>;
    template<typename element_type>
    using tensor = typename type::tensor<element_type>;
    template<typename element_type>
    using value_type = typename tensor<element_type>::value_type;

    template<typename element_type>
    NuclearInt<element_type>::NuclearInt() : sde::ModuleBase(this) {
        description("Computes nuclear integrals with Libint");
        satisfies_property_type<nuclear_type<element_type>>();
        satisfies_property_type<libint_type<element_type>>();
    }

    template<typename element_type>
    sde::type::result_map NuclearInt<element_type>::run_(sde::type::input_map inputs,
                                                         sde::type::submodule_map submods) const {
        auto [bra, ket, mol, deriv] = nuclear_type<element_type>::unwrap_inputs(inputs);
        auto [thresh, tile_size, cs_thresh] = libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        auto fill = nwx_TA::FillNDFunctor<value_type<element_type>, libint2::Operator::nuclear, 2>();

        fill.LIBasis_sets = nwx_libint::make_basis_sets({bra, ket});

        fill.factory = nwx_libint::LibintFactory<2, libint2::Operator::nuclear>();
        fill.factory.max_nprims = nwx_libint::sets_max_nprims(fill.LIBasis_sets);
        fill.factory.max_l = nwx_libint::sets_max_l(fill.LIBasis_sets);
        fill.factory.thresh = thresh;
        fill.factory.deriv = deriv;
        fill.factory.mol = mol;

        if (cs_thresh != 0.0) {
            fill.cs_thresh = cs_thresh;
            fill.screen.initialize(fill.LIBasis_sets, fill.factory);
        }

        auto trange = nwx_TA::make_trange(fill.LIBasis_sets, tile_size);

        auto V = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

        auto rv = results();
        return nuclear_type<element_type>::wrap_results(rv, V);
    }

    template class NuclearInt<double>;

} // namespace integrals
