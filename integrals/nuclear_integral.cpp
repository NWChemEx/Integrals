#include "nuclear_integral.hpp"
#include "nwx_libint/nwx_libint.hpp"
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
        auto [thresh, tile_size, cs_thresh, atom_ranges] = libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        std::vector<std::pair<double, std::array<double, 3>>> qs;
        for(const auto& ai : mol)
            qs.emplace_back(static_cast<const double&>(ai.Z()), ai.coords());

        auto fill = nwx_TA::FillNDFunctor<value_type<element_type>, libint2::Operator::nuclear, 2>();
        fill.initialize(nwx_libint::make_basis_sets({bra, ket}), deriv, thresh, cs_thresh);
        fill.factory.qs = qs;

        TA::TiledRange trange;
        if (atom_ranges.empty()) {
            trange = nwx_TA::make_trange({bra, ket}, tile_size);
        } else {
            trange = nwx_TA::make_trange({bra, ket}, atom_ranges);
        }

        auto V = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

        auto rv = results();
        return nuclear_type<element_type>::wrap_results(rv, V);
    }

    template class NuclearInt<double>;

} // namespace integrals
