#include "integrals/emultipole_integrals.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_libint/nwx_libint_factory.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_ND_functor.hpp"
#include <property_types/ao_integrals/emultipole.hpp>

namespace integrals {

    template<typename element_type>
    using eDipole_type = property_types::EDipoleIntegral<element_type>;
    template<typename element_type>
    using tensor = typename type::tensor<element_type>;

    template<typename element_type>
    EDipoleInt<element_type>::EDipoleInt() : sde::ModuleBase(this) {
        description("Computes dipole integrals with Libint");
        satisfies_property_type<eDipole_type<element_type>>();

        add_input<element_type>("Threshold")
                .set_description("Convergence threshold of integrals")
                .set_default(1.0E-16);

        add_input<std::vector<type::size>>("Tile Size")
                .set_description("Size threshold for tiling tensors by atom blocks")
                .set_default(std::vector<type::size>{180});

        add_input<element_type>("Screening Threshold")
                .set_description("Threshold for Cauchy-Schwarz screening")
                .set_default(0.0);
    }

    template<typename element_type>
    sde::type::result_map EDipoleInt<element_type>::run_(sde::type::input_map inputs,
                                                         sde::type::submodule_map submods) const {
        auto [bra, ket, deriv, origin] = eDipole_type<element_type>::unwrap_inputs(inputs);
        auto thresh = inputs.at("Threshold").value<element_type>();
        auto tile_size = inputs.at("Tile Size").value<std::vector<type::size>>();
        auto cs_thresh = inputs.at("Screening Threshold").value<element_type>();
        auto& world = TA::get_default_world();

        auto fill = nwx_TA::FillNDFunctor<typename tensor<element_type>::value_type,
                                          libint2::Operator::emultipole1, 2>();

        fill.LIBasis_sets = nwx_libint::make_basis_sets({bra, ket});

        fill.factory = nwx_libint::LibintFactory<2, libint2::Operator::emultipole1>();
        fill.factory.max_nprims = nwx_libint::sets_max_nprims(fill.LIBasis_sets);
        fill.factory.max_l = nwx_libint::sets_max_l(fill.LIBasis_sets);
        fill.factory.thresh = thresh;
        fill.factory.deriv = deriv;

        fill.cs_thresh = cs_thresh;

        auto nopers = libint2::operator_traits<libint2::Operator::emultipole1>::nopers;
        auto component_range = nwx_TA::make_tiled_range(nopers, nopers);
        auto trange = nwx_TA::make_trange(fill.LIBasis_sets, tile_size, {component_range});

        auto S = TiledArray::make_array<tensor<element_type>>(world, trange, fill);

        auto rv = results();
        return eDipole_type<element_type>::wrap_results(rv, S);
    }

    template class EDipoleInt<double>;

}
