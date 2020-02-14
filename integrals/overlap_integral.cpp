#include "integrals/overlap_integral.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_libint/nwx_libint_factory.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_2D_functor.hpp"
#include "integrals/world.hpp"

namespace integrals {

    template<typename element_type>
    using overlap_type = property_types::OverlapIntegral<element_type>;
    template<typename element_type>
    using tensor = typename integrals::type::tensor<element_type>;

    template<typename element_type>
    OverlapInt<element_type>::OverlapInt() : sde::ModuleBase(this) {
        description("Computes overlap integrals with Libint");
        satisfies_property_type<overlap_type<element_type>>();

        add_input<element_type>("Threshold")
                .set_description("Convergence threshold of integrals")
                .set_default(1.0E-16);

        add_input<std::vector<type::size>>("Tile Size")
                .set_description("Size threshold for tiling tensors by atom blocks")
                .set_default(std::vector<type::size>{180});
    }

    template<typename element_type>
    sde::type::result_map OverlapInt<element_type>::run_(sde::type::input_map inputs,
                                                         sde::type::submodule_map submods) const {
        auto [bra, ket, deriv] = overlap_type<element_type>::unwrap_inputs(inputs);
        auto thresh = inputs.at("Threshold").value<element_type>();
        auto tile_size = inputs.at("Tile Size").value<std::vector<type::size>>();
        auto& world = *pworld; // cf. world.hpp

        auto fill = nwx_TA::Fill2DFunctor<typename tensor<element_type>::value_type, libint2::Operator::overlap>();

        fill.LIBasis_sets = nwx_libint::make_basis_sets({bra, ket});

        fill.factory = nwx_libint::LibintFactory<2, libint2::Operator::overlap>();
        fill.factory.max_nprims = nwx_libint::sets_max_nprims(fill.LIBasis_sets);
        fill.factory.max_l = nwx_libint::sets_max_l(fill.LIBasis_sets);
        fill.factory.thresh = thresh;
        fill.factory.deriv = deriv;

        auto trange = nwx_TA::make_trange(fill.LIBasis_sets, tile_size);

        libint2::initialize();
        auto S = TiledArray::make_array<tensor<element_type>>(world, trange, fill);
        libint2::finalize();

        auto rv = results();
        return overlap_type<element_type>::wrap_results(rv, S);
    }

    template class OverlapInt<double>;

} // namespace integrals