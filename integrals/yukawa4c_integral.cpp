#include "integrals/yukawa_integrals.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_libint/nwx_libint_factory.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_4D_functor.hpp"
#include "integrals/world.hpp"

namespace integrals {

    template<typename element_type>
    Yukawa4CInt<element_type>::Yukawa4CInt() : sde::ModuleBase(this) {
        description("Computes 2-center Slater geminal integrals with Libint");
        satisfies_property_type<yukawa4c_type>();

        add_input<element_type>("Threshold")
                .set_description("Convergence threshold of integrals")
                .set_default(1.0E-16);

        add_input<size_vec>("Tile Size")
                .set_description("Size threshold for tiling tensors by atom blocks")
                .set_default(size_vec{180});
    }

    template<typename element_type>
    sde::type::result_map Yukawa4CInt<element_type>::run_(sde::type::input_map inputs,
                                                          sde::type::submodule_map submods) const {
        auto [bra1, bra2, ket1, ket2, deriv, stg_exponent] = yukawa4c_type::unwrap_inputs(inputs);
        auto thresh = inputs.at("Threshold").value<element_type>();
        auto tile_size = inputs.at("Tile Size").value<size_vec>();
        auto& world = *pworld; // cf. world.hpp

        auto fill = nwx_TA::Fill4DFunctor<typename tensor::value_type, libint2::Operator::yukawa>();

        fill.LIBasis_sets = nwx_libint::make_basis_sets({bra1, bra2, ket1, ket2});

        fill.factory = nwx_libint::LibintFactory<4, libint2::Operator::yukawa>();
        fill.factory.max_nprims = nwx_libint::sets_max_nprims(fill.LIBasis_sets);
        fill.factory.max_l = nwx_libint::sets_max_l(fill.LIBasis_sets);
        fill.factory.thresh = thresh;
        fill.factory.deriv = deriv;
        fill.factory.stg_exponent = stg_exponent;

        auto trange = nwx_TA::make_trange(fill.LIBasis_sets, tile_size);

        libint2::initialize();
        auto I = TiledArray::make_array<tensor>(world, trange, fill);
        libint2::finalize();

        auto rv = results();
        return yukawa4c_type::wrap_results(rv, I);
    }

    template class Yukawa4CInt<double>;

} // namespace integrals
