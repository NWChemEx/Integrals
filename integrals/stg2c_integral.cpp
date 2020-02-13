#include "integrals/stg_integrals.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_libint/nwx_libint_factory.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_2D_functor.hpp"
#include "integrals/world.hpp"

namespace integrals {

    template<typename element_type>
    STG2CInt<element_type>::STG2CInt() : sde::ModuleBase(this) {
        description("Computes 2-center Slater geminal integrals with Libint");
        satisfies_property_type<stg2c_type>();

        add_input<element_type>("Threshold")
                .set_description("Convergence threshold of integrals")
                .set_default(1.0E-16);

        add_input<size_vec>("Tile Size")
                .set_description("Size threshold for tiling tensors by atom blocks")
                .set_default(size_vec{180});
    }

    template<typename element_type>
    sde::type::result_map STG2CInt<element_type>::run_(sde::type::input_map inputs,
                                                       sde::type::submodule_map submods) const {
        auto [bra, ket, deriv, stg_exponent] = stg2c_type::unwrap_inputs(inputs);
        auto thresh = inputs.at("Threshold").value<element_type>();
        auto tile_size = inputs.at("Tile Size").value<size_vec>();
        auto& world = *pworld; // cf. world.hpp

        std::vector<basis_set> basis_sets{bra, ket};

        std::vector<libint2::BasisSet> LIBasis_sets{};
        size_type max_nprim = 0;
        int max_l = 0;

        for (auto i = 0; i < basis_sets.size(); ++i) {
            LIBasis_sets.push_back(nwx_libint::make_basis(basis_sets[i]));

            auto max_nprim_i = libint2::max_nprim(LIBasis_sets[i]);
            auto max_l_i = libint2::max_l(LIBasis_sets[i]);
            max_nprim = std::max(max_nprim, max_nprim_i);
            max_l = std::max(max_l, max_l_i);
        }

        auto trange = nwx_TA::make_trange(LIBasis_sets, tile_size);

        libint2::initialize();

        nwx_libint::LibintFactory<2, libint2::Operator::stg> factory(max_nprim, max_l, thresh, deriv);
        factory.stg_exponent = stg_exponent;
        nwx_TA::Fill2DFunctor<typename tensor::value_type, libint2::Operator::stg> fill(LIBasis_sets, factory);

        auto I = TiledArray::make_array<tensor>(world, trange, fill);

        libint2::finalize();

        auto rv = results();
        return stg2c_type::wrap_results(rv, I);
    }

    template class STG2CInt<double>;

} // namespace integrals