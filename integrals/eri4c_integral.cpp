#include "integrals/eri_integrals.hpp"
#include "nwx_libint/nwx_libint.hpp"
#include "nwx_libint/nwx_libint_factory.hpp"
#include "nwx_TA/nwx_TA_utils.hpp"
#include "nwx_TA/fill_4D_functor.hpp"
#include "integrals/world.hpp"

namespace integrals {

    template<typename element_type>
    ERI4CInt<element_type>::ERI4CInt() : sde::ModuleBase(this) {
        description("Computes 4-center electron repulsion integrals with Libint");
        satisfies_property_type<eri4c_type>();

        add_input<element_type>("Threshold")
                .set_description("Convergence threshold of integrals")
                .set_default(1.0E-16);

        add_input<size_type>("Tile Size")
                .set_description("Size threshold for tiling tensors by atom blocks")
                .set_default(size_type{180});
    }

    template<typename element_type>
    sde::type::result_map ERI4CInt<element_type>::run_(sde::type::input_map inputs,
                                                       sde::type::submodule_map submods) const {
        auto [bra1, bra2, ket1, ket2, deriv] = eri4c_type::unwrap_inputs(inputs);
        auto thresh = inputs.at("Threshold").value<element_type>();
        auto tile_size = inputs.at("Tile Size").value<size_type>();
        auto& world = *pworld; // cf. world.hpp

        std::vector<basis_set> basis_sets{bra1, bra2, ket1, ket2};

        std::vector<TA::TiledRange1> ranges{};
        std::vector<libint2::BasisSet> LIBasis_sets{};
        size_type max_nprim = 0;
        int max_l = 0;

        for (auto i = 0; i < basis_sets.size(); ++i) {
            LIBasis_sets.push_back(nwx_libint::make_basis(basis_sets[i]));
            ranges.push_back(nwx_TA::make_tiled_range(LIBasis_sets[i], tile_size));

            auto max_nprim_i = libint2::max_nprim(LIBasis_sets[i]);
            auto max_l_i = libint2::max_l(LIBasis_sets[i]);
            max_nprim = std::max(max_nprim, max_nprim_i);
            max_l = std::max(max_l, max_l_i);
        }

        TA::TiledRange trange(ranges.begin(), ranges.end());

        libint2::initialize();

        nwx_libint::LibintFactory<4, libint2::Operator::coulomb> factory(max_nprim, max_l, thresh, deriv);
        nwx_TA::Fill4DFunctor<typename tensor::value_type, libint2::Operator::coulomb> fill(LIBasis_sets, factory);

        auto I = TiledArray::make_array<TiledArray::TSpArrayD>(world, trange, fill);

        libint2::finalize();

        auto rv = results();
        return eri4c_type::wrap_results(rv, I);
    }

    template class ERI4CInt<double>;

} // namespace integrals
