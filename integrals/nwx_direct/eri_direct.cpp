#include "integrals/nwx_direct/eri_direct.hpp"
#include "integrals/nwx_libint/nwx_libint.hpp"
#include "integrals/nwx_TA/nwx_TA_utils.hpp"
#include "integrals/nwx_TA/fill_ND_functor.hpp"
#include "integrals/libint_integral.hpp"
#include "integrals/nwx_direct/eri_direct_type.hpp"
#include "integrals/nwx_TA/builder_factory.hpp"

namespace integrals {

    template<typename element_type>
    using eri3c_type = property_types::ERI3CDirect<element_type>;
    template<typename element_type>
    using eri4c_type = property_types::ERI4CDirect<element_type>;
    template<typename element_type>
    using libint_type = property_types::LibIntIntegral<element_type>;
    template<typename element_type>
    using tile_type = TA::Tensor<element_type>;

    template<typename element_type>
    ERI3CIntDirect<element_type>::ERI3CIntDirect() : sde::ModuleBase(this) {
        description("Computes 3-center electron repulsion integrals with Libint");
        satisfies_property_type<eri3c_type<element_type>>();
        satisfies_property_type<libint_type<element_type>>();
    }

    template<typename element_type>
    sde::type::result_map ERI3CIntDirect<element_type>::run_(sde::type::input_map inputs,
                                                             sde::type::submodule_map submods) const {
        using builder_type = nwx_TA::FillNDFunctor<tile_type<element_type>, libint2::Operator::coulomb, 3>;
        using bfactory_type = nwx_TA::BuilderFactory<tile_type<element_type>, libint2::Operator::coulomb, 3>;
        using direct_type = DirectTile<tile_type<element_type>, builder_type>;
        using tensor = TA::DistArray<direct_type, TA::SparsePolicy>;

        auto [bra, ket1, ket2, deriv] = eri3c_type<element_type>::unwrap_inputs(inputs);
        auto [thresh, tile_size, cs_thresh] = libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        auto master = builder_type();
        master.initialize(nwx_libint::make_basis_sets({bra, ket1, ket2}), deriv, thresh, cs_thresh);
        auto bfactory = bfactory_type(master);

        auto trange = nwx_TA::make_trange(master.LIBasis_sets, tile_size);
        auto I = tensor(world, trange);

        auto initer = [=](TA::Range& range) { return direct_type(range, bfactory(range)); };

        for (const auto& idx : I) {
            auto range = trange.make_tile_range(idx.index());
            auto tile_future = world.taskq.add(initer, range);
            I.set(idx.index(), tile_future);
        }

        auto rv = results();
        return eri3c_type<element_type>::wrap_results(rv, I);
    }

    template<typename element_type>
    ERI4CIntDirect<element_type>::ERI4CIntDirect() : sde::ModuleBase(this) {
        description("Computes 4-center electron repulsion integrals with Libint");
        satisfies_property_type<eri4c_type<element_type>>();
        satisfies_property_type<libint_type<element_type>>();
    }

    template<typename element_type>
    sde::type::result_map ERI4CIntDirect<element_type>::run_(sde::type::input_map inputs,
                                                             sde::type::submodule_map submods) const {
        using builder_type = nwx_TA::FillNDFunctor<tile_type<element_type>, libint2::Operator::coulomb, 4>;
        using bfactory_type = nwx_TA::BuilderFactory<tile_type<element_type>, libint2::Operator::coulomb, 4>;
        using direct_type = DirectTile<tile_type<element_type>, builder_type>;
        using tensor = TA::DistArray<direct_type, TA::SparsePolicy>;

        auto [bra1, bra2, ket1, ket2, deriv] = eri4c_type<element_type>::unwrap_inputs(inputs);
        auto [thresh, tile_size, cs_thresh] = libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        auto master = builder_type();
        master.initialize(nwx_libint::make_basis_sets({bra1, bra2, ket1, ket2}), deriv, thresh, cs_thresh);
        auto bfactory = bfactory_type(master);

        auto trange = nwx_TA::make_trange(master.LIBasis_sets, tile_size);
        auto I = tensor(world, trange);

        auto initer = [=](TA::Range& range) { return direct_type(range, bfactory(range)); };

        for (const auto& idx : I) {
            auto range = trange.make_tile_range(idx.index());
            auto tile_future = world.taskq.add(initer, range);
            I.set(idx.index(), tile_future);
        }

        auto rv = results();
        return eri4c_type<element_type>::wrap_results(rv, I);
    }

    template class ERI3CIntDirect<double>;
    template class ERI4CIntDirect<double>;

} // namespace integrals