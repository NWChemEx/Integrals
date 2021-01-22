#include "integrals/libint_integral.hpp"
#include "integrals/nwx_TA/builder_factory.hpp"
#include "integrals/nwx_TA/fill_ND_functor.hpp"
#include "integrals/nwx_TA/nwx_TA_utils.hpp"
#include "integrals/nwx_direct/eri_direct.hpp"
#include "integrals/nwx_direct/eri_direct_type.hpp"
#include "integrals/nwx_libint/nwx_libint.hpp"
#include <property_types/cauchy_schwarz_approximation.hpp>

namespace integrals {

template<typename ElementType>
using eri3_pt = property_types::ao_integrals::ERI3C<ElementType>;

template<typename ElementType>
using eri4_pt = property_types::ao_integrals::ERI4C<ElementType>;

template<typename element_type>
using libint_type = property_types::LibIntIntegral<element_type>;

template<typename element_type>
using tile_type = TA::Tensor<element_type>;

template<typename element_type>
using cs_approx_type = property_types::CauchySchwarzApprox<element_type>;

namespace detail_ {

template<typename PT>
struct DirectImpl;

// TODO: Condense the implementations down (bare minimum make them work for any
// 3c/4c property type)

template<typename T>
struct DirectImpl<eri3_pt<T>> {
    static auto run(sde::type::input_map inputs,
                    sde::type::submodule_map submods) {
        using base_type = eri3_pt<T>;
        using builder_type =
          nwx_TA::FillNDFunctor<tile_type<element_type>,
                                libint2::Operator::coulomb, 3>;
        using bfactory_type =
          nwx_TA::BuilderFactory<tile_type<element_type>,
                                 libint2::Operator::coulomb, 3>;
        using direct_type = DirectTile<tile_type<element_type>, builder_type>;
        using tensor      = TA::DistArray<direct_type, TA::SparsePolicy>;

        auto [bra_space, ket1_space, ket2_space] =
          base_type::unwrap_inputs(inputs);

        auto& bra         = bra_space.basis_set();
        auto& ket1        = ket1_space.basis_set();
        auto& ket2        = ket2_space.basis_set();
        std::size_t deriv = 0;

        auto [thresh, tile_size, cs_thresh, atom_ranges] =
          libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        auto master = builder_type();
        master.initialize(nwx_libint::make_basis_sets({bra, ket1, ket2}), deriv,
                          thresh, cs_thresh);
        if(cs_thresh > 0.0) {
            auto [cs_mat] =
              submods.at("Cauchy-Schwarz")
                .run_as<cs_approx_type<element_type>>(ket1, ket2, deriv);
            master.screen.cs_mat2 = cs_mat;
        }
        auto bfactory = bfactory_type(master);

        auto trange =
          nwx_TA::select_tiling({bra, ket1, ket2}, tile_size, atom_ranges);
        auto I = tensor(world, trange);

        auto initer = [=](TA::Range& range) {
            return direct_type(range, bfactory(range));
        };

        for(const auto& idx : I) {
            auto range       = trange.make_tile_range(idx.index());
            auto tile_future = world.taskq.add(initer, range);
            I.set(idx.index(), tile_future);
        }

        auto rv = results();
        return base_type::wrap_results(rv, I);
    }
};

template<typename T>
struct DirectImpl<eri4_pt<T>> {
    static auto run(sde::type::input_map inputs,
                    sde::type::submodule_map submods) {
        using base_type = eri4_pt<T>;
        using builder_type =
          nwx_TA::FillNDFunctor<tile_type<element_type>,
                                libint2::Operator::coulomb, 4>;
        using bfactory_type =
          nwx_TA::BuilderFactory<tile_type<element_type>,
                                 libint2::Operator::coulomb, 4>;
        using direct_type = DirectTile<tile_type<element_type>, builder_type>;
        using tensor      = TA::DistArray<direct_type, TA::SparsePolicy>;

        auto [bra1_space, bra2_space, ket1_space, ket2_space] =
          base_type::unwrap_inputs(inputs);

        auto& bra1        = bra1_space.basis_set();
        auto& bra2        = bra2_space.basis_set();
        auto& ket1        = ket1_space.basis_set();
        auto& ket2        = ket2_space.basis_set();
        std::size_t deriv = 0;

        auto [thresh, tile_size, cs_thresh, atom_ranges] =
          libint_type<element_type>::unwrap_inputs(inputs);
        auto& world = TA::get_default_world();

        auto master = builder_type();
        master.initialize(nwx_libint::make_basis_sets({bra1, bra2, ket1, ket2}),
                          deriv, thresh, cs_thresh);
        if(cs_thresh > 0.0) {
            auto [cs_mat1] =
              submods.at("Cauchy-Schwarz")
                .run_as<cs_approx_type<element_type>>(bra1, bra2, deriv);
            auto [cs_mat2] =
              submods.at("Cauchy-Schwarz")
                .run_as<cs_approx_type<element_type>>(ket1, ket2, deriv);
            master.screen.cs_mat1 = cs_mat1;
            master.screen.cs_mat2 = cs_mat2;
        }
        auto bfactory = bfactory_type(master);

        auto trange = nwx_TA::select_tiling({bra1, bra2, ket1, ket2}, tile_size,
                                            atom_ranges);
        auto I      = tensor(world, trange);

        auto initer = [=](TA::Range& range) {
            return direct_type(range, bfactory(range));
        };

        for(const auto& idx : I) {
            auto range       = trange.make_tile_range(idx.index());
            auto tile_future = world.taskq.add(initer, range);
            I.set(idx.index(), tile_future);
        }

        auto rv = results();
        return base_type::wrap_results(rv, I);
    }
};

} // namespace detail_

template<typename BaseType>
TEMPLATED_MODULE_CTOR(LibintDirect, BaseType) {
    using element_type = double;

    description("Builds and forgets an integral");
    satisfies_property_type<BaseType>();
    satisfies_property_type<libint_type<element_type>>;

    add_submodule<cs_approx_type<element_type>>("Cauchy-Schwarz")
      .set_description(
        "Computes the Cauchy-Schwarz Matrix for a pair of basis sets");
}

template<typename BaseType>
TEMPLATED_MODULE_RUN(LibintDirect, BaseType) {
    return detail_::DirectImpl<BaseType>::run(inputs, submods);
}

template class LibintDirect<property_types::ao_integrals::ERI3C<double>>;
template class LibintDirect<property_types::ao_integrals::ERI4C<double>>;

} // namespace integrals