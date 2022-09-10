/*
 * Copyright 2022 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "integrals/libint_integral.hpp"
#include "integrals/nwx_TA/builder_factory.hpp"
#include "integrals/nwx_TA/fill_ND_functor.hpp"
#include "integrals/nwx_TA/nwx_TA_utils.hpp"
#include "integrals/nwx_direct/yukawa_direct.hpp"
#include "integrals/nwx_direct/yukawa_direct_type.hpp"
#include "integrals/nwx_libint/nwx_libint.hpp"
#include <property_types/cauchy_schwarz_approximation.hpp>

namespace integrals {

template<typename element_type>
using eri3c_type = property_types::Yukawa3CDirect<element_type>;
template<typename element_type>
using eri4c_type = property_types::Yukawa4CDirect<element_type>;
template<typename element_type>
using libint_type = property_types::LibIntIntegral<element_type>;
template<typename element_type>
using tile_type = TA::Tensor<element_type>;
template<typename element_type>
using cs_approx_type = property_types::CauchySchwarzApprox<element_type>;

template<typename element_type>
Yukawa3CIntDirect<element_type>::Yukawa3CIntDirect() : sde::ModuleBase(this) {
    description("Computes 3-center electron repulsion integrals with Libint");
    satisfies_property_type<eri3c_type<element_type>>();
    satisfies_property_type<libint_type<element_type>>();

    add_submodule<cs_approx_type<element_type>>("Cauchy-Schwarz")
      .set_description(
        "Computes the Cauchy-Schwarz Matrix for a pair of basis sets");
}

template<typename element_type>
sde::type::result_map Yukawa3CIntDirect<element_type>::run_(
  sde::type::input_map inputs, sde::type::submodule_map submods) const {
    using builder_type  = nwx_TA::FillNDFunctor<tile_type<element_type>,
                                               libint2::Operator::yukawa, 3>;
    using bfactory_type = nwx_TA::BuilderFactory<tile_type<element_type>,
                                                 libint2::Operator::yukawa, 3>;
    using direct_type   = DirectTile<tile_type<element_type>, builder_type>;
    using tensor        = TA::DistArray<direct_type, TA::SparsePolicy>;

    auto [bra, ket1, ket2, deriv, stg_exponent] =
      eri3c_type<element_type>::unwrap_inputs(inputs);
    auto [thresh, tile_size, cs_thresh, atom_ranges] =
      libint_type<element_type>::unwrap_inputs(inputs);
    auto& world = TA::get_default_world();

    auto master = builder_type();
    master.initialize(nwx_libint::make_basis_sets({bra, ket1, ket2}), deriv,
                      thresh, cs_thresh);
    master.factory.stg_exponent = stg_exponent;
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
    return eri3c_type<element_type>::wrap_results(rv, I);
}

template<typename element_type>
Yukawa4CIntDirect<element_type>::Yukawa4CIntDirect() : sde::ModuleBase(this) {
    description("Computes 4-center electron repulsion integrals with Libint");
    satisfies_property_type<eri4c_type<element_type>>();
    satisfies_property_type<libint_type<element_type>>();

    add_submodule<cs_approx_type<element_type>>("Cauchy-Schwarz")
      .set_description(
        "Computes the Cauchy-Schwarz Matrix for a pair of basis sets");
}

template<typename element_type>
sde::type::result_map Yukawa4CIntDirect<element_type>::run_(
  sde::type::input_map inputs, sde::type::submodule_map submods) const {
    using builder_type  = nwx_TA::FillNDFunctor<tile_type<element_type>,
                                               libint2::Operator::yukawa, 4>;
    using bfactory_type = nwx_TA::BuilderFactory<tile_type<element_type>,
                                                 libint2::Operator::yukawa, 4>;
    using direct_type   = DirectTile<tile_type<element_type>, builder_type>;
    using tensor        = TA::DistArray<direct_type, TA::SparsePolicy>;

    auto [bra1, bra2, ket1, ket2, deriv, stg_exponent] =
      eri4c_type<element_type>::unwrap_inputs(inputs);
    auto [thresh, tile_size, cs_thresh, atom_ranges] =
      libint_type<element_type>::unwrap_inputs(inputs);
    auto& world = TA::get_default_world();

    auto master = builder_type();
    master.initialize(nwx_libint::make_basis_sets({bra1, bra2, ket1, ket2}),
                      deriv, thresh, cs_thresh);
    master.factory.stg_exponent = stg_exponent;
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

    auto trange =
      nwx_TA::select_tiling({bra1, bra2, ket1, ket2}, tile_size, atom_ranges);
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
    return eri4c_type<element_type>::wrap_results(rv, I);
}

template class Yukawa3CIntDirect<double>;
template class Yukawa4CIntDirect<double>;

} // namespace integrals