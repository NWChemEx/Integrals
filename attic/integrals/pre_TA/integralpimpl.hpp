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

#pragma once
#include "integrals/libint_functor.hpp"
#include "integrals/libint_integral.hpp"

namespace integrals::libint::detail_ {
// Defines the API for the LibintIntegral PIMPL (move to header if needed)
template<libint2::Operator op, std::size_t NBases, typename element_type>
struct IntegralPIMPL {
    using size_type = std::size_t;
    // For integrals with multiple components
    constexpr static size_type extra =
      (libint2::operator_traits<op>::nopers > 1) ? 1 : 0;
    // Typedef of the main class
    using main_type        = Integral<op, NBases, element_type>;
    using tensor_type      = typename main_type::tensor_type;
    using basis_array_type = typename main_type::basis_array_type;
    using tiled_AO         = std::vector<tamm::TiledIndexSpace>;
    using fxn_type         = LibintFunctor<NBases>;

    // Public API to PIMPL
    virtual tensor_type run_impl(
      const tiled_AO& tAOs,
      const std::array<std::vector<size_type>, NBases>& atom_blocks,
      const basis_array_type& bases, fxn_type&& fxn,
      const element_type schwarz_thresh) {
        return run_impl_(tAOs, atom_blocks, bases, std::move(fxn),
                         schwarz_thresh);
    }

private:
    // Implemented by derived class
    virtual tensor_type run_impl_(
      const tiled_AO& tAOs,
      const std::array<std::vector<size_type>, NBases>& atom_blocks,
      const basis_array_type& bases, fxn_type&& fxn,
      const element_type schwarz_thresh) = 0;
};

} // namespace integrals::libint::detail_
