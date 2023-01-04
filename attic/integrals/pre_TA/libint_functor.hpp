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
#include <chemist/molecule.hpp>
#include <libint2.hpp>
#include <sde/types.hpp>
#include <sde/detail_/sde_any.hpp>

namespace integrals::libint::detail_ {

// Functor that wraps the call to libint in a uniform API, used by PIMPLs
template<std::size_t NBases>
struct LibintFunctor {
    using size_type = std::size_t;

    // The type of the index to a shell block
    using shell_index = std::array<size_type, NBases>;

    // The object that actually computes the integrals, made by LibIntIntegral
    libint2::Engine engine;

    // The basis sets, from left to right, in the integral
    std::array<libint2::BasisSet, NBases> bs;

    // initializes libint
    LibintFunctor() { libint2::initialize(); }
    // finalizes libint
    ~LibintFunctor() { libint2::finalize(); }

    // Takes a set of shell indices returns the shell block
    void operator()(shell_index shells) {
        call_libint_(shells, std::make_index_sequence<NBases>());
    }

private:
    /**
     * @brief The function that operator() dispatches to.
     *
     * We need to unpack the indices given to operator() into distinct arguments
     * and not an array.  That's what this function does via the usual
     * std::index_sequence trick.
     *
     * @tparam Is A variadic parameter pack of integers from [0,NBases) to
     *         expand.
     * @param shells The index of the requested shell block
     * @return An std::vector filled with the requested block per operator
     *         component
     * @throws std::bad_alloc if there is insufficient memory to copy the
     * pointers over.  Strong throw guarantee.
     *
     * @par Data Races:
     * Calls to this function modify the internal state and data races may occur
     * if multiple threads call this function concurrently.
     */
    template<size_type... Is>
    void call_libint_(shell_index shells,
                      std::index_sequence<Is...>) {
        engine.compute((bs[Is][shells[Is]])...);
    }

}; // Class LibIntFunctor

// Factors out the building of a Libint2 engine.
template<libint2::Operator op, std::size_t NBases>
static auto make_engine(
  const chemist::Molecule& molecule,
  const typename LibintFunctor<NBases>::size_type max_prims,
  const typename LibintFunctor<NBases>::size_type max_l, const double thresh,
  const typename LibintFunctor<NBases>::size_type deriv,
  const sde::type::any& op_params = sde::type::any{}) {
    libint2::Engine engine(op, max_prims, max_l, deriv, thresh);
    // Take care of any special set-up
    if constexpr(libint2::rank(op) == 2 && NBases == 2) {
        engine.set(libint2::BraKet::xs_xs);
    }
    if constexpr(libint2::rank(op) == 2 && NBases == 3) {
        engine.set(libint2::BraKet::xs_xx);
    }
    if constexpr(op == libint2::Operator::nuclear) {
      std::vector<std::pair<double, std::array<double, 3>>> qs;
      for(const auto& ai : molecule)
        qs.push_back({static_cast<const double&>(ai.Z()), ai.coords()});
      engine.set_params(qs);
    } else if constexpr (op == libint2::Operator::stg || op == libint2::Operator::yukawa) {
      auto stg_exponent = op_params.has_value() ? sde::detail_::SDEAnyCast<double>(op_params) : 1.0;
      if (std::isnan(stg_exponent)) {
        stg_exponent = 1.0;
      }
      engine.set_params(stg_exponent);
    }
    return engine;
}

} // namespace integrals::libint::detail_
