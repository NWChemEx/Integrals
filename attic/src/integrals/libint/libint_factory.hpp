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

#include "detail_/make_engine.hpp"
#include <chemist/integral_factory/integral_factory_pimpl.hpp>
#include <libint2.hpp>

namespace integrals::libint {

template<std::size_t NCenters, typename OpType>
class LibintFactory : public chemist::IntegralFactoryPIMPL {
public:
    /// Base type and this type
    using base_type = chemist::IntegralFactoryPIMPL;
    using my_type   = LibintFactory<NCenters, OpType>;

    /// Types from the base type
    using pimpl_ptr_type = typename base_type::pimpl_ptr_t;
    using indices_type   = typename base_type::indices_t;
    using buffer_type    = typename base_type::buffer_t;

    /// Types unique to this class
    using libint_basis_type       = libint2::BasisSet;
    using libint_basis_vector     = std::vector<libint_basis_type>;
    using const_lbv_reference     = const libint_basis_vector&;
    using operator_type           = OpType;
    using threshold_type          = double;
    using derivative_order_type   = int;
    using const_buffer_reference  = const buffer_type&;
    using const_indices_reference = const indices_type&;

    LibintFactory(libint_basis_vector bases, operator_type op,
                  threshold_type thesh, derivative_order_type deriv);

private:
    const_buffer_reference compute_(const_indices_reference idx) override;
    pimpl_ptr_type clone_() const override;
    bool are_equal_(const base_type& rhs) const override;

    /** @brief Wrap the call of LibInt2 engine so it can take a variable
     * number of shell inputs.
     *
     *  Libint requires the engine for an `n`-center integral to be called
     * with `n` positional arguments. The API of `IntegralFactory` relies on
     *  `std::vector` to avoid needing to know the number of centers at
     * compile time. This function wraps the process of unpacking the
     * `std::vector` of indices.
     *
     * @tparam Is A variadic parameter pack of integers from [0,NCenters) to
     *         expand.
     *
     * @param[in] shells The indices of the requested shell block.
     * @param[in] <anonymous> An `index_sequence` used to get the compile time
     *                        indices for expanding @p shells.
     */
    template<std::size_t... Is>
    void run_engine_(const_indices_reference shells,
                     std::index_sequence<Is...>);

    libint_basis_vector m_bases_;
    operator_type m_op_;
    threshold_type m_thresh_;
    derivative_order_type m_deriv_;
    libint2::Engine m_engine_;
    buffer_type m_buffer_;
};

#define LIBINT_FACTORY LibintFactory<NCenters, OpType>

template<std::size_t NCenters, typename OpType>
LIBINT_FACTORY::LibintFactory(libint_basis_vector bases, OpType op,
                              threshold_type thresh,
                              derivative_order_type deriv) :
  m_bases_(std::move(bases)),
  m_op_(op),
  m_thresh_(thresh),
  m_deriv_(deriv),
  m_engine_(detail_::make_engine(m_bases_, m_op_, m_thresh_, m_deriv_)) {}

template<std::size_t NCenters, typename OpType>
typename LIBINT_FACTORY::const_buffer_reference LIBINT_FACTORY::compute_(
  const_indices_reference idx) {
    auto idx_seq = std::make_index_sequence<NCenters>();
    run_engine_(idx, idx_seq);
    const auto& results = m_engine_.results();
    m_buffer_           = buffer_type(results.begin(), results.end());
    return m_buffer_;
}

template<std::size_t NCenters, typename OpType>
typename LIBINT_FACTORY::pimpl_ptr_type LIBINT_FACTORY::clone_() const {
    return std::make_unique<my_type>(m_bases_, m_op_, m_thresh_, m_deriv_);
}

template<std::size_t NCenters, typename OpType>
bool LIBINT_FACTORY::are_equal_(const base_type& rhs) const {
    auto ptr = dynamic_cast<const my_type*>(&rhs);
    if(ptr == nullptr) return false;
    return (m_bases_ == ptr->m_bases_) && (m_op_ == ptr->m_op_) &&
           (m_thresh_ == ptr->m_thresh_) && (m_deriv_ == ptr->m_deriv_) &&
           (m_buffer_ == ptr->m_buffer_);
}

template<std::size_t NCenters, typename OpType>
template<std::size_t... Is>
void LIBINT_FACTORY::run_engine_(const_indices_reference shells,
                                 std::index_sequence<Is...>) {
    m_engine_.compute(m_bases_[Is][shells[Is]]...);
}

#undef LIBINT_FACTORY

} // namespace integrals::libint
