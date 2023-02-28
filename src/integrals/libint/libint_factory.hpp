#include <chemist/integral_factory/integral_factory_pimpl.hpp>
#include <libint2.hpp>

namespace integrals::libint {

template<std::size_t NCenters>
class LibintFactory : public chemist::IntegralFactoryPIMPL {
public:
    using libint_basis_type     = libint2::BasisSet;
    using libint_basis_vector   = std::vector<libint_basis_type>;
    using const_lbv_reference   = const libint_basis_vector&;
    using libint_op_type        = libint2::Operator;
    using threshold_type        = double;
    using derivative_order_type = int;

    using const_buffer_reference  = const buffer_t&;
    using const_indices_reference = const indices_t&;

    LibintFactory(libint_basis_vector bases, libint_op_type op,
                  threshold_type thesh, derivative_order_type deriv = 0);

private:
    const_buffer_reference compute_(const_indices_reference idx) const override;

    /** @brief Wrap the call of LibInt2 engine so it can take a variable number
     *         of shell inputs.
     *
     *  Libint requires the engine for an `n`-center integral to be called with
     *  `n` positional arguments. The API of `IntegralFactory` relies on
     *  `std::vector` to avoid needing to know the number of centers at compile
     *  time. This function wraps the process of unpacking the `std::vector` of
     *  indices.
     *
     * @tparam Is A variadic parameter pack of integers from [0,NCenters) to
     *         expand.
     *
     * @param[in] shells The indices of the requested shell block.
     * @param[in] <anonymous> An `index_sequence` used to get the compile time
     *                        indices for expanding @p shells.
     */
    template<std::size_t... Is>
    void run_engine_(const_indices_reference& shells,
                     std::index_sequence<Is...>);

    std::vector<libint2::BasisSet> m_bases_;
    libint_op_type m_op_;
    threshold_type m_thresh_;
    derivative_order_type m_deriv_;
    libint2::Engine m_engine_;
};

#define LIBINT_FACTORY LibintFactory<NCenters>

template<std::size_t NCenters>
LIBINT_FACTORY::LibintFactory(libint_basis_vector bases, libint_op_type op,
                              threshold_type thresh,
                              derivative_order_type deriv) :
  m_bases_(std::move(bases)),
  m_op_(op),
  m_thresh_(thresh),
  m_deriv_(deriv),
  m_engine_(make_engine(m_bases_, m_op_, m_thresh_, m_deriv_)) {}

template<std::size_t NCenters>
typename LIBINT_FACTORY::const_buffer_reference LIBINT_FACTORY::compute_(
  const_indices_reference idx) const {
    auto idx_seq = std::make_index_sequence<NCenters>();
    run_engine_(m_engine_, m_bases_, idx, idx_seq);
    return m_engine_.results();
}

template<std::size_t NCenters>
template<std::size_t... Is>
void LIBINT_FACTORY::run_engine_(const_indices_reference shells,
                                 std::index_sequence<Is...>) {
    m_engine_.compute(m_bases_[Is][shells[Is]]...);
}

#undef LIBINT_FACTORY

} // namespace integrals::libint
