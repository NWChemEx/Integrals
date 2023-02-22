#include <chemist/integral_factory/pimpl.hpp>
#include <libint2.hpp>

namespace integrals::libint {

template<std::size_t NCenters>
class LibintFactory : public chemist::IntegralFactoryPIMPL {
public:
    using libint_basis_vector = std::vector<libint2::BasisSet>;
    using const_lbv_reference = const libint_basis_vector&;

    using const_buffer_reference  = const buffer_t&;
    using const_indices_reference = const indices_t&;

private:
    const_buffer_reference compute_(const_indices_reference idx) const override;

    std::vector<libint2::BasisSet> m_bases_;
    libint2::Engine m_engine_;
};

#define LIBINT_FACTORY LibintFactory<NCenters>

template<std::size_t NCenters>
LIBINT_FACTORY::LibintFactory(libint_basis_vector bases, op, thresh, deriv) :
  m_bases_(std::move(bases)),
  m_engine_(make_engine(m_bases_, op, thresh, deriv)) {}

template<std::size_t NCenters>
typename LIBINT_FACTORY::const_buffer_reference LIBINT_FACTORY::compute_(
  const_indices_reference idx) const {
    auto idx_seq = std::make_index_sequence<NCenters>();
    run_engine_(m_engine_, m_bases_, idx, idx_seq);
    return engine.results();
}

#undef LIBINT_FACTORY

} // namespace integrals::libint