#include <chemist/integral_factory/pimpl.hpp>
#include <libint2.hpp>

namespace integrals::libint {

template<typename Operator>
class LibintFactory : public chemist::IntegralFactoryPIMPL {
private:
    libint2::Engine m_engine_;
};

template<typename Operator>
LibintFactory::LibintFactory(linbint2::Engine engine) :
  m_engine_(std::move(engine)) {}

} // namespace integrals::libint