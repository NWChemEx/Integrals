#include "integrals/nwx_libint/nwx_libint_factory.hpp"

namespace nwx_libint {

    template<size_type NBases, libint2::Operator op>
    libint2::Engine LibintFactory<NBases, op>::operator()() {
        libint2::Engine engine(op, max_nprims, max_l, deriv, thresh);

        // Take care of any special set-up
        if constexpr(libint2::rank(op) == 2 && NBases == 2) {
            engine.set(libint2::BraKet::xs_xs);
        }

        if constexpr(libint2::rank(op) == 2 && NBases == 3) {
            engine.set(libint2::BraKet::xs_xx);
        }

        if constexpr(op == libint2::Operator::nuclear) {
            std::vector<std::pair<double, std::array<double, 3>>> qs;

            for(const auto& ai : mol)
                qs.push_back({static_cast<const double&>(ai.Z()), ai.coords()});

            engine.set_params(qs);

        } else if constexpr (op == libint2::Operator::stg ||
                             op == libint2::Operator::yukawa) {
            engine.set_params(stg_exponent);
        } else if constexpr (op == libint2::Operator::emultipole1 ||
                             op == libint2::Operator::emultipole2 ||
                             op == libint2::Operator::emultipole3 ) {
            engine.set_params(origin);
        }

        return engine;
    }

    template class LibintFactory<2, libint2::Operator::overlap>;
    template class LibintFactory<2, libint2::Operator::kinetic>;
    template class LibintFactory<2, libint2::Operator::nuclear>;
    template class LibintFactory<2, libint2::Operator::coulomb>;
    template class LibintFactory<3, libint2::Operator::coulomb>;
    template class LibintFactory<4, libint2::Operator::coulomb>;
    template class LibintFactory<2, libint2::Operator::stg>;
    template class LibintFactory<3, libint2::Operator::stg>;
    template class LibintFactory<4, libint2::Operator::stg>;
    template class LibintFactory<2, libint2::Operator::yukawa>;
    template class LibintFactory<3, libint2::Operator::yukawa>;
    template class LibintFactory<4, libint2::Operator::yukawa>;
    template class LibintFactory<2, libint2::Operator::emultipole1>;
    template class LibintFactory<2, libint2::Operator::emultipole2>;
    template class LibintFactory<2, libint2::Operator::emultipole3>;
    template class LibintFactory<4, libint2::Operator::delta>;
    
} // namespace nwx_libint