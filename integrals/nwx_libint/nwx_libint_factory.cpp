#include "integrals/nwx_libint/nwx_libint_factory.hpp"

namespace nwx_libint {

libint2::Engine LibintFactory::operator()(size_type NBases,
                                          libint2::Operator op) {
    libint2::Engine engine(op, max_nprims, max_l, deriv, thresh);

    // Take care of any special set-up
    if(libint2::rank(op) == 2 && NBases == 2) {
        engine.set(libint2::BraKet::xs_xs);
    }

    if(libint2::rank(op) == 2 && NBases == 3) {
        engine.set(libint2::BraKet::xs_xx);
    }

    if(op == libint2::Operator::nuclear) {
        engine.set_params(qs);

    } else if(op == libint2::Operator::stg || op == libint2::Operator::yukawa) {
        engine.set_params(stg_exponent);
    } else if(op == libint2::Operator::emultipole1 ||
              op == libint2::Operator::emultipole2 ||
              op == libint2::Operator::emultipole3) {
        engine.set_params(origin);
    }

    return engine;
}

} // namespace nwx_libint