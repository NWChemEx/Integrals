#pragma once
#include "integrals/types.hpp"
#include <libint2.hpp>

namespace nwx_libint {

// Factory class that produces Libint2 engines.
struct LibintFactory {
    using size_type = integrals::type::size;

    // general parameters of the Libint engine
    size_type max_nprims = 0;
    size_type max_l      = 0;
    double thresh        = 0;
    size_type deriv      = 0;

    // parameters for specific integral types
    double stg_exponent = 1.0;
    std::array<double, 3> origin{0, 0, 0};
    std::vector<std::pair<double, std::array<double, 3>>> qs;

    // produce a LibInt2 engine, given the current parameters
    libint2::Engine operator()(size_type NBases, libint2::Operator op);

}; // Class LibIntFactory

} // namespace nwx_libint
