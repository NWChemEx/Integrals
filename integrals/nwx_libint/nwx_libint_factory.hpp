#pragma once
#include <libint2.hpp>
#include "integrals/types.hpp"

namespace nwx_libint {

    using size_type = integrals::type::size;
    using mol_type = integrals::type::molecule;

    // Factory class that produces Libint2 engines.
    template<size_type NBases, libint2::Operator op>
    struct LibintFactory {

        // general parameters of the Libint engine
        size_type max_nprims = 0;
        size_type max_l = 0;
        double thresh = 0;
        size_type deriv = 0;

        // parameters for specific integral types
        double stg_exponent = 1.0;
        std::array<double, 3> origin{0, 0, 0};
        mol_type mol;

        LibintFactory() = default;

        LibintFactory(size_type max_nprims, size_type max_l, double thresh, size_type deriv) :
                      max_nprims(max_nprims), max_l(max_l), thresh(thresh), deriv(deriv) {}

        // produce a LibInt2 engine, given the current parameters
        libint2::Engine operator()();

    }; // Class LibIntFactory

    extern template class LibintFactory<2, libint2::Operator::overlap>;
    extern template class LibintFactory<2, libint2::Operator::kinetic>;
    extern template class LibintFactory<2, libint2::Operator::nuclear>;
    extern template class LibintFactory<2, libint2::Operator::coulomb>;
    extern template class LibintFactory<3, libint2::Operator::coulomb>;
    extern template class LibintFactory<4, libint2::Operator::coulomb>;
    extern template class LibintFactory<2, libint2::Operator::stg>;
    extern template class LibintFactory<3, libint2::Operator::stg>;
    extern template class LibintFactory<4, libint2::Operator::stg>;
    extern template class LibintFactory<2, libint2::Operator::yukawa>;
    extern template class LibintFactory<3, libint2::Operator::yukawa>;
    extern template class LibintFactory<4, libint2::Operator::yukawa>;
    extern template class LibintFactory<2, libint2::Operator::emultipole1>;
    extern template class LibintFactory<2, libint2::Operator::emultipole2>;
    extern template class LibintFactory<2, libint2::Operator::emultipole3>;
    extern template class LibintFactory<4, libint2::Operator::delta>;


} // namespace nwx_libint
