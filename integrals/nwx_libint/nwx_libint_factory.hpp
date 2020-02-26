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
        libint2::Engine operator()() {
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

    }; // Class LibIntFactory

} // namespace nwx_libint
