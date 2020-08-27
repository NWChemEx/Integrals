#include "integrals/cauchy_schwarz.hpp"
#include "integrals/nwx_libint/nwx_libint_factory.hpp"
#include "integrals/nwx_libint/nwx_libint.hpp"
#include <property_types/cauchy_schwarz_approximation.hpp>
#include <tiledarray.h>

namespace integrals {

    template<typename element_type>
    using cs_approx_type = property_types::CauchySchwarzApprox<element_type>;

    template<typename element_type, libint2::Operator op>
    CauchySchwarz<element_type, op>::CauchySchwarz() : sde::ModuleBase(this) {
        description("Calculates the Cauchy-Schwarz screening matrix for a pair of basis sets");
        satisfies_property_type<cs_approx_type<element_type>>();

        add_input<element_type>("Threshold")
                .set_description("Convergence threshold of integrals")
                .set_default(1.0E-16);
    }

    template<typename element_type, libint2::Operator op>
    sde::type::result_map CauchySchwarz<element_type, op>::run_(sde::type::input_map inputs,
                                                                sde::type::submodule_map submods) const {
        using elem_vec = typename std::vector<element_type>;
        using return_vec = typename std::vector<elem_vec>;
        using basis_type = libint2::BasisSet;
        using basis_vec = std::vector<basis_type>;

        // Get inputs
        auto [basis1, basis2, deriv] = cs_approx_type<element_type>::unwrap_inputs(inputs);
        auto thresh = inputs.at("Threshold").value<element_type>();

        // Set up
        auto bs1 = nwx_libint::make_basis(basis1);
        auto bs2 = nwx_libint::make_basis(basis2);
        auto factory = nwx_libint::LibintFactory();
        factory.max_nprims = nwx_libint::sets_max_nprims({bs1, bs2});
        factory.max_l = nwx_libint::sets_max_l({bs1, bs2});
        factory.thresh = thresh;
        factory.deriv = deriv;

        // In case it was finalized
        if (not libint2::initialized()) { libint2::initialize(); }

        // Place for values to go
        return_vec mat(bs1.size(), elem_vec(bs2.size(), 0.0));

        // Check if the basis sets are the same
        bool same_bs = (bs1 == bs2);

        // Lambda to fill in the values
        auto into_mat = [&] (int i, int j) {
            auto engine = factory(4, op);
            const auto& buf = engine.results();

            engine.compute(bs1[i], bs2[j], bs1[i], bs2[j]);

            auto vals = buf[0];

            // Determine the number of compute values
            std::size_t nvals = (bs1[i].size() * bs1[i].size());
            nvals *= (bs2[j].size() * bs2[j].size());

            // Find the norm and take the square root
            double infinity_norm = 0.0;
            if (vals != nullptr) {
                for (int a = 0; a < nvals; ++a) { infinity_norm = std::max(infinity_norm, std::abs(vals[a])); }
            }
            mat[i][j] = std::sqrt(infinity_norm);
            if (same_bs && (i!=j)) { mat[j][i] = mat[i][j]; } // cut down on work
        };

        // Calculate values
        auto& world = TA::get_default_world();
        for (int i = 0; i < bs1.size(); ++i) {
            auto len = (same_bs) ? i : bs2.size() - 1; // only do lower triangle, since it's mirrored
            for (int j = 0; j <= len; ++j) {
                world.taskq.add(into_mat, i, j);
            }
        }
        world.gop.fence();

        auto rv = results();
        return cs_approx_type<element_type>::wrap_results(rv, mat);
    }

    template class CauchySchwarz<double>;
    template class CauchySchwarz<double, libint2::Operator::stg>;
    template class CauchySchwarz<double, libint2::Operator::yukawa>;

} // namespace integrals
