#include "detail_/bases_helper.hpp"
#include "detail_/make_engine.hpp"
#include "shellnorms.hpp"
#include <simde/cauchy_schwarz_approximation.hpp>

namespace integrals {

template<typename OperatorType>
TEMPLATED_MODULE_CTOR(ShellNorms, OperatorType) {
    description("Calculates the Cauchy-Schwarz screening matrix for a pair of "
                "basis sets");

    using my_pt = simde::ShellNorms<OperatorType>;
    satisfies_property_type<my_pt>();

    add_input<double>("Threshold")
      .set_description("Convergence threshold of integrals")
      .set_default(1.0E-16);
}

template<typename OperatorType>
TEMPLATED_MODULE_RUN(ShellNorms, OperatorType) {
    using my_pt      = simde::ShellNorms<OperatorType>;
    using elem_vec   = typename std::vector<double>;
    using return_vec = typename std::vector<elem_vec>;

    // Get inputs
    auto bases  = detail_::unpack_bases<2>(inputs);
    auto op_str = OperatorType().as_string();
    auto op     = inputs.at(op_str).template value<const OperatorType&>();
    auto thresh = inputs.at("Threshold").value<double>();

    // Check if the basis sets are the same
    bool same_bs = (bases[0] == bases[1]);

    // Lambda to fill in the values
    return_vec mat(bases[0].size(), elem_vec(bases[1].size(), 0.0));
    auto into_mat = [&](int i, int j) {
        auto engine     = detail_::make_engine(bases, op, thresh);
        const auto& buf = engine.results();

        engine.compute(bases[0][i], bases[1][j], bases[0][i], bases[1][j]);
        auto vals = buf[0];

        // Determine the number of compute values
        std::size_t nvals = (bases[0][i].size() * bases[0][i].size());
        nvals *= (bases[1][j].size() * bases[1][j].size());

        // Find the norm and take the square root
        double infinity_norm = 0.0;
        if(vals != nullptr) {
            for(int a = 0; a < nvals; ++a) {
                infinity_norm = std::max(infinity_norm, std::abs(vals[a]));
            }
        }
        mat[i][j] = std::sqrt(infinity_norm);
        if(same_bs && (i != j)) { mat[j][i] = mat[i][j]; } // cut down on work
    };

    // Calculate values
    auto& my_runtime = get_runtime();
    auto& world      = my_runtime.madness_world();
    for(int i = 0; i < bases[0].size(); ++i) {
        // only do lower triangle if basis sets are the same
        auto len = (same_bs) ? i : bases[1].size() - 1;
        for(int j = 0; j <= len; ++j) { world.taskq.add(into_mat, i, j); }
    }
    world.gop.fence();

    auto rv = results();
    return my_pt::wrap_results(rv, mat);
}

template class ShellNorms<simde::type::el_el_coulomb>;
template class ShellNorms<simde::type::el_el_stg>;
template class ShellNorms<simde::type::el_el_yukawa>;

} // namespace integrals
