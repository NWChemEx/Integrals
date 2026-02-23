#include "utils.hpp"
#include <integrals/integrals.hpp>
namespace integrals::utils {

namespace {

const auto desc = R"(
Screen Primitive Pairs
======================

This module returns a set of primitive pairs whose estimated contribution to
the overall ERI is greater than a user specified threshold.
)";

}

using pair_estimate_pt = integrals::property_types::PrimitivePairEstimator;
using pt               = integrals::property_types::PairScreener;

MODULE_CTOR(ScreenPrimitivePairs) {
    description(desc);
    satisfies_property_type<pt>();

    add_submodule<pair_estimate_pt>("Primitive Pair Estimator")
      .set_description("The module used to estimate the contributions of "
                       "primitive pairs to the overall integral values");
}

MODULE_RUN(ScreenPrimitivePairs) {
    const auto&& [bra, ket, tol] = pt::unwrap_inputs(inputs);

    auto& screener  = submods.at("Primitive Pair Estimator");
    const auto& Kij = screener.run_as<pair_estimate_pt>(bra, ket);

    const auto n_bra_prims = bra.n_primitives();
    const auto n_ket_prims = ket.n_primitives();

    using float_type = double;
    using tensorwrapper::buffer::make_contiguous;
    const auto& buffer = make_contiguous(Kij.buffer());

    using index_vector = std::vector<std::size_t>;
    index_vector prim{0, 0};
    std::vector<index_vector> rv;
    using wtf::fp::float_cast;

    for(prim[0] = 0; prim[0] < n_bra_prims; ++prim[0]) {
        for(prim[1] = 0; prim[1] < n_ket_prims; ++prim[1]) {
            const auto val = float_cast<float_type>(buffer.get_elem(prim));
            if(val > tol) rv.push_back(prim);
        }
    }

    auto result = results();
    return pt::wrap_results(result, rv);
}

} // namespace integrals::utils