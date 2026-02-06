#include "libint.hpp"
#include <integrals/integrals.hpp>

namespace integrals::libint {
namespace {

const auto desc = "Uses the error in the ERI4 as the uncertainty.";
}

using eri4_pt = simde::ERI4;
using pt      = integrals::property_types::Uncertainty<eri4_pt>;

MODULE_CTOR(AnalyticError) {
    satisfies_property_type<pt>();
    description(desc);

    add_submodule<eri4_pt>("ERI4s");
}

MODULE_RUN(AnalyticError) {
    const auto& [braket, tol] = pt::unwrap_inputs(inputs);

    auto& eri_mod = submods.at("ERI4s");

    auto normal_mod = eri_mod.value().unlocked_copy();
    normal_mod.change_input("Threshold", tol);

    // N.b., t_0 is the benchmark value
    const auto& t_0 = eri_mod.run_as<eri4_pt>(braket);
    const auto& t   = normal_mod.run_as<eri4_pt>(braket);

    simde::type::tensor error;
    error("m,n,l,s") = t("m,n,l,s") - t_0("m,n,l,s");

    // Wrap and return the results
    auto rv = results();
    return pt::wrap_results(rv, error);
}

} // namespace integrals::libint