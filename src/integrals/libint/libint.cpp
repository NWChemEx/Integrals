/*
 * Copyright 2024 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "detail_/fill_tensor.hpp"
#include "detail_/get_basis_sets.hpp"
#include "libint.hpp"
namespace integrals::libint {

template<typename BraKetType>
TEMPLATED_MODULE_CTOR(Libint, BraKetType) {
    using my_pt = simde::EvaluateBraKet<BraKetType>;
    satisfies_property_type<my_pt>();
    description("Driver for computing integrals with Libint");

    add_input<std::string>("UQ Type").set_default("none");
    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");
}

template<typename BraKetType>
TEMPLATED_MODULE_RUN(Libint, BraKetType) {
    using my_pt = simde::EvaluateBraKet<BraKetType>;

    const auto& [braket] = my_pt::unwrap_inputs(inputs);
    auto thresh          = inputs.at("Threshold").value<double>();
    auto uq_type         = inputs.at("UQ Type").value<std::string>();
    auto bra             = braket.bra();
    auto ket             = braket.ket();
    auto& op             = braket.op();
    auto& rv             = this->get_runtime();

    // Gather information from Bra, Ket, and Op
    auto basis_sets = detail_::get_basis_sets(bra, ket);
    constexpr int N = detail_::get_n(bra, ket);
    simde::type::tensor t;
    using float_type = double;
    if(uq_type == "none") {
        t = detail_::fill_tensor<N, float_type>(basis_sets, op, rv, thresh);
    } else if(uq_type == "uncertain") {
        using uncertain_type = tensorwrapper::types::udouble;
        t = detail_::fill_tensor<N, uncertain_type>(basis_sets, op, rv, thresh);
    } else if(uq_type == "interval") {
        using interval_type = tensorwrapper::types::interval_type<float_type>;
        t = detail_::fill_tensor<N, interval_type>(basis_sets, op, rv, thresh);
    } else {
        throw std::runtime_error("Invalid UQ type");
    }
    auto result = results();
    return my_pt::wrap_results(result, t);
}

#define LIBINT(bra, op, ket) Libint<braket<bra, op, ket>>
#define EXTERN_LIBINT(bra, op, ket) template struct LIBINT(bra, op, ket)

EXTERN_LIBINT(aos, op_base_type, aos);
EXTERN_LIBINT(aos, op_base_type, aos_squared);
EXTERN_LIBINT(aos_squared, op_base_type, aos_squared);
EXTERN_LIBINT(aos, s_e_type, aos);
EXTERN_LIBINT(aos, t_e_type, aos);
EXTERN_LIBINT(aos, v_en_type, aos);
EXTERN_LIBINT(aos, v_ee_type, aos);
EXTERN_LIBINT(aos, v_ee_type, aos_squared);
EXTERN_LIBINT(aos_squared, v_ee_type, aos_squared);

#undef EXTERN_LIBINT

void set_defaults(pluginplay::ModuleManager& mm) {
    mm.change_submod("CauchySchwarz Estimator", "Decontract Basis Set",
                     "Decontract Basis Set");
    mm.copy_module("ERI4", "Benchmark ERI4");
    mm.change_input("Benchmark ERI4", "Threshold", 1.0E-16);
    mm.change_submod("CauchySchwarz Estimator", "ERI4", "Benchmark ERI4");
    mm.change_submod("Analytic Error", "ERI4s", "Benchmark ERI4");
    mm.change_submod("Raw Primitive ERI4", "Decontract Basis Set",
                     "Decontract Basis Set");
    mm.change_submod("Primitive Contractor ERI4", "Raw Primitive ERI4",
                     "Raw Primitive ERI4");
    mm.change_submod("Primitive Contractor ERI4", "Primitive Normalization",
                     "Primitive Normalization");
}

#define LOAD_LIBINT(bra, op, ket, key) mm.add_module<LIBINT(bra, op, ket)>(key)

void load_modules(pluginplay::ModuleManager& mm) {
    LOAD_LIBINT(aos, op_base_type, aos, "Evaluate 2-Index BraKet");
    LOAD_LIBINT(aos, op_base_type, aos_squared, "Evaluate 3-Index BraKet");
    LOAD_LIBINT(aos_squared, op_base_type, aos_squared,
                "Evaluate 4-Index BraKet");
    LOAD_LIBINT(aos, s_e_type, aos, "Overlap");
    LOAD_LIBINT(aos, t_e_type, aos, "Kinetic");
    LOAD_LIBINT(aos, v_en_type, aos, "Nuclear");
    LOAD_LIBINT(aos, v_ee_type, aos, "ERI2");
    LOAD_LIBINT(aos, v_ee_type, aos_squared, "ERI3");
    LOAD_LIBINT(aos_squared, v_ee_type, aos_squared, "ERI4");
    mm.add_module<BlackBoxPrimitiveEstimator>(
      "Black Box Primitive Pair Estimator");
    mm.add_module<CauchySchwarzPrimitiveEstimator>("CauchySchwarz Estimator");
    mm.add_module<PrimitiveErrorModel>("Primitive Error Model");
    mm.add_module<AnalyticError>("Analytic Error");
    mm.add_module<RawPrimitiveERIs>("Raw Primitive ERI4");
    mm.add_module<PrimitiveNormalization>("Primitive Normalization");
    mm.add_module<PrimitiveContractor>("Primitive Contractor ERI4");
}

#undef LOAD_LIBINT
#undef LIBINT

} // namespace integrals::libint
