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
#include "../uncertain_types.hpp"
#include "detail_/get_basis_sets.hpp"
#include "detail_/libint_op.hpp"
#include "detail_/make_engine.hpp"
#include "detail_/make_libint_basis_set.hpp"
#include "detail_/shells2ord.hpp"
#include "libint.hpp"
#include "libint_visitor.hpp"
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace integrals::libint {
namespace {

#ifdef _OPENMP
int get_num_threads() {
    int num_threads;
#pragma omp parallel
    { num_threads = omp_get_num_threads(); }
    return num_threads;
}

int get_thread_num() { return omp_get_thread_num(); }
#else

int get_num_threads() { return 1; }

int get_thread_num() { return 0; }

#endif

template<typename FloatType>
auto build_eigen_buffer(const std::vector<libint2::BasisSet>& basis_sets,
                        parallelzone::runtime::RuntimeView& rv, double thresh) {
    FloatType initial_value;
    if constexpr(std::is_same_v<FloatType, double>) {
        initial_value = 0.0;
    } else { // Presumably sigma::UDouble
        initial_value = FloatType(0.0, thresh);
    }
    auto N = basis_sets.size();
    std::vector<decltype(N)> dims(N);
    for(decltype(N) i = 0; i < N; ++i) dims[i] = basis_sets[i].nbf();

    using namespace tensorwrapper;
    using shape_t = shape::Smooth;

    shape_t s{dims.begin(), dims.end()};
    auto buffer = buffer::make_contiguous<FloatType>(s, initial_value);
    return std::make_unique<buffer::Contiguous>(std::move(buffer));
}

template<std::size_t N, typename FloatType>
auto fill_tensor(const std::vector<libint2::BasisSet>& basis_sets,
                 const chemist::qm_operator::OperatorBase& op,
                 parallelzone::runtime::RuntimeView& rv, double thresh) {
    using size_type = decltype(N);

    // Dimensional information
    std::vector<size_type> dim_stepsizes(N, 1);
    size_type num_shell_combinations = 1;

    for(size_type i = 0; i < N; ++i) {
        num_shell_combinations *= basis_sets[i].size();
        for(size_type j = i; j < N - 1; ++j) {
            dim_stepsizes[i] *= basis_sets[j].size();
        }
    }

    // Make an engine for each thread
    int num_threads = get_num_threads();
    std::vector<libint2::Engine> engines(num_threads);
    LibintVisitor visitor(basis_sets, thresh);
    op.visit(visitor);
    for(int i = 0; i != num_threads; ++i) { engines[i] = visitor.engine(); }

    // Fill in values
    auto pbuffer = build_eigen_buffer<FloatType>(basis_sets, rv, thresh);
    auto data    = pbuffer->get_mutable_data();
    auto span    = wtf::buffer::contiguous_buffer_cast<FloatType>(data);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(size_type i_pair = 0; i_pair != num_shell_combinations; ++i_pair) {
        auto thread_id = get_thread_num();

        std::vector<size_type> shells(N);
        auto shell_ord = i_pair;
        for(size_type i = 0; i < N; ++i) {
            shells[i] = shell_ord / dim_stepsizes[i];
            shell_ord = shell_ord % dim_stepsizes[i];
        }

        detail_::run_engine_(engines[thread_id], basis_sets, shells,
                             std::make_index_sequence<N>());

        const auto& buf = engines[thread_id].results();
        auto vals       = buf[0];
        if(vals) {
            auto ord   = detail_::shells2ord(basis_sets, shells);
            auto n_ord = ord.size();
            for(decltype(n_ord) i_ord = 0; i_ord < n_ord; ++i_ord) {
                auto update      = span[ord[i_ord]] + vals[i_ord];
                span[ord[i_ord]] = update;
            }
        }
    }

    auto pshape = pbuffer->layout().shape().clone();
    return simde::type::tensor(std::move(pshape), std::move(pbuffer));
}

} // namespace

template<typename BraKetType>
TEMPLATED_MODULE_CTOR(Libint, BraKetType) {
    using my_pt = simde::EvaluateBraKet<BraKetType>;
    satisfies_property_type<my_pt>();
    description("Driver for computing integrals with Libint");

    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");

    add_input<bool>("With UQ?").set_default(false);
}

template<typename BraKetType>
TEMPLATED_MODULE_RUN(Libint, BraKetType) {
    using my_pt = simde::EvaluateBraKet<BraKetType>;

    const auto& [braket] = my_pt::unwrap_inputs(inputs);
    auto thresh          = inputs.at("Threshold").value<double>();
    auto with_uq         = inputs.at("With UQ?").value<bool>();
    auto bra             = braket.bra();
    auto ket             = braket.ket();
    auto& op             = braket.op();
    auto& rv             = this->get_runtime();

    // Gather information from Bra, Ket, and Op
    auto basis_sets = detail_::get_basis_sets(bra, ket);
    constexpr int N = detail_::get_n(bra, ket);

    simde::type::tensor t;
    if(with_uq) {
        if constexpr(integrals::type::has_sigma()) {
            t = fill_tensor<N, type::uncertain_double>(basis_sets, op, rv,
                                                       thresh);
        } else {
            throw std::runtime_error("Sigma support not enabled!");
        }
    } else {
        t = fill_tensor<N, double>(basis_sets, op, rv, thresh);
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
    mm.change_submod("Primitive Error Model", "Primitive Pair Estimator",
                     "Black Box Primitive Pair Estimator");
    mm.change_submod("CauchySchwarz Estimator", "Decontract Basis Set",
                     "Decontract Basis Set");
    mm.change_submod("CauchySchwarz Estimator", "ERI4", "ERI4");
    mm.change_submod("Analytic Error", "ERI4s", "ERI4");
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
}

#undef LOAD_LIBINT
#undef LIBINT

} // namespace integrals::libint
