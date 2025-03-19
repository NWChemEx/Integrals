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
#include "libint_uq.hpp"
#include "libint_visitor.hpp"
#include <type_traits>

namespace integrals::libint::uq {
namespace {

auto build_eigen_buffer(const std::vector<libint2::BasisSet>& basis_sets,
                        parallelzone::runtime::RuntimeView& rv) {
    type::uncertain_double initial_value(0.0, 0.0);
    auto N = basis_sets.size();
    std::vector<decltype(N)> dims(N);
    for(decltype(N) i = 0; i < N; ++i) dims[i] = basis_sets[i].nbf();

    using shape_t  = tensorwrapper::shape::Smooth;
    using layout_t = tensorwrapper::layout::Physical;

    shape_t s{dims.begin(), dims.end()};
    layout_t l(s);
    tensorwrapper::allocator::Eigen<type::uncertain_double> alloc(rv);
    return alloc.construct(l, initial_value);
}

template<std::size_t N>
auto fill_tensor(const std::vector<libint2::BasisSet>& basis_sets,
                 const chemist::qm_operator::OperatorBase& op,
                 parallelzone::runtime::RuntimeView& rv, double thresh,
                 double bench_thresh) {
    // Dimensional information
    std::vector<std::size_t> dims_shells(N);
    for(decltype(N) i = 0; i < N; ++i) dims_shells[i] = basis_sets[i].size();

    auto pbuffer = build_eigen_buffer(basis_sets, rv);

    // Make libint engine
    LibintVisitor visitor(basis_sets, thresh);
    op.visit(visitor);
    auto engine     = visitor.engine();
    const auto& buf = engine.results();

    // Make benchmark engine
    LibintVisitor bench_visitor(basis_sets, bench_thresh);
    op.visit(bench_visitor);
    auto bench_engine     = bench_visitor.engine();
    const auto& bench_buf = bench_engine.results();

    // Fill in values
    std::vector<std::size_t> shells(N, 0);
    while(shells[0] < dims_shells[0]) {
        detail_::run_engine_(engine, basis_sets, shells,
                             std::make_index_sequence<N>());
        detail_::run_engine_(bench_engine, basis_sets, shells,
                             std::make_index_sequence<N>());

        auto vals       = buf[0];
        auto benchmarks = bench_buf[0];
        if(vals) {
            auto ord   = detail_::shells2ord(basis_sets, shells);
            auto n_ord = ord.size();
            for(decltype(n_ord) i_ord = 0; i_ord < n_ord; ++i_ord) {
                double value = vals[i_ord];
                // Assuming benchmarks return as zeroes, the error defaults to
                // the current value. Otherwise, the error is the difference
                // between the current value and the benchmark value.
                double error = vals[i_ord];
                if(benchmarks) { error = value - benchmarks[i_ord]; }
                pbuffer->data()[ord[i_ord]] =
                  type::uncertain_double(value, error);
            }
        }

        // Increment index
        shells[N - 1] += 1;
        for(decltype(N) i = 1; i < N; ++i) {
            if(shells[N - i] >= dims_shells[N - i]) {
                // Reset this dimension and increment the next one
                // shells[0] accumulates until we reach the end
                shells[N - i] = 0;
                shells[N - i - 1] += 1;
            }
        }
    }

    auto pshape = pbuffer->layout().shape().clone();
    return simde::type::tensor(std::move(pshape), std::move(pbuffer));
}

} // namespace

template<typename BraKetType>
TEMPLATED_MODULE_CTOR(LibintUQ, BraKetType) {
    using my_pt = simde::EvaluateBraKet<BraKetType>;
    satisfies_property_type<my_pt>();
    description("Driver for computing integrals with Libint");

    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");

    add_input<double>("Benchmark Threshold")
      .set_default(1.0e-16)
      .set_description(
        "The target precision with which the uncertainty will be computed");
}

template<typename BraKetType>
TEMPLATED_MODULE_RUN(LibintUQ, BraKetType) {
    using my_pt = simde::EvaluateBraKet<BraKetType>;

    const auto& [braket] = my_pt::unwrap_inputs(inputs);
    auto thresh          = inputs.at("Threshold").value<double>();
    auto bench_thresh    = inputs.at("Benchmark Threshold").value<double>();
    auto bra             = braket.bra();
    auto ket             = braket.ket();
    auto& op             = braket.op();
    auto& rv             = this->get_runtime();

    // Gather information from Bra, Ket, and Op
    auto basis_sets = detail_::get_basis_sets(bra, ket);
    constexpr int N = detail_::get_n(bra, ket);

    simde::type::tensor t;
    if constexpr(integrals::type::has_sigma()) {
        t = fill_tensor<N>(basis_sets, op, rv, thresh, bench_thresh);
    } else {
        throw std::runtime_error("Sigma support not enabled!");
    }

    auto result = results();
    return my_pt::wrap_results(result, t);
}

#define LIBINT(bra, op, ket) LibintUQ<braket<bra, op, ket>>
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
    // Set any default associations
}

#define LOAD_LIBINT(bra, op, ket, key) mm.add_module<LIBINT(bra, op, ket)>(key)

void load_modules(pluginplay::ModuleManager& mm) {
    LOAD_LIBINT(aos, op_base_type, aos, "Evaluate 2-Index BraKet (UQ)");
    LOAD_LIBINT(aos, op_base_type, aos_squared, "Evaluate 3-Index BraKet (UQ)");
    LOAD_LIBINT(aos_squared, op_base_type, aos_squared,
                "Evaluate 4-Index BraKet (UQ)");
    LOAD_LIBINT(aos, s_e_type, aos, "Overlap (UQ)");
    LOAD_LIBINT(aos, t_e_type, aos, "Kinetic (UQ)");
    LOAD_LIBINT(aos, v_en_type, aos, "Nuclear (UQ)");
    LOAD_LIBINT(aos, v_ee_type, aos, "ERI2 (UQ)");
    LOAD_LIBINT(aos, v_ee_type, aos_squared, "ERI3 (UQ)");
    LOAD_LIBINT(aos_squared, v_ee_type, aos_squared, "ERI4 (UQ)");
    set_defaults(mm);
}

#undef LOAD_LIBINT
#undef LIBINT

} // namespace integrals::libint::uq