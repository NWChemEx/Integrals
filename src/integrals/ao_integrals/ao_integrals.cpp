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

#include "ao_integrals.hpp"
#include "detail_/get_basis_sets.hpp"
#include "detail_/libint_op.hpp"
#include "detail_/make_engine.hpp"
#include "detail_/make_libint_basis_set.hpp"
#include "detail_/shells2ord.hpp"
#include <type_traits>

namespace integrals::ao_integrals {

class LibIntVisitor : public chemist::qm_operator::OperatorVisitor {
public:
    using t_e_type  = simde::type::t_e_type;
    using v_ee_type = simde::type::v_ee_type;
    using v_en_type = simde::type::v_en_type;

    LibIntVisitor(const std::vector<libint2::BasisSet>& bases, double thresh,
                  std::size_t deriv = 0) :
      m_bases(bases), m_thresh(thresh), m_deriv(deriv) {};

    void run(const t_e_type& T_e) {
        m_engine = detail_::make_engine(m_bases, T_e, m_thresh, m_deriv);
    }

    void run(const v_en_type& V_en) {
        m_engine = detail_::make_engine(m_bases, V_en, m_thresh, m_deriv);
    }

    void run(const v_ee_type& V_ee) {
        m_engine = detail_::make_engine(m_bases, V_ee, m_thresh, m_deriv);
    }

    libint2::Engine& engine() { return m_engine; }

private:
    const std::vector<libint2::BasisSet>& m_bases;
    double m_thresh;
    std::size_t m_deriv;
    libint2::Engine m_engine;
};

template<typename BraKetType>
TEMPLATED_MODULE_CTOR(AOIntegral, BraKetType) {
    using my_pt = simde::EvaluateBraKet<BraKetType>;
    satisfies_property_type<my_pt>();
    description("Computes integrals with Libint");

    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");
}

template<typename BraKetType>
TEMPLATED_MODULE_RUN(AOIntegral, BraKetType) {
    using my_pt = simde::EvaluateBraKet<BraKetType>;

    const auto& [braket] = my_pt::unwrap_inputs(inputs);
    auto thresh          = inputs.at("Threshold").value<double>();
    auto bra             = braket.bra();
    auto ket             = braket.ket();
    auto op              = braket.op();

    // Gather information from Bra, Ket, and Op
    auto basis_sets = detail_::get_basis_sets(bra, ket);
    constexpr int N = detail_::get_n(bra, ket);

    // Dimensional information
    std::vector<std::size_t> dims_shells(N);
    Eigen::array<Eigen::Index, N> dims_bfs;
    for(auto i = 0; i < N; ++i) {
        dims_shells[i] = basis_sets[i].size();
        dims_bfs[i]    = basis_sets[i].nbf();
    }

    // Build tensor inputs
    using tensor_t = simde::type::tensor;
    using shape_t  = tensorwrapper::shape::Smooth;
    using layout_t = tensorwrapper::layout::Physical;
    using buffer_t = tensorwrapper::buffer::Eigen<double, N>;
    using data_t   = typename buffer_t::data_type;

    shape_t s{dims_bfs.begin(), dims_bfs.end()};
    layout_t l(s);
    data_t d(dims_bfs);
    buffer_t b{d, l};
    b.value().setZero();

    // Make libint engine
    LibIntVisitor visitor(basis_sets, thresh);
    op.visit(visitor);
    auto engine = visitor.engine();
    // auto engine     = detail_::make_engine(basis_sets, op, thresh);
    const auto& buf = engine.results();

    // Fill in values
    std::vector<std::size_t> shells(N, 0);
    while(shells[0] < dims_shells[0]) {
        detail_::run_engine_(engine, basis_sets, shells,
                             std::make_index_sequence<N>());

        auto vals = buf[0];
        if(vals) {
            auto ord   = detail_::shells2ord(basis_sets, shells);
            auto n_ord = ord.size();
            for(decltype(n_ord) i_ord = 0; i_ord < n_ord; ++i_ord) {
                b.value().data()[ord[i_ord]] = vals[i_ord];
            }
        }

        // Increment index
        shells[N - 1] += 1;
        for(auto i = 1; i < N; ++i) {
            if(shells[N - i] >= dims_shells[N - i]) {
                // Reset this dimension and increment the next one
                // shells[0] accumulates until we reach the end
                shells[N - i] = 0;
                shells[N - i - 1] += 1;
            }
        }
    }

    tensor_t t({s, b});
    auto rv = results();
    return my_pt::wrap_results(rv, t);
}

#define AOI(bra, op, ket) AOIntegral<braket<bra, op, ket>>
#define EXTERN_AOI(bra, op, ket) template struct AOI(bra, op, ket)
#define LOAD_AOI(bra, op, ket, key) mm.add_module<AOI(bra, op, ket)>(key)

// EXTERN_AOI(aos, op_type, aos);
// EXTERN_AOI(aos, op_type, aos_squared);
// EXTERN_AOI(aos_squared, op_type, aos_squared);

EXTERN_AOI(aos, t_e_type, aos);
EXTERN_AOI(aos, v_en_type, aos);
EXTERN_AOI(aos, v_ee_type, aos);
EXTERN_AOI(aos, v_ee_type, aos_squared);
EXTERN_AOI(aos_squared, v_ee_type, aos_squared);

void ao_integrals_set_defaults(pluginplay::ModuleManager& mm) {
    // Set any default associations
}

void load_ao_integrals(pluginplay::ModuleManager& mm) {
    // LOAD_AOI(aos, op_type, aos, "Evaluate 2-Index BraKet");
    // LOAD_AOI(aos, op_type, aos_squared, "Evaluate 3-Index BraKet");
    // LOAD_AOI(aos_squared, op_type, aos_squared, "Evaluate 4-Index BraKet");

    LOAD_AOI(aos, t_e_type, aos, "Kinetic");
    LOAD_AOI(aos, v_en_type, aos, "Nuclear");
    LOAD_AOI(aos, v_ee_type, aos, "ERI2");
    LOAD_AOI(aos, v_ee_type, aos_squared, "ERI3");
    LOAD_AOI(aos_squared, v_ee_type, aos_squared, "ERI4");
    ao_integrals_set_defaults(mm);
}

#undef AOI
#undef ADD_AOI

} // namespace integrals::ao_integrals