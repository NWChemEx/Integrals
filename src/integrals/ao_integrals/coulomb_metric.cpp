/*
 * Copyright 2025 NWChemEx-Project
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
#include <Eigen/Eigen>

namespace integrals::ao_integrals {

using pt = simde::ERI2;

namespace {

auto desc = R"(
Inverse Coulomb Metric
---------------------
)";

}

MODULE_CTOR(CoulombMetric) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<pt>("Two-center ERI");
}

MODULE_RUN(CoulombMetric) {
    const auto& [braket] = pt::unwrap_inputs(inputs);
    auto& eri2_mod       = submods.at("Two-center ERI");

    // Get ERI2
    const auto& I = eri2_mod.run_as<pt>(braket);

    // Cholesky Decomp
    tensorwrapper::allocator::Eigen<double> allocator(get_runtime());
    const auto& eigen_I = allocator.rebind(I.buffer());
    const auto* pI      = eigen_I.data();
    const auto& shape_I = eigen_I.layout().shape().as_smooth();
    auto rows           = shape_I.extent(0);
    auto cols           = shape_I.extent(1);

    using emat_t = Eigen::MatrixXd;
    Eigen::Map<const emat_t> I_map(pI, rows, cols);
    Eigen::LLT<emat_t> lltOfI(I_map);
    emat_t U    = lltOfI.matrixU();
    emat_t Linv = U.inverse().transpose();

    // Wrap result
    tensorwrapper::shape::Smooth matrix_shape{rows, cols};
    tensorwrapper::layout::Physical matrix_layout(matrix_shape);
    auto pM_buffer = allocator.allocate(matrix_layout);

    for(auto i = 0; i < rows; ++i) {
        for(auto j = 0; j < cols; ++j) { pM_buffer->at(i, j) = Linv(i, j); }
    }
    simde::type::tensor M(matrix_shape, std::move(pM_buffer));

    auto rv = results();
    return pt::wrap_results(rv, M);
}

} // namespace integrals::ao_integrals