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

struct Kernel {
    template<typename FloatType>
    auto run(const tensorwrapper::buffer::BufferBase& I,
             parallelzone::runtime::RuntimeView& rv) {
        // Unwrap buffer info
        tensorwrapper::allocator::Eigen<FloatType> allocator(rv);
        const auto& eigen_I = allocator.rebind(I);
        const auto* pI      = eigen_I.get_immutable_data();
        const auto& shape_I = eigen_I.layout().shape().as_smooth();
        auto rows           = shape_I.extent(0);
        auto cols           = shape_I.extent(1);

        // Cholesky Decomp
        constexpr auto rmajor = Eigen::RowMajor;
        constexpr auto edynam = Eigen::Dynamic;
        using matrix_type = Eigen::Matrix<FloatType, edynam, edynam, rmajor>;
        using map_type    = Eigen::Map<const matrix_type>;
        map_type I_map(pI, rows, cols);
        Eigen::LLT<matrix_type> lltOfI(I_map);
        matrix_type U    = lltOfI.matrixU();
        matrix_type Linv = U.inverse().transpose();

        // Wrap result
        tensorwrapper::shape::Smooth matrix_shape{rows, cols};
        tensorwrapper::layout::Physical matrix_layout(matrix_shape);
        auto pM_buffer = allocator.allocate(matrix_layout);
        for(decltype(rows) i = 0; i < rows; ++i) {
            for(decltype(cols) j = 0; j < cols; ++j) {
                pM_buffer->set_elem({i, j}, Linv(i, j));
            }
        }
        return simde::type::tensor(matrix_shape, std::move(pM_buffer));
    }
};

auto desc = R"(
Inverse Coulomb Metric
---------------------
)";

} // namespace

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

    // Compute metric
    using tensorwrapper::utilities::floating_point_dispatch;
    auto r = get_runtime();
    Kernel k;
    const auto& I_buffer = I.buffer();
    auto M               = floating_point_dispatch(k, I_buffer, r);

    auto rv = results();
    return pt::wrap_results(rv, M);
}

} // namespace integrals::ao_integrals
