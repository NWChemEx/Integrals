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

#include "../water_sto3g.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <integrals/integrals.hpp>

namespace test {

template<std::size_t N, typename FloatType = double>
auto eigen_buffer(const tensorwrapper::buffer::BufferBase& buffer) {
    return static_cast<const tensorwrapper::buffer::Eigen<FloatType, N>&>(
      buffer);
}

template<typename FloatType, unsigned short Rank>
auto trace(const tensorwrapper::buffer::Eigen<FloatType, Rank>& t) {
    Eigen::Tensor<FloatType, 0, Eigen::RowMajor> trace = t.value().trace();
    return trace.coeff();
}

template<typename FloatType, unsigned short Rank>
auto norm(const tensorwrapper::buffer::Eigen<FloatType, Rank>& t) {
    Eigen::Tensor<FloatType, 0, Eigen::RowMajor> norm =
      t.value().square().sum().sqrt();
    return norm.coeff();
}
} // namespace test