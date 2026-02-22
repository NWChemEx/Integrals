/* Copyright 2024 NWChemEx - Project
 * Licensed under the Apache License, Version 2.0(the "License");
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
#pragma once
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <integrals/integrals.hpp>
#include <simde/simde.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

#include "ao_bases.hpp"
#include "molecules.hpp"
#include "shell_quartets.hpp"

namespace integrals::testing {

inline auto initialize_integrals() {
    auto mm = pluginplay::ModuleManager();
    integrals::load_modules(mm);
    integrals::set_defaults(mm);
    return mm;
}

template<typename FloatType, std::size_t N, std::size_t... Is>
auto eigen_tensor_(const tensorwrapper::buffer::BufferBase& buffer,
                   std::array<int, N> extents, std::index_sequence<Is...>) {
    using namespace tensorwrapper;
    const auto pdata = buffer::get_raw_data<FloatType>(buffer);
    using eigen_type = Eigen::Tensor<const FloatType, N, Eigen::RowMajor>;
    return Eigen::TensorMap<eigen_type>(pdata.data(), extents[Is]...);
}

// Checking eigen outputs
template<std::size_t N, typename FloatType = double>
auto eigen_tensor(const tensorwrapper::buffer::BufferBase& buffer) {
    std::array<int, N> extents;
    auto shape = buffer.layout().shape().as_smooth();
    for(std::size_t i = 0; i < N; ++i) extents[i] = shape.extent(i);
    return eigen_tensor_<FloatType>(buffer, extents,
                                    std::make_index_sequence<N>());
}

template<unsigned short Rank, typename FloatType = double>
auto trace(const tensorwrapper::buffer::BufferBase& buffer) {
    auto t = eigen_tensor<Rank, FloatType>(buffer);
    Eigen::Tensor<FloatType, 0, Eigen::RowMajor> trace = t.trace();
    return trace.coeff();
}

template<unsigned short Rank, typename FloatType = double>
auto norm(const tensorwrapper::buffer::BufferBase& buffer) {
    auto t = eigen_tensor<Rank, FloatType>(buffer);
    Eigen::Tensor<FloatType, 0, Eigen::RowMajor> norm = t.square().sum().sqrt();
    return norm.coeff();
}

inline auto h2_mos() {
    using mos_type    = simde::type::mos;
    using tensor_type = typename mos_type::transform_type;
    tensor_type c({{-0.565516, -1.07019}, {-0.565516, 1.07019}});
    return mos_type(simde::type::aos(h2_sto3g_basis_set()), std::move(c));
}

inline auto h2_density() {
    using density_type = simde::type::decomposable_e_density;
    typename density_type::value_type rho(
      {{0.31980835, 0.31980835}, {0.31980835, 0.31980835}});
    return density_type(rho, h2_mos());
}

} // namespace integrals::testing