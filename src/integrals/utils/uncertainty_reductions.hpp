/*
 * Copyright 2026 NWChemEx-Project
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

#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace integrals::utils {

enum class mean_type { none, max, geometric, harmonic };

template<typename ContainerType>
auto max_mean(const ContainerType& values) {
    using float_type = std::decay_t<decltype(values[0])>;
    if(values.empty()) return float_type{0};
    return *std::max_element(values.begin(), values.end());
}

template<typename ContainerType>
auto geometric_mean(const ContainerType& values) {
    using float_type           = std::decay_t<decltype(values[0])>;
    float_type log_sum         = 0.0;
    std::size_t non_zero_count = 0;
    for(const auto& val : values) {
        if(val == 0) continue;
        log_sum += std::log(std::fabs(val));
        ++non_zero_count;
    }
    if(non_zero_count == 0 || log_sum == 0) return float_type{0};
    return std::exp(log_sum / non_zero_count);
}

template<typename ContainerType>
auto harmonic_mean(const ContainerType& values) {
    // N.b. If val is very small then 1 / val can be very large and be infinity.
    // Then non_zero_count / reciprocal_sum can be NaN.
    using float_type           = std::decay_t<decltype(values[0])>;
    float_type reciprocal_sum  = 0.0;
    std::size_t non_zero_count = 0;
    for(const auto& val : values) {
        if(val == 0 || std::isnan(1 / val)) continue;
        reciprocal_sum += 1 / val;
        ++non_zero_count;
    }
    if(non_zero_count == 0 || reciprocal_sum == 0) return float_type{0};
    if(std::isnan(non_zero_count / reciprocal_sum)) return float_type{0};
    return non_zero_count / reciprocal_sum;
}

template<typename ContainerType>
auto compute_mean(mean_type mean, const ContainerType& span) {
    if(mean == mean_type::max) {
        return max_mean(span);
    } else if(mean == mean_type::harmonic) {
        return harmonic_mean(span);
    } else if(mean == mean_type::geometric) {
        return geometric_mean(span);
    } else {
        throw std::runtime_error("Mean type not supported");
    }
}

inline auto mean_from_string(const std::string& mean_str) {
    if(mean_str == std::string("max")) {
        return mean_type::max;
    } else if(mean_str == std::string("geometric")) {
        return mean_type::geometric;
    } else if(mean_str == std::string("harmonic")) {
        return mean_type::harmonic;
    } else if(mean_str == std::string("none")) {
        return mean_type::none;
    } else {
        throw std::runtime_error("Mean type not supported: " + mean_str);
    }
}

} // namespace integrals::utils
