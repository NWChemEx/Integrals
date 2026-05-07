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

#include "../utils/uncertainty_reductions.hpp"
#include "ao_integrals.hpp"
#include <integrals/integrals.hpp>
#ifdef ENABLE_SIGMA
#include <sigma/sigma.hpp>
#endif

using namespace tensorwrapper;

namespace integrals::ao_integrals {
namespace {

template<template<typename> typename UQType>
struct Kernel {
    using shape_type = buffer::Contiguous::shape_type;
    Kernel(shape_type shape, utils::mean_type mean) :
      m_shape(std::move(shape)), m_mean(mean) {}

    template<typename FloatType0, typename FloatType1>
    Tensor operator()(const std::span<FloatType0> t,
                      const std::span<FloatType1> error) {
        throw std::runtime_error(
          "UQ Integrals Driver kernel only supports same float types");
    }

    template<typename FloatType>
    auto operator()(const std::span<FloatType> t,
                    const std::span<FloatType> error) {
        Tensor rv;

        using float_type = std::decay_t<FloatType>;
        if constexpr(types::is_uncertain_v<float_type> ||
                     types::is_interval_v<float_type>) {
            throw std::runtime_error("Did not expect an uncertain type");
        } else {
#ifdef ENABLE_SIGMA
            using uq_type  = UQType<float_type>;
            auto rv_buffer = buffer::make_contiguous<uq_type>(m_shape);
            auto rv_data   = buffer::get_raw_data<uq_type>(rv_buffer);

            bool use_mean = (m_mean != utils::mean_type::none);
            if(!use_mean) {
                for(std::size_t i = 0; i < t.size(); ++i) {
                    const auto elem = t[i];
                    if constexpr(types::is_interval_v<uq_type>) {
                        auto ei    = std::fabs(error[i]);
                        rv_data[i] = uq_type{elem - ei, elem + ei};
                    } else if constexpr(types::is_uncertain_v<uq_type>) {
                        rv_data[i] = uq_type{elem, error[i]};
                    } else {
                        throw std::runtime_error("Invalid UQ type");
                    }
                }
            } else {
                float_type max_error = utils::compute_mean(m_mean, error);
                uq_type max_uq;
                if constexpr(types::is_interval_v<uq_type>) {
                    max_uq = uq_type(-max_error, max_error);
                } else if constexpr(types::is_uncertain_v<uq_type>) {
                    max_uq = uq_type(0.0, max_error);
                } else {
                    throw std::runtime_error("Invalid UQ type");
                }
                for(std::size_t i = 0; i < t.size(); ++i) {
                    const auto elem = t[i];
                    rv_data[i]      = uq_type(elem) + max_uq;
                }
            }
            rv = tensorwrapper::Tensor(m_shape, std::move(rv_buffer));
#else
            throw std::runtime_error("Sigma support not enabled!");
#endif
        }

        return rv;
    }
    shape_type m_shape;
    utils::mean_type m_mean;
};

const auto desc = R"(
UQ Integrals Driver
-------------------

)";

} // namespace

using eri_pt   = simde::ERI4;
using error_pt = integrals::property_types::Uncertainty<eri_pt>;

MODULE_CTOR(UQDriver) {
    satisfies_property_type<eri_pt>();
    description(desc);
    add_submodule<eri_pt>("ERIs");
    add_submodule<error_pt>("ERI Error");
    add_input<std::string>("UQ Type").set_default("uncertain");
    add_input<std::string>("Mean Type").set_default("none");
}

MODULE_RUN(UQDriver) {
    const auto& [braket] = eri_pt::unwrap_inputs(inputs);
    auto uq_type         = inputs.at("UQ Type").value<std::string>();
    auto mean_str        = inputs.at("Mean Type").value<std::string>();
    auto mean            = utils::mean_from_string(mean_str);

    auto& eri_mod = submods.at("ERIs").value();
    auto tol      = eri_mod.inputs().at("Threshold").value<double>();

    const auto& t     = eri_mod.run_as<eri_pt>(braket);
    const auto& error = submods.at("ERI Error").run_as<error_pt>(braket, tol);

    using buffer::visit_contiguous_buffer;
    shape::Smooth shape = t.buffer().layout().shape().as_smooth().make_smooth();
    auto t_buffer       = make_contiguous(t.buffer());
    auto error_buffer   = make_contiguous(error.buffer());

    simde::type::tensor t_w_error;
    if(uq_type == "uncertain") {
        Kernel<tensorwrapper::types::uncertain_type> k(shape, mean);
        t_w_error = visit_contiguous_buffer(k, t_buffer, error_buffer);
    } else if(uq_type == "interval") {
        Kernel<tensorwrapper::types::interval_type> k(shape, mean);
        t_w_error = visit_contiguous_buffer(k, t_buffer, error_buffer);
    } else {
        throw std::runtime_error("Invalid UQ type");
    }
    auto rv = results();
    return eri_pt::wrap_results(rv, t_w_error);
}
} // namespace integrals::ao_integrals
