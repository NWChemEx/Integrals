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
#include <integrals/integrals.hpp>

using namespace tensorwrapper;

namespace integrals::ao_integrals {
namespace {

struct Kernel {
    using shape_type = buffer::Contiguous::shape_type;
    Kernel(shape_type shape) : m_shape(std::move(shape)) {}

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
        if constexpr(types::is_uncertain_v<float_type>) {
            throw std::runtime_error("Did not expect an uncertain type");
        } else {
            using uq_type  = sigma::Uncertain<float_type>;
            auto rv_buffer = buffer::make_contiguous<uq_type>(m_shape);
            auto rv_data   = buffer::get_raw_data<uq_type>(rv_buffer);
            for(std::size_t i = 0; i < t.size(); ++i) {
                const auto elem       = t[i];
                const auto elem_error = error[i];
                rv_data[i]            = uq_type(elem, elem_error);
            }
            rv = tensorwrapper::Tensor(m_shape, std::move(rv_buffer));
        }

        return rv;
    }
    shape_type m_shape;
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
}

MODULE_RUN(UQDriver) {
    const auto& [braket] = eri_pt::unwrap_inputs(inputs);

    auto& eri_mod = submods.at("ERIs").value();
    auto tol      = eri_mod.inputs().at("Threshold").value<double>();

    const auto& t     = eri_mod.run_as<eri_pt>(braket);
    const auto& error = submods.at("ERI Error").run_as<error_pt>(braket, tol);

    using buffer::visit_contiguous_buffer;
    shape::Smooth shape = t.buffer().layout().shape().as_smooth().make_smooth();
    Kernel k(shape);
    auto t_buffer     = make_contiguous(t.buffer());
    auto error_buffer = make_contiguous(error.buffer());
    auto t_w_error    = visit_contiguous_buffer(k, t_buffer, error_buffer);

    auto rv = results();
    return eri_pt::wrap_results(rv, t_w_error);
}
} // namespace integrals::ao_integrals
