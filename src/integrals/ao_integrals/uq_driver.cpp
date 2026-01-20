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

        if constexpr(types::is_uncertain_v<FloatType>) {
            auto rv_buffer = buffer::make_contiguous<FloatType>(m_shape);
            auto rv_data   = buffer::get_raw_data<FloatType>(rv_buffer);
            for(std::size_t i = 0; i < t.size(); ++i) {
                const auto elem       = t[i].mean();
                const auto elem_error = error[i].mean();
                rv_data[i]            = FloatType(elem, elem_error);
            }

            rv = tensorwrapper::Tensor(m_shape, std::move(rv_buffer));
        } else {
            throw std::runtime_error("Expected an uncertain type");
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

using eri_pt = simde::ERI4;

MODULE_CTOR(UQDriver) {
    satisfies_property_type<eri_pt>();
    description(desc);
    add_submodule<eri_pt>("ERIs");
    add_input<double>("benchmark precision").set_default(1.0e-16);
    add_input<double>("precision").set_default(1.0e-16);
}

MODULE_RUN(UQDriver) {
    auto tau_0 = inputs.at("benchmark precision").value<double>();
    auto tau   = inputs.at("precision").value<double>();

    auto& eri_mod = submods.at("ERIs").value();

    auto benchmark_mod = eri_mod.unlocked_copy();
    benchmark_mod.change_input("Threshold", tau_0);
    benchmark_mod.change_input("With UQ?", true);

    auto normal_mod = eri_mod.unlocked_copy();
    normal_mod.change_input("Threshold", tau);
    normal_mod.change_input("With UQ?", true);

    const auto& [t_0] = eri_pt::unwrap_results(benchmark_mod.run(inputs));
    const auto& [t]   = eri_pt::unwrap_results(normal_mod.run(inputs));

    simde::type::tensor error;
    error("m,n,l,s") = t("m,n,l,s") - t_0("m,n,l,s");

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
