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

namespace integrals::ao_integrals {
namespace {

struct Kernel {
    using buffer_base_type = tensorwrapper::buffer::BufferBase;
    template<typename FloatType>
    auto run(const buffer_base_type& t, const buffer_base_type& error) {
        tensorwrapper::Tensor rv;

        if constexpr(tensorwrapper::types::is_uncertain_v<FloatType>) {
            using alloc_type = tensorwrapper::allocator::Eigen<FloatType>;
            alloc_type alloc(t.allocator().runtime());

            const auto& t_eigen     = alloc.rebind(t);
            const auto& error_eigen = alloc.rebind(error);

            auto rv_buffer = alloc.allocate(t_eigen.layout());
            for(std::size_t i = 0; i < t_eigen.size(); ++i) {
                const auto elem          = (t_eigen.data() + i)->mean();
                const auto elem_error    = (error_eigen.data() + i)->mean();
                *(rv_buffer->data() + i) = FloatType(elem, elem_error);
            }

            const auto& shape = t_eigen.layout().shape();
            rv = tensorwrapper::Tensor(shape, std::move(rv_buffer));
        } else {
            throw std::runtime_error("Expected an uncertain type");
        }
        return rv;
    }
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

    using tensorwrapper::utilities::floating_point_dispatch;
    Kernel k;
    auto t_w_error = floating_point_dispatch(k, t.buffer(), error.buffer());

    auto rv = results();
    return eri_pt::wrap_results(rv, t_w_error);
}
} // namespace integrals::ao_integrals