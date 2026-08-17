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

/** @file taylor_model_factory.hpp
 *
 * Defines integrals::property_types::TaylorModelFactory, the result type of
 * the UQInitializer property type (declared in property_types.hpp).
 */
#pragma once
#include <cstddef>
#include <tensorwrapper/tensorwrapper.hpp>

namespace integrals::property_types {

/** @brief A copyable, comparable factory for constructing Taylor-model UQ
 *         values at a caller-configured truncation order.
 *
 *  sigma::TaylorModel's truncation order can only be set at construction
 *  time (there is no in-place mutator, and sweep_to_order can only lower an
 *  already-built model's order). UQAtomSymmBlockedDriver constructs one UQ
 *  scalar per ERI tensor element, so routing every construction through a
 *  PluginPlay module call would be prohibitively slow; instead the
 *  UQInitializer module is run once per UQAtomSymmBlockedDriver::run() to
 *  produce one of these factories, which is then called directly in the
 *  per-element hot loop. PluginPlay's AnyField type erasure requires results
 *  to be copyable, equality-comparable, and less-than-comparable (a raw
 *  std::function does not satisfy this), so the factory is this small
 *  hand-rolled value type rather than a std::function.
 */
class TaylorModelFactory {
public:
    using order_type = std::size_t;

    TaylorModelFactory(order_type order = 2) : m_order_(order) {}

    // Templated on the underlying floating-point type so this factory can be
    // used generically wherever tensorwrapper::types::taylor_model_type<T> is
    // instantiated (e.g. UQAtomSymmBlockedDriver's Kernel is instantiated for
    // both float and double element types, even though real ERI tensors are
    // always double).
    template<typename T>
    tensorwrapper::types::taylor_model_type<T> operator()(T center,
                                                          T radius) const {
        using uq_t = tensorwrapper::types::taylor_model_type<T>;
        return uq_t(center - radius, center + radius,
                    typename uq_t::Order(m_order_));
    }

    order_type order() const noexcept { return m_order_; }

    bool operator==(const TaylorModelFactory& rhs) const noexcept {
        return m_order_ == rhs.m_order_;
    }
    bool operator<(const TaylorModelFactory& rhs) const noexcept {
        return m_order_ < rhs.m_order_;
    }

private:
    order_type m_order_;
};

} // namespace integrals::property_types
