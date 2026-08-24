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

/** @file uq_factory.hpp
 *
 * Defines integrals::property_types::UQFactory, the result type of the
 * UQInitializer property type (declared in property_types.hpp).
 */
#pragma once
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <tensorwrapper/tensorwrapper.hpp>

namespace integrals::property_types {

/// The UQ representations UQInitializer knows how to build a factory for.
enum class UQKind {
    uncertain,
    interval,
    affine,
    thresholded_affine,
    taylor_model
};

/** @brief Parses a "UQ Type" module-input string into a UQKind.
 *
 *  This is the single place the accepted "UQ Type" string literals are
 *  spelled out; UQInitializer and UQAtomSymmBlockedDriver's tests both rely
 *  on this mapping staying in sync with UQFactory's construction logic.
 *
 *  @throw std::runtime_error if @p uq_type does not match a known UQ kind.
 */
inline UQKind uq_kind_from_string(const std::string& uq_type) {
    if(uq_type == "uncertain") return UQKind::uncertain;
    if(uq_type == "interval") return UQKind::interval;
    if(uq_type == "affine") return UQKind::affine;
    if(uq_type == "thresholded affine") return UQKind::thresholded_affine;
    if(uq_type == "taylor model") return UQKind::taylor_model;
    throw std::runtime_error(
      "integrals::property_types::uq_kind_from_string: Invalid UQ type name " +
      uq_type);
}

/// Inverse of uq_kind_from_string, kept alongside it for symmetry/debugging.
inline std::string to_string(UQKind kind) {
    switch(kind) {
        case UQKind::uncertain: return "uncertain";
        case UQKind::interval: return "interval";
        case UQKind::affine: return "affine";
        case UQKind::thresholded_affine: return "thresholded affine";
        case UQKind::taylor_model: return "taylor model";
    }
    throw std::runtime_error(
      "integrals::property_types::to_string(UQKind): Unrecognized UQKind");
}

/** @brief Polymorphic base class for constructing a UQ scalar from a center
 *         and a radius.
 *
 *  The UQFactorBase class is used as a PIMPL for UQFactory (below) so that
 *
 *  Each UQ representation (sigma::Uncertain, sigma::Interval, sigma::Affine,
 *  sigma::ThresholdedAffine, sigma::TaylorModel) is an independent template
 *  class with no shared base in `sigma`/`tensorwrapper` (they're unified only
 *  by compile-time traits). This class supplies the runtime-polymorphic base
 *  UQInitializer needs so it can hand back "a factory for whichever UQ kind
 *  the caller asked for" as a single result type, and so downstream code
 *  (UQAtomSymmBlockedDriver's Kernel) can construct UQ values without a
 *  compile-time dependence on which kind was selected.
 *
 *  N.b., only `double` centers/radii are supported at present. If/when needed
 *  the API could be extended to other types.
 */
class UQFactoryBase {
public:
    /// Virtual dtor
    virtual ~UQFactoryBase() = default;

    /// Constructs a UQ scalar of *this's kind, centered at @p center with
    /// radius @p radius, type-erased into a wtf::fp::Float.
    virtual wtf::fp::Float operator()(double center, double radius) const = 0;

    /// Identifies which UQ representation *this constructs.
    virtual UQKind kind() const noexcept = 0;

    /// Deep-copies *this. Kept for forward compatibility even though
    /// UQFactory (below) uses shared, not owning, copies.
    virtual std::unique_ptr<UQFactoryBase> clone() const = 0;

    /// True if *this and @p rhs would construct value-equal UQ scalars.
    virtual bool equal(const UQFactoryBase& rhs) const noexcept {
        return kind() == rhs.kind();
    }

    /// Imposes a strict weak ordering, for PluginPlay's AnyField.
    virtual bool less(const UQFactoryBase& rhs) const noexcept {
        return kind() < rhs.kind();
    }
};

namespace detail {

/// Implements UQFactory for all UQ kinds except TaylorModel.
template<UQKind Kind, template<typename> typename UQType>
class SimpleFactoryImpl final : public UQFactoryBase {
public:
    /// Dispatches to tensorwrapper::types::construct_uq_type<UQType>
    wtf::fp::Float operator()(double center, double radius) const override {
        return wtf::fp::Float(
          tensorwrapper::types::construct_uq_type<UQType<double>>(center,
                                                                  radius));
    }

    /// Returns the UQKind template parameter.
    UQKind kind() const noexcept override { return Kind; }

    /// Returns a deep-copy of *this.
    std::unique_ptr<UQFactoryBase> clone() const override {
        return std::make_unique<SimpleFactoryImpl>(*this);
    }
};

/// Convenience aliases for the 4 SimpleFactoryImpl specializations.
using UncertainFactoryImpl =
  SimpleFactoryImpl<UQKind::uncertain, tensorwrapper::types::uncertain_type>;
using IntervalFactoryImpl =
  SimpleFactoryImpl<UQKind::interval, tensorwrapper::types::interval_type>;
using AffineFactoryImpl =
  SimpleFactoryImpl<UQKind::affine, tensorwrapper::types::affine_type>;
using ThresholdedAffineFactoryImpl =
  SimpleFactoryImpl<UQKind::thresholded_affine,
                    tensorwrapper::types::thresholded_affine_type>;

/** @brief Constructs a TaylorModel scalar spanning [center - radius,
 *         center + radius] with a maximum order of @p order.
 *
 *  When Sigma is disabled tensorwrapper::types::taylor_model_type<T> is just
 *  @p T, which has no order to set, so we defer to TensorWrapper's generic
 *  construction (which simply returns a scalar carrying no uncertainty). The
 *  dispatch lives in a function template because a discarded `if constexpr`
 *  branch is only left uninstantiated inside a template.
 */
template<typename T>
T make_taylor_model(double center, double radius, std::size_t order) {
    if constexpr(tensorwrapper::types::is_taylor_model_v<T>) {
        return T(center - radius, center + radius, typename T::Order(order));
    } else {
        return tensorwrapper::types::construct_uq_type<T>(center, radius);
    }
}

/** @brief Specializes UQFactoryBase to construct TaylorModel UQ scalars.
 *
 *  TaylorModel objects need an order parameter. This class stores that
 *  parameter and passes it to each TaylorModel upon construction.
 */
class TaylorModelFactoryImpl final : public UQFactoryBase {
public:
    /// Takes the TaylorModel order to use when constructing UQ scalars.
    explicit TaylorModelFactoryImpl(std::size_t order = 2) : m_order_(order) {}

    /// Calls ctor for tensorwrapper::types::taylor_model_type<double>
    wtf::fp::Float operator()(double center, double radius) const override {
        using uq_t = tensorwrapper::types::taylor_model_type<double>;
        return wtf::fp::Float(
          make_taylor_model<uq_t>(center, radius, m_order_));
    }

    /// Returns UQKind::taylor_model.
    UQKind kind() const noexcept override { return UQKind::taylor_model; }

    /// Returns the TaylorModel order stored in *this.
    std::size_t order() const noexcept { return m_order_; }

    // Deep-copy *this
    std::unique_ptr<UQFactoryBase> clone() const override {
        return std::make_unique<TaylorModelFactoryImpl>(*this);
    }

    /// Additionally compares the order
    bool equal(const UQFactoryBase& rhs) const noexcept override {
        if(rhs.kind() != UQKind::taylor_model) return false;
        return m_order_ ==
               static_cast<const TaylorModelFactoryImpl&>(rhs).m_order_;
    }

    bool less(const UQFactoryBase& rhs) const noexcept override {
        if(kind() != rhs.kind()) return kind() < rhs.kind();
        return m_order_ <
               static_cast<const TaylorModelFactoryImpl&>(rhs).m_order_;
    }

private:
    std::size_t m_order_;
};

} // namespace detail

/** @brief A copyable, comparable factory for constructing UQ values of a
 *         caller-configured kind.
 *
 *  This is the public API for creating UQ Forms without needing to know what
 *  kind of UQ Form is being created
 *
 *  This class wraps one in a std::shared_ptr<const UQFactoryBase>: copies are
 *  shallow (the underlying Impl is immutable after construction), giving
 *  cheap value semantics without deep-cloning on every copy.
 */
class UQFactory {
public:
    /// Defaults to an "uncertain" factory so PluginPlay's requirement that
    /// results be default constructible is met without a null state.
    UQFactory() : m_impl_(std::make_shared<detail::UncertainFactoryImpl>()) {}

    /// Constructs a factory for the UQ kind specified by @p kind, with optional
    /// TaylorModel order @p order (ignored for non-TaylorModel kinds).
    explicit UQFactory(UQKind kind, std::size_t order = 2) :
      m_impl_(make_impl_(kind, order)) {}

    /// Main API for building the UQ value
    wtf::fp::Float operator()(double center, double radius) const {
        return (*m_impl_)(center, radius);
    }

    /// Returns the UQKind *this will build
    UQKind kind() const noexcept { return m_impl_->kind(); }

    /// True if *this and @p rhs would build value-equal UQ scalars.
    bool operator==(const UQFactory& rhs) const noexcept {
        return m_impl_->equal(*rhs.m_impl_);
    }

    /// Imposes a strict weak ordering
    bool operator<(const UQFactory& rhs) const noexcept {
        return m_impl_->less(*rhs.m_impl_);
    }

private:
    static std::shared_ptr<const UQFactoryBase> make_impl_(UQKind kind,
                                                           std::size_t order) {
        switch(kind) {
            case UQKind::uncertain:
                return std::make_shared<detail::UncertainFactoryImpl>();
            case UQKind::interval:
                return std::make_shared<detail::IntervalFactoryImpl>();
            case UQKind::affine:
                return std::make_shared<detail::AffineFactoryImpl>();
            case UQKind::thresholded_affine:
                return std::make_shared<detail::ThresholdedAffineFactoryImpl>();
            case UQKind::taylor_model:
                return std::make_shared<detail::TaylorModelFactoryImpl>(order);
        }
        throw std::runtime_error(
          "integrals::property_types::UQFactory: Unrecognized UQKind");
    }

    std::shared_ptr<const UQFactoryBase> m_impl_;
};

} // namespace integrals::property_types
