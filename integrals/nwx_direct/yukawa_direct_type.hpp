#pragma once
#include "integrals/nwx_TA/fill_ND_functor.hpp"
#include "integrals/nwx_direct/direct_tile.hpp"
#include <property_types/types.hpp>
#include <sde/property_type.hpp>

namespace property_types {

/**
 * @brief The property type for modules that build tensors filled with
 * Yukawa (STG times Coulomb) integrals in the AO basis set.
 *
 * @tparam ElementType The type of the element in the tensor. Defaults to
 *                     `double`.
 */
template<typename ElementType = double>
struct Yukawa3CDirect : public sde::PropertyType<Yukawa3CDirect<ElementType>> {
    /// The type of an std::array of basis sets
    using basis_type = type::basis_set<ElementType>;
    /// The type of a direct tile builder
    using builder_type = nwx_TA::FillNDFunctor<TA::Tensor<ElementType>,
                                               libint2::Operator::yukawa, 3>;
    /// The type of a direct tile
    using tile_type = DirectTile<TA::Tensor<ElementType>, builder_type>;
    /// The type of a tensor accounting for ElementType
    using tensor_type = TA::DistArray<tile_type, TA::SparsePolicy>;
    /// Generates the input fields required by this property type
    auto inputs_();
    /// Generates the result fields required by this property type
    auto results_();
}; // class Yukawa3CIntegral

template<typename ElementType = double>
struct Yukawa4CDirect : public sde::PropertyType<Yukawa4CDirect<ElementType>> {
    /// The type of an std::array of basis sets
    using basis_type = type::basis_set<ElementType>;
    /// The type of a direct tile builder
    using builder_type = nwx_TA::FillNDFunctor<TA::Tensor<ElementType>,
                                               libint2::Operator::yukawa, 4>;
    /// The type of a direct tile
    using tile_type = DirectTile<TA::Tensor<ElementType>, builder_type>;
    /// The type of a tensor accounting for ElementType
    using tensor_type = TA::DistArray<tile_type, TA::SparsePolicy>;
    /// Generates the input fields required by this property type
    auto inputs_();
    /// Generates the result fields required by this property type
    auto results_();
}; // class Yukawa4CIntegral

//------------------------Implementations---------------------------------------

template<typename ElementType>
auto Yukawa3CDirect<ElementType>::inputs_() {
    auto rv =
      sde::declare_input()
        .add_field<const basis_type&>("Bra")
        .template add_field<const basis_type&>("Ket1")
        .template add_field<const basis_type&>("Ket2")
        .template add_field<type::size>("Derivative", type::size{0})
        .template add_field<ElementType>("STG Exponent", ElementType{1.0});
    rv["Bra"].set_description("The basis set for the bra");
    rv["Ket1"].set_description("The first basis set for the ket");
    rv["Ket2"].set_description("The second basis set for the ket");
    rv["Derivative"].set_description(
      "The derivative order of the integrals to be computed");
    rv["STG Exponent"].set_description(
      "The exponent for the Slate type geminal");
    return rv;
}

template<typename ElementType>
auto Yukawa3CDirect<ElementType>::results_() {
    auto rv = sde::declare_result().add_field<tensor_type>("Yukawa Integrals");
    rv["Yukawa Integrals"].set_description("The requested Yukawa integrals");
    return rv;
}

template<typename ElementType>
auto Yukawa4CDirect<ElementType>::inputs_() {
    auto rv =
      sde::declare_input()
        .add_field<const basis_type&>("Bra1")
        .template add_field<const basis_type&>("Bra2")
        .template add_field<const basis_type&>("Ket1")
        .template add_field<const basis_type&>("Ket2")
        .template add_field<type::size>("Derivative", type::size{0})
        .template add_field<ElementType>("STG Exponent", ElementType{1.0});
    rv["Bra1"].set_description("The first basis set for the bra");
    rv["Bra2"].set_description("The second basis set for the bra");
    rv["Ket1"].set_description("The first basis set for the ket");
    rv["Ket2"].set_description("The second basis set for the ket");
    rv["Derivative"].set_description(
      "The derivative order of the integrals to be computed");
    rv["STG Exponent"].set_description(
      "The exponent for the Slate type geminal");
    return rv;
}

template<typename ElementType>
auto Yukawa4CDirect<ElementType>::results_() {
    auto rv = sde::declare_result().add_field<tensor_type>("Yukawa Integrals");
    rv["Yukawa Integrals"].set_description("The requested Yukawa integrals");
    return rv;
}

extern template class Yukawa3CDirect<double>;
extern template class Yukawa3CDirect<float>;
extern template class Yukawa4CDirect<double>;
extern template class Yukawa4CDirect<float>;

} // namespace property_types