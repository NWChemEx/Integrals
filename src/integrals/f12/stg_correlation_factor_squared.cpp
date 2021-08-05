// #include "f12.hpp"
// #include "integrals/types.hpp"
// #include <property_types/ao_integrals/type_traits.hpp>

// using namespace property_types::ao_integrals;

// static constexpr auto module_desc = R"""(

// For a Slater-type geminal the correlation factor, :math:`f_{12}(r_{12})`,
// is defined as:

// .. math::

//    f_{12}\left(r_{12}\right) = -\frac{1}{\gamma}e^{-\gamma r_{12}}

// and its square is thus given by:

// .. math::

//    f_{12}^2\left(r_{12}\right) = \frac{1}{\gamma^2}e^{-2\gamma r_{12}}

// This module calls a submodule to compute the integral of
// :math:`\exp(-2\gamma r_{12})` and then returns the result scaled
// appropriately.

// )""";

// namespace integrals::f12 {

// template<typename PropType>
// TEMPLATED_MODULE_CTOR(STGCorrelationFactorSquared, PropType) {
//     using element_type       = element_t<PropType>;
//     constexpr auto n_centers = n_centers_v<PropType>;
//     using n_center_type = NCenter<n_centers, type::ao_space_t<element_type>>;
//     using kernel_type   = STG<n_center_type>;

//     satisfies_property_type<PropType>();
//     description(module_desc);

//     add_submodule<kernel_type>("STG kernel");
//     add_input<element_type>("STG Exponent").set_default(1.2);
// }

// template<typename PropType>
// TEMPLATED_MODULE_RUN(STGCorrelationFactorSquared, PropType) {
//     using element_type       = element_t<PropType>;
//     using tensor_type        = type::tensor<element_type>;
//     constexpr auto n_centers = n_centers_v<PropType>;
//     using n_center_type = NCenter<n_centers, type::ao_space_t<element_type>>;
//     using kernel_type   = STG<n_center_type>;

//     auto gamma = inputs.at("STG Exponent").value<element_type>();
//     inputs.at("STG Exponent").change(2.0 * gamma);
//     auto kernel_output = submods.at("STG kernel").value().run(inputs);
//     auto [X]           = kernel_type::unwrap_results(kernel_output);
//     const auto c0      = -element_type{1.0} / gamma;
//     tensor_type c0X;
//     auto idx = TA::detail::dummy_annotation(n_centers);
//     c0X(idx) = c0 * c0 * X(idx);

//     auto rv = results();
//     return PropType::wrap_results(rv, c0X);
// }

// template<typename T>
// using cf2_2c = pt::correlation_factor_squared_2c<T>;

// template<typename T>
// using cf2_3c = pt::correlation_factor_squared_3c<T>;

// template<typename T>
// using cf2_4c = pt::correlation_factor_squared_4c<T>;

// template class STGCorrelationFactorSquared<cf2_2c<double>>;
// template class STGCorrelationFactorSquared<cf2_3c<double>>;
// template class STGCorrelationFactorSquared<cf2_4c<double>>;

// } // namespace integrals::f12
