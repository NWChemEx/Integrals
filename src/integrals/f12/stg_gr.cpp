// #include "f12.hpp"
// #include "integrals/types.hpp"
// #include <property_types/ao_integrals/type_traits.hpp>

// using namespace property_types::ao_integrals;

// static constexpr auto module_desc = R"""(
// This module computes the quantity :math:`\frac{f_{12}(r_{12})}{r_{12}}` for a
// Slater-type geminal correlation factor according to:

// .. math::

//    \newcommand{\r_one_two}{r_{12}}
//    \newcommand{\f_one_two}{f_{12}\left(\r_one_two\right)}
//    \newcommand{\stg}{e^{-\gamma \r_one_two}}

//    \frac{\f_one_two}{\r_one_two} = \frac{-\stg}{\gamma \r_one_two}

// In practice this is computed by scaling the Yukawa integral by
// :math:`-\frac{1}{\gamma}`.
// )""";

// namespace integrals::f12 {

// template<typename PropType>
// TEMPLATED_MODULE_CTOR(STGGR, PropType) {
//     using element_type       = element_t<PropType>;
//     constexpr auto n_centers = n_centers_v<PropType>;
//     using n_center_type = NCenter<n_centers, type::ao_space_t<element_type>>;
//     using kernel_type   = Yukawa<n_center_type>;

//     satisfies_property_type<PropType>();
//     description(module_desc);

//     add_submodule<kernel_type>("Yukawa kernel");
//     add_input<element_type>("STG Exponent").set_default(1.2);
// }

// template<typename PropType>
// TEMPLATED_MODULE_RUN(STGGR, PropType) {
//     using element_type       = element_t<PropType>;
//     using tensor_type        = type::tensor<element_type>;
//     constexpr auto n_centers = n_centers_v<PropType>;
//     using n_center_type = NCenter<n_centers, type::ao_space_t<element_type>>;
//     using kernel_type   = Yukawa<n_center_type>;

//     auto gamma         = inputs.at("STG Exponent").value<element_type>();
//     auto kernel_output = submods.at("Yukawa kernel").value().run(inputs);
//     auto [X]           = kernel_type::unwrap_results(kernel_output);
//     const auto c0      = -element_type{1.0} / gamma;
//     tensor_type c0X;
//     auto idx = TA::detail::dummy_annotation(n_centers);
//     c0X(idx) = c0 * X(idx);

//     auto rv = results();
//     return PropType::wrap_results(rv, c0X);
// }

// template class STGGR<pt::gr2c<double>>;
// template class STGGR<pt::gr3c<double>>;
// template class STGGR<pt::gr4c<double>>;

// } // namespace integrals::f12
