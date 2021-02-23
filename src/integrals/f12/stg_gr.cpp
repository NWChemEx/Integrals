#include "f12.hpp"
#include "integrals/types.hpp"
#include <property_types/ao_integrals/type_traits.hpp>

using namespace property_types::ao_integrals;

static constexpr auto module_desc = R"""(
    TODO: Write me!!!!
)""";

namespace integrals::f12 {

template<typename PropType>
TEMPLATED_MODULE_CTOR(STGGR, PropType) {
    using element_type       = element_t<PropType>;
    constexpr auto n_centers = n_centers_v<PropType>;
    using kernel_type        = Yukawa<NCenter<n_centers, element_type>>;

    satisfies_property_type<PropType>();

    add_submodule<kernel_type>("Yukawa kernel");
    add_input<element_type>("STG Exponent").set_default(1.2);
}

template<typename PropType>
TEMPLATED_MODULE_RUN(STGGR, PropType) {
    using element_type       = element_t<PropType>;
    using tensor_type        = type::tensor<element_type>;
    constexpr auto n_centers = n_centers_v<PropType>;
    using kernel_type        = Yukawa<NCenter<n_centers, element_type>>;

    auto gamma         = inputs.at("STG Exponent").value<element_type>();
    auto kernel_output = submods.at("Yukawa kernel").value().run(inputs);
    auto [X]           = kernel_type::unwrap_results(kernel_output);
    const auto c0      = -element_type{1.0} / gamma;
    tensor_type c0X;
    auto idx = TA::detail::dummy_annotation(n_centers);
    c0X(idx) = c0 * X(idx);

    auto rv = results();
    return PropType::wrap_results(rv, c0X);
}

template class STGGR<pt::gr2c<double>>;
template class STGGR<pt::gr3c<double>>;
template class STGGR<pt::gr4c<double>>;

} // namespace integrals::f12