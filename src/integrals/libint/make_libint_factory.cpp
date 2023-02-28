// #include "detail_/make_libint_basis_set.hpp"
// #include "libint_factory.hpp"
// #include "libint_op.hpp"
// #include <simde/integral_factory.hpp>

// namespace integrals::libint {

// template<typename std::size_t N, typename OperatorType>
// TEMPLATED_MODULE_CTOR(MakeLibintFactory, N, OperatorType) {
//     using my_pt = simde::IntegralFactory<OperatorType>;

//     satisfies_property_type<my_pt>();

//     add_input<double>("Threshold")
//       .set_default(1.0E-16)
//       .set_description(
//         "The target precision with which the integrals will be computed");
// }

// template<typename std::size_t N, typename OperatorType>
// TEMPLATED_MODULE_RUN(MakeLibintFactory, N, OperatorType) {
//     using my_pt               = simde::IntegralFactory<OperatorType>;
//     using libint_factory      = LibintFactory<N>;
//     using libint_basis_vector = typename libint_factory::libint_basis_vector;

//     const auto& [bases, op] = my_pt::unwrap_inputs(inputs);
//     auto thresh             = inputs.at("Threshold").value<double>();

//     // Convert from NWX bases to libint
//     libint_basis_vector libint_bases;
//     for(const auto& basis_i : bases)
//         libint_bases.push_back(detail_::make_libint_basis_set(basis_i));

//     auto deriv = 0; // TODO: Get from TMP

//     constexpr auto libint_op = op_v<OperatorType>;

//     auto pfactory = std::make_unique<libint_factory>(std::move(libint_bases),
//                                                      libint_op, thresh, deriv);

//     IntegralFactory fac(std::move(pfactory));
//     auto rv = results();
//     return my_pt
// }

// } // namespace integrals::libint
