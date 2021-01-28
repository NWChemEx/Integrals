// #include "transformed.hpp"

// namespace integrals {

// template<typename BaseType>
// TEMPLATED_MODULE_CTOR(Transformed, BaseType) {
//     using my_pt = pt::transformed<BaseType>;
//     satisfies_property_type<my_pt>();
//     add_submodule<BaseType>("AO integral");
// }

// template<typename BaseType>
// TEMPLATED_MODULE_RUN(Transformed, BaseType) {
//     using my_pt = pt::transformed<BaseType> :

//       constexpr auto n_centers =
//         property_types::ao_integrals::n_centers_v<BaseType>;

//     auto& ao_submod = submods.at("AO integral");
//     auto results    = ao_submod.run(inputs);
//     auto [X]        = BaseType::unwrap_results(results);

//     if constexpr(n_centers == 2) {
//         const auto& final_bra = inputs.at("final bra");
//         const auto& final_ket = inputs.at("final ket");
//         const auto& bra_space = inputs.at("bra");
//         const auto& ket_space = inputs.at("ket");

//         const bool transform_bra = final_bra != bra_space;
//         const bool transform_ket = final_ket != ket_space;
//     }
// }

// } // namespace integrals