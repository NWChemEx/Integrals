/*
 * Copyright 2024 NWChemEx-Project
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
#include "detail_/make_libint_basis_set.hpp"
#include <type_traits>

namespace integrals::ao_integrals {

template<typename BraType, typename KetType>
constexpr int get_n(const BraType& bra, const KetType& ket) {
    constexpr auto bra_is_aos         = std::is_same_v<BraType, aos>;
    constexpr auto bra_is_aos_squared = std::is_same_v<BraType, aos_squared>;
    constexpr auto ket_is_aos         = std::is_same_v<KetType, aos>;
    constexpr auto ket_is_aos_squared = std::is_same_v<KetType, aos_squared>;
    if constexpr(bra_is_aos && ket_is_aos) {
        return 2;
    } else if constexpr(bra_is_aos && ket_is_aos_squared) {
        return 3;
    } else if constexpr(bra_is_aos_squared && ket_is_aos_squared) {
        return 4;
    }
}

template<typename BraType, typename KetType>
std::vector<libint2::BasisSet> get_basis_sets(const BraType& bra,
                                              const KetType& ket) {
    using detail_::make_libint_basis_set;

    std::vector<libint2::BasisSet> basis_sets;

    if constexpr(std::is_same_v<BraType, aos>) {
        basis_sets.push_back(make_libint_basis_set(bra.ao_basis_set()));
    } else if constexpr(std::is_same_v<BraType, aos_squared>) {
        basis_sets.push_back(make_libint_basis_set(bra.lhs().ao_basis_set()));
        basis_sets.push_back(make_libint_basis_set(bra.rhs().ao_basis_set()));
    }

    if constexpr(std::is_same_v<KetType, aos>) {
        basis_sets.push_back(make_libint_basis_set(ket.ao_basis_set()));
    } else if constexpr(std::is_same_v<KetType, aos_squared>) {
        basis_sets.push_back(make_libint_basis_set(ket.lhs().ao_basis_set()));
        basis_sets.push_back(make_libint_basis_set(ket.rhs().ao_basis_set()));
    }

    return basis_sets;
}

template<typename BraKetType>
TEMPLATED_MODULE_CTOR(AOIntegral, BraKetType) {
    using my_pt = simde::EvaluateBraKet<BraKetType>;
    satisfies_property_type<my_pt>();
    description("Computes integrals with Libint");
}

template<typename BraKetType>
TEMPLATED_MODULE_RUN(AOIntegral, BraKetType) {
    using my_pt = simde::EvaluateBraKet<BraKetType>;

    const auto& [braket] = my_pt::unwrap_inputs(inputs);
    auto bra             = braket.bra();
    auto ket             = braket.ket();
    auto op              = braket.op();

    auto basis_sets = get_basis_sets(bra, ket);
    constexpr int N = get_n(bra, ket);

    using tensor_t = simde::type::tensor;
    using shape_t  = tensorwrapper::shape::Smooth;
    using layout_t = tensorwrapper::layout::Physical;
    using buffer_t = tensorwrapper::buffer::Eigen<double, N>;
    using data_t   = typename buffer_t::data_type;

    shape_t s{3, 3};
    layout_t l(s);
    data_t m(3, 3);
    buffer_t b{m, l};
    for(std::size_t i = 0; i < 3; ++i) {
        for(std::size_t j = 0; j < 3; ++j) {
            b.value()(i, j) = (i + 1) * (j + 1);
        }
    }
    tensor_t t({s, b});

    auto rv = results();
    return my_pt::wrap_results(rv, t);
}

#define AOI(bra, op, ket) AOIntegral<braket<bra, op, ket>>
#define EXTERN_AOI(bra, op, ket) template struct AOI(bra, op, ket)
#define LOAD_AOI(bra, op, ket, key) mm.add_module<AOI(bra, op, ket)>(key)

EXTERN_AOI(aos, t_e_type, aos);
EXTERN_AOI(aos, v_en_type, aos);
EXTERN_AOI(aos, v_ee_type, aos);
EXTERN_AOI(aos, v_ee_type, aos_squared);
EXTERN_AOI(aos_squared, v_ee_type, aos_squared);

void ao_integrals_set_defaults(pluginplay::ModuleManager& mm) {
    // Set any default associations
}

void load_ao_integrals(pluginplay::ModuleManager& mm) {
    LOAD_AOI(aos, t_e_type, aos, "Kinetic");
    LOAD_AOI(aos, v_en_type, aos, "Nuclear");
    LOAD_AOI(aos, v_ee_type, aos, "ERI2");
    LOAD_AOI(aos, v_ee_type, aos_squared, "ERI3");
    LOAD_AOI(aos_squared, v_ee_type, aos_squared, "ERI4");
    ao_integrals_set_defaults(mm);
}

#undef AOI
#undef ADD_AOI

} // namespace integrals::ao_integrals