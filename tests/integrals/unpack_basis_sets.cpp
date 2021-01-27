#include "../test_common_TA.hpp"
#include "integrals/unpack_basis_sets.hpp"
#include <catch2/catch.hpp>

using namespace integrals;

namespace {

auto generate_inputs() {
    sde::type::input_map inputs;
    auto sto_3g      = std::get<1>(make_molecule());
    auto cc_pvdz     = std::get<1>(make_molecule("cc-pvdz"));
    auto aug_cc_pvdz = std::get<1>(make_molecule("aug-cc-pvdz"));
    auto cc_pvtz     = std::get<1>(make_molecule("cc-pvtz"));
    using ao_space_t = decltype(sto_3g);
    inputs["bra"].set_type<ao_space_t>().change(sto_3g);
    inputs["bra 1"].set_type<ao_space_t>().change(sto_3g);
    inputs["bra 2"].set_type<ao_space_t>().change(cc_pvtz);
    inputs["ket"].set_type<ao_space_t>().change(cc_pvdz);
    inputs["ket 1"].set_type<ao_space_t>().change(cc_pvdz);
    inputs["ket 2"].set_type<ao_space_t>().change(aug_cc_pvdz);
    return inputs;
}

} // namespace

TEMPLATE_LIST_TEST_CASE("unpack_spaces", "", all_2c) {
    using space_type = integrals::type::ao_space_t<double>;
    auto inputs      = generate_inputs();
    auto rv          = unpack_spaces<TestType, space_type>(inputs);
    REQUIRE(rv[0] == inputs.at("bra").value<space_type>());
    REQUIRE(rv[1] == inputs.at("ket").value<space_type>());
}

TEMPLATE_LIST_TEST_CASE("unpack_spaces", "", all_3c) {
    using space_type = integrals::type::ao_space_t<double>;
    auto inputs      = generate_inputs();
    auto rv          = unpack_spaces<TestType, space_type>(inputs);
    REQUIRE(rv[0] == inputs.at("bra").value<space_type>());
    REQUIRE(rv[1] == inputs.at("ket 1").value<space_type>());
    REQUIRE(rv[2] == inputs.at("ket 2").value<space_type>());
}

TEMPLATE_LIST_TEST_CASE("unpack_spaces", "", all_4c) {
    using space_type = integrals::type::ao_space_t<double>;
    auto inputs      = generate_inputs();
    auto rv          = unpack_spaces<TestType, space_type>(inputs);
    if constexpr(property_types::ao_integrals::is_doi_v<TestType>) {
        REQUIRE(rv[0] == inputs.at("bra").value<space_type>());
        REQUIRE(rv[1] == inputs.at("bra").value<space_type>());
        REQUIRE(rv[2] == inputs.at("ket").value<space_type>());
        REQUIRE(rv[3] == inputs.at("ket").value<space_type>());
    } else {
        REQUIRE(rv[0] == inputs.at("bra 1").value<space_type>());
        REQUIRE(rv[1] == inputs.at("bra 2").value<space_type>());
        REQUIRE(rv[2] == inputs.at("ket 1").value<space_type>());
        REQUIRE(rv[3] == inputs.at("ket 2").value<space_type>());
    }
}

TEMPLATE_LIST_TEST_CASE("unpack_basis_sets", "", all_2c) {
    using space_type = integrals::type::ao_space_t<double>;
    auto inputs      = generate_inputs();
    auto rv          = unpack_basis_sets<TestType>(inputs);
    REQUIRE(rv[0] == inputs.at("bra").value<space_type>().basis_set());
    REQUIRE(rv[1] == inputs.at("ket").value<space_type>().basis_set());
}

TEMPLATE_LIST_TEST_CASE("unpack_basis_sets", "", all_3c) {
    using space_type = integrals::type::ao_space_t<double>;
    auto inputs      = generate_inputs();
    auto rv          = unpack_basis_sets<TestType>(inputs);
    REQUIRE(rv[0] == inputs.at("bra").value<space_type>().basis_set());
    REQUIRE(rv[1] == inputs.at("ket 1").value<space_type>().basis_set());
    REQUIRE(rv[2] == inputs.at("ket 2").value<space_type>().basis_set());
}

TEMPLATE_LIST_TEST_CASE("unpack_basis_sets", "", all_4c) {
    using space_type = integrals::type::ao_space_t<double>;
    auto inputs      = generate_inputs();
    auto rv          = unpack_basis_sets<TestType>(inputs);
    if constexpr(property_types::ao_integrals::is_doi_v<TestType>) {
        REQUIRE(rv[0] == inputs.at("bra").value<space_type>().basis_set());
        REQUIRE(rv[1] == inputs.at("bra").value<space_type>().basis_set());
        REQUIRE(rv[2] == inputs.at("ket").value<space_type>().basis_set());
        REQUIRE(rv[3] == inputs.at("ket").value<space_type>().basis_set());
    } else {
        REQUIRE(rv[0] == inputs.at("bra 1").value<space_type>().basis_set());
        REQUIRE(rv[1] == inputs.at("bra 2").value<space_type>().basis_set());
        REQUIRE(rv[2] == inputs.at("ket 1").value<space_type>().basis_set());
        REQUIRE(rv[3] == inputs.at("ket 2").value<space_type>().basis_set());
    }
}