#include "../test_common_TA.hpp"
#include "integrals/unpack_basis_sets.hpp"

namespace {
template<typename T>
auto generate_inputs() {
    sde::type::input_map inputs;
    using bs_t       = libchemist::AOBasisSet<T>;
    using ao_space_t = integrals::type::ao_space_t<T>;
    bs_t bs0, bs1, bs2, bs3;
    libchemist::Center<T> c0{0.0, 0.0, 1.0}, c1{0.0, 1.0, 0.0};
    bs1.add_center(c0);
    bs2.add_center(c1);
    bs3.add_center(c0);
    bs3.add_center(c1);
    inputs["bra"].set_type<ao_space_t>().change(ao_space_t{bs0});
    inputs["bra 1"].set_type<ao_space_t>().change(ao_space_t{bs0});
    inputs["bra 2"].set_type<ao_space_t>().change(ao_space_t{bs1});
    inputs["ket"].set_type<ao_space_t>().change(ao_space_t{bs2});
    inputs["ket 1"].set_type<ao_space_t>().change(ao_space_t{bs2});
    inputs["ket 2"].set_type<ao_space_t>().change(ao_space_t{bs3});
    return inputs;
}

} // namespace

TEMPLATE_LIST_TEST_CASE("unpack_basis_sets", "", all_2c) {
    using element_type = double;
    using space_type   = property_types::type::ao_space_t<element_type>;
    auto inputs        = generate_inputs<element_type>();
    auto rv            = integrals::unpack_basis_sets<TestType>(inputs);
    REQUIRE(rv.size() == 2);
    REQUIRE(rv[0] == inputs.at("bra").value<space_type>().basis_set());
    REQUIRE(rv[1] == inputs.at("ket").value<space_type>().basis_set());
}

TEMPLATE_LIST_TEST_CASE("unpack_basis_sets", "", all_3c) {
    using element_type = double;
    using space_type   = property_types::type::ao_space_t<element_type>;
    auto inputs        = generate_inputs<element_type>();
    auto rv            = integrals::unpack_basis_sets<TestType>(inputs);
    REQUIRE(rv.size() == 3);
    REQUIRE(rv[0] == inputs.at("bra").value<space_type>().basis_set());
    REQUIRE(rv[1] == inputs.at("ket 1").value<space_type>().basis_set());
    REQUIRE(rv[2] == inputs.at("ket 2").value<space_type>().basis_set());
}

TEMPLATE_LIST_TEST_CASE("unpack_basis_sets", "", all_4c) {
    using element_type = double;
    using space_type   = property_types::type::ao_space_t<element_type>;
    auto inputs        = generate_inputs<element_type>();
    auto rv            = integrals::unpack_basis_sets<TestType>(inputs);
    if constexpr(property_types::ao_integrals::is_doi_v<TestType>) {
        REQUIRE(rv.size() == 2);
        REQUIRE(rv[0] == inputs.at("bra").value<space_type>().basis_set());
        REQUIRE(rv[1] == inputs.at("ket").value<space_type>().basis_set());
    } else {
        REQUIRE(rv.size() == 4);
        REQUIRE(rv[0] == inputs.at("bra 1").value<space_type>().basis_set());
        REQUIRE(rv[1] == inputs.at("bra 2").value<space_type>().basis_set());
        REQUIRE(rv[2] == inputs.at("ket 1").value<space_type>().basis_set());
        REQUIRE(rv[3] == inputs.at("ket 2").value<space_type>().basis_set());
    }
}