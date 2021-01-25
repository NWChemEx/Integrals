#pragma once
#include <catch2/catch.hpp>
#include <integrals/types.hpp>
#include <libchemist/libchemist.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <libint2.hpp>
#include <tiledarray.h>

using TensorType  = integrals::type::tensor<double>;
using IndexType   = std::vector<std::size_t>;
using string_list = std::initializer_list<std::string>;

inline auto make_molecule(const std::string& bs_name = "sto-3g") {
    using libchemist::Atom;
    using c_t = typename Atom::coord_type;
    Atom O{8ul, c_t{0.000000000000000, -0.143222342980786, 0.000000000000000}};
    Atom H1{1ul, c_t{1.638033502034240, 1.136556880358410, 0.000000000000000}};
    Atom H2{1ul, c_t{-1.638033502034240, 1.136556880358410, 0.000000000000000}};
    libchemist::Molecule water(O, H1, H2);
    auto bs = libchemist::apply_basis(bs_name, water);
    return std::make_tuple(water,
                           libchemist::orbital_space::AOSpace<double>(bs));
}

/**@brief Factors out the boilerplate required to test a property type
 *
 * @tparam T The type of the property type
 * @param input_fields An initializer list of the property type's input fields
 * @param result_fields An initializer list of the property type's returns
 */
template<typename T>
inline static void test_property_type(string_list input_fields,
                                      string_list result_fields) {
    SECTION("inputs") {
        auto inputs = T::inputs();
        REQUIRE(inputs.size() == input_fields.size());
        for(const auto& field : input_fields) SECTION(field) {
                REQUIRE(inputs.count(field) == 1);
            }
    }
    SECTION("results") {
        auto results = T::results();
        REQUIRE(results.size() == result_fields.size());
        for(const auto& field : result_fields) SECTION(field) {
                REQUIRE(results.count(field) == 1);
            }
    }
}