#pragma once
#include <catch2/catch.hpp>
#include <tiledarray.h>
#include <libint2.hpp>
#include <libchemist/libchemist.hpp>
#include <integrals/types.hpp>

using TensorType = integrals::type::tensor<double>;
using IndexType = std::vector<std::size_t>;
using BlockTensor = std::map<IndexType, std::vector<double>>;
using string_list = std::initializer_list<std::string>;

inline auto make_molecule(const std::string& bs_name = "sto-3g") {
    using libchemist::Atom;
    using c_t = typename Atom::coord_type;
    Atom O{8ul, c_t{0.000000000000000, -0.143222342980786, 0.000000000000000}};
    Atom H1{1ul, c_t{1.638033502034240, 1.136556880358410, 0.000000000000000}};
    Atom H2{1ul, c_t{-1.638033502034240, 1.136556880358410, 0.000000000000000}};
    libchemist::Molecule water(O, H1, H2);
    return std::make_tuple(water, libchemist::apply_basis(bs_name, water));
}


inline void compare_integrals(TensorType& calc, BlockTensor& corr,
                              const double eps  = 10000 * std::numeric_limits<double>::epsilon(),
                              const double marg = 100 * std::numeric_limits<double>::epsilon()) {
    REQUIRE(calc.size() == corr.size());

    for (const auto& elem : corr) {
        auto idx = elem.first;
        auto corr_block = elem.second;
        auto is_zero = calc.shape().is_zero(idx);

        if (is_zero) {
            for (double val : corr_block) {
                REQUIRE(0.0 == Approx(val).epsilon(eps).margin(marg));
            }
        } else {
            auto calc_block = calc.find(idx).get();

            REQUIRE(calc_block.size() <= corr_block.size());
            for (int i = 0; i < calc_block.size(); ++i) {
                REQUIRE(calc_block[i] == Approx(corr_block[i]).epsilon(eps).margin(marg));
            }
        }
    }
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
    SECTION("inputs"){
        auto inputs = T::inputs();
        REQUIRE(inputs.size() == input_fields.size());
        for(const auto& field : input_fields)
            SECTION(field){ REQUIRE(inputs.count(field) == 1); }

    }
    SECTION("results"){
        auto results = T::results();
        REQUIRE(results.size() == result_fields.size());
        for(const auto& field : result_fields)
            SECTION(field){ REQUIRE(results.count(field) == 1); }
    }
}