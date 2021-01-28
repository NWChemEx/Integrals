#pragma once
#include <catch2/catch.hpp>
#include <integrals/types.hpp>
#include <libchemist/libchemist.hpp>
#include <libchemist/ta_helpers/ta_helpers.hpp>
#include <libint2.hpp>
#include <property_types/ao_integrals/ao_integrals.hpp>
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

// The following types will be used to loop over template parameters in some
// unit tests

// Tuple of the 2-centered integral property types templated on element type
template<typename T>
using two_center = std::tuple<property_types::ao_integrals::ERI2C<T>,
                              property_types::ao_integrals::EDipole<T>,
                              property_types::ao_integrals::EQuadrupole<T>,
                              property_types::ao_integrals::EOctopole<T>,
                              property_types::ao_integrals::Kinetic<T>,
                              property_types::ao_integrals::Nuclear<T>,
                              property_types::ao_integrals::Overlap<T>,
                              property_types::ao_integrals::STG2C<T>,
                              property_types::ao_integrals::Yukawa2C<T>>;

// Tuple of the 3-centered integral property types templated on element type
template<typename T>
using three_center = std::tuple<property_types::ao_integrals::ERI3C<T>,
                                property_types::ao_integrals::STG3C<T>,
                                property_types::ao_integrals::Yukawa3C<T>>;

// Tuple of the 4-centered integral property types templated on element type
template<typename T>
using four_center = std::tuple<property_types::ao_integrals::DOI<T>,
                               property_types::ao_integrals::ERI4C<T>,
                               property_types::ao_integrals::STG4C<T>,
                               property_types::ao_integrals::Yukawa4C<T>>;

// Tuple with all 2-centered integral types in it
using all_2c =
  decltype(std::tuple_cat(two_center<float>{}, two_center<double>{}));

using all_3c =
  decltype(std::tuple_cat(three_center<float>{}, three_center<double>{}));

using all_4c =
  decltype(std::tuple_cat(four_center<float>{}, four_center<double>{}));
