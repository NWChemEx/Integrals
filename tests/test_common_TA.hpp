#pragma once
#include <catch2/catch.hpp>
#include <tiledarray.h>
#include <libint2.hpp>
#include <libchemist/libchemist.hpp>
#include <integrals/types.hpp>

using TensorType = integrals::type::tensor<double>;
using IndexType = std::vector<std::size_t>;
using BlockTensor = std::map<IndexType, std::vector<double>>;

inline auto make_molecule() {
    using libchemist::Atom;
    using c_t = typename Atom::coord_type;
    Atom O{8ul, c_t{0.000000000000000, -0.143222342980786, 0.000000000000000}};
    Atom H1{1ul, c_t{1.638033502034240, 1.136556880358410, 0.000000000000000}};
    Atom H2{1ul, c_t{-1.638033502034240, 1.136556880358410, 0.000000000000000}};
    libchemist::Molecule water(O, H1, H2);
    return std::make_tuple(water, libchemist::apply_basis("sto-3g", water));
}


inline void compare_integrals(TensorType& calc, BlockTensor& corr,
                              const double eps  = 10000 * std::numeric_limits<double>::epsilon(),
                              const double marg = 100 * std::numeric_limits<double>::epsilon()) {
    REQUIRE(calc.size() == corr.size());

    for (const auto& it : calc) {
        auto calc_block = it.get();
        auto& corr_block = corr.at(it.index());

        REQUIRE(calc_block.size() <= corr_block.size());
        auto check_len = std::min(calc_block.size(), corr_block.size());

        for (int i = 0; i < check_len; ++i) {
            REQUIRE(calc_block[i] == Approx(corr_block[i]).epsilon(eps).margin(marg));
        }
    }
}