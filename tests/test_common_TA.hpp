#pragma once
#include <catch2/catch.hpp>
#include <tiledarray.h>
#include <libint2.hpp>
#include <libchemist/libchemist.hpp>

extern TA::World* pworld;

// Make tile range based on BS
template<typename T = double>
TA::TiledRange1 make_tiled_range(libchemist::AOBasisSet<T> basis_set, std::size_t tile_size) {
    std::vector<std::size_t> bounds{0};

    auto span = 0;
    for (auto ishell = 0; ishell < basis_set.n_shells(); ++ishell) {
        span += basis_set.shell(ishell).size();

        if (span < tile_size) continue;

        auto bound = bounds.back() + span;
        bounds.push_back(bound);
        span = 0;
    }
    if (span != 0) bounds.push_back(bounds.back() + span);

    return TA::TiledRange1(bounds.begin(), bounds.end());
}

// Given AO bounds, return shells containing those AOS
std::vector<std::size_t> aos2shells(libint2::BasisSet basis_set, std::size_t lower, std::size_t upper) {
    std::vector<std::size_t> return_vec;

    for (auto ishell = 0, offset = 0; ishell < basis_set.size(); ++ishell) {
        if (offset >= upper) break;
        if (offset >= lower) return_vec.push_back(ishell);
        offset += basis_set[ishell].size();
    }

    return return_vec;
}

inline auto make_molecule() {
    using libchemist::Atom;
    using c_t = typename Atom::coord_type;
    Atom O{8ul, c_t{0.000000000000000, -0.143222342980786, 0.000000000000000}};
    Atom H1{1ul, c_t{1.638033502034240, 1.136556880358410, 0.000000000000000}};
    Atom H2{1ul, c_t{-1.638033502034240, 1.136556880358410, 0.000000000000000}};
    libchemist::Molecule water(O, H1, H2);
    return std::make_tuple(water, libchemist::apply_basis("sto-3g", water));
}
