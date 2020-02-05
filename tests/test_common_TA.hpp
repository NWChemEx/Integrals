#pragma once
#include <catch2/catch.hpp>
#include "integrals/nwx_libint/nwx_libint.hpp"
#include <tiledarray.h>
#include <libint2.hpp>
#include <fstream>

extern TA::World* pworld;

// Sample Basis Set
template<typename T = double>
libchemist::AOBasisSet<T> water_basis() {
    libchemist::AOBasisSet<T> bs;

    // Atomic basis set for the oxygen atom
    libchemist::Center<T> O(0.0, -0.1432223429807816, 0.0);
    O.add_shell(libchemist::ShellType::pure, 0,
                std::vector<T>{0.154329, 0.535328, 0.444635},
                std::vector<T>{130.709320, 23.808861, 6.443608});
    O.add_shell(libchemist::ShellType::pure, 0,
                std::vector<T>{-0.099967, 0.399513, 0.700115},
                std::vector<T>{5.033151, 1.169596, 0.380389});
    O.add_shell(libchemist::ShellType::pure, 1,
                std::vector<T>{0.155916, 0.607684, 0.391957},
                std::vector<T>{5.033151, 1.169596, 0.380389});

    // Atomic basis set for the first hydrogen atom
    libchemist::Center<T> H1(1.6380335020342418, 1.1365568803584036, 0.0);
    H1.add_shell(libchemist::ShellType::pure, 0,
                 std::vector<T>{0.15432897, 0.53532814, 0.44463454},
                 std::vector<T>{3.42525091, 0.62391373, 0.16885540});

    // Atomic basis set for the second hydrogen atom
    libchemist::Center<T> H2(-1.6380335020342418, 1.1365568803584036, 0.0);
    H2.add_shell(libchemist::ShellType::pure, 0,
                 std::vector<T>{0.15432897, 0.53532814, 0.44463454},
                 std::vector<T>{3.42525091, 0.62391373, 0.16885540});

    bs.add_center(O);
    bs.add_center(H1);
    bs.add_center(H2);

    return bs;
}

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
