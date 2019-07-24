#pragma once
#include <catch/catch.hpp>
#include <iomanip>
#include <iostream>
#include <libchemist/libchemist.hpp>
#include <map>
#include <tamm/tamm.hpp>

using IndexType   = std::vector<unsigned int>;
using BlockTensor = std::map<IndexType, std::vector<double>>;
using TensorType  = tamm::Tensor<double>;

inline void compare_impl(const TensorType& calc, const BlockTensor& corr,
                         IndexType idx,
                         const double eps,
                         const double marg,
                         size_t depth = 0) {
    auto idx_spaces   = calc.tiled_index_spaces();
    const auto nmodes = idx_spaces.size();

    if(depth == nmodes) { // End recursion
        // Get buffer size
        size_t block_size = 1;
        for(size_t i = 0; i < nmodes; ++i)
            block_size *= idx_spaces[i].tile_size(idx[i]);

        std::vector<double> buffer(block_size);
        // Avoids warning for size_t to long int
        long int dim = buffer.size();
        calc.get(idx, {buffer.data(), dim});

        for(auto x = 0; x < block_size; ++x) {
            const auto abs_error = std::abs(buffer[x] - corr.at(idx)[x]);
            if (abs_error > marg && abs_error/std::abs(buffer[x]) > eps) {
              std::cerr << std::setprecision(20) << std::scientific << "value=" << buffer[x] << " ref_value=" << corr.at(idx)[x] << " abs_error=" << abs_error << " rel_error=" << (abs_error/std::abs(buffer[x])) << std::endl;
            }
            REQUIRE(buffer[x] ==
                    Approx(corr.at(idx)[x]).epsilon(eps).margin(marg));
        }
    } else { // Adjust the depth-th mode's index
        const auto ntiles = idx_spaces[depth].num_tiles();
        for(size_t blocki = 0; blocki < ntiles; ++blocki) {
            idx[depth] = blocki;
            // Reset indices after depth
            for(size_t dimi = depth + 1; dimi < nmodes; ++dimi) idx[dimi] = 0;
            compare_impl(calc, corr, idx, eps, marg, depth + 1);
        }
    }
}

/// compares @c values to @c ref_values
/// @param[in] values the values to be checked
/// @param[in] ref_values the reference values against which @c values will be tested
/// @param[in] epsilon the relative error margin
/// @param[in] margin the absolute error margin
/// @note two values are deemed equal if either they differ by less than @c margin or
///       their relative difference is less than @c epsilon. To check absolute errors only
///       set @c epsilon to zero. To check relative errors only set @c margin to zero.
inline void compare_integrals(const TensorType& values, const BlockTensor& ref_values,
                              const double epsilon  = 10000 * std::numeric_limits<double>::epsilon(),
                              const double margin = 100 * std::numeric_limits<double>::epsilon()) {
    compare_impl(calc, corr, IndexType(calc.tiled_index_spaces().size()), epsilon, margin);
}

// prints the integrals in the format used for correctness checks
inline void print_impl(const TensorType& calc, IndexType idx,
                       size_t depth = 0) {
    auto idx_spaces   = calc.tiled_index_spaces();
    const auto nmodes = idx_spaces.size();

    if(depth == nmodes) { // End recursion
        // Get buffer size
        size_t block_size = 1;
        for(size_t i = 0; i < nmodes; ++i)
            block_size *= idx_spaces[i].tile_size(idx[i]);

        std::vector<double> buffer(block_size);
        // Avoids warning for size_t to long int
        long int dim = buffer.size();
        calc.get(idx, {buffer.data(), dim});

        std::cout << "{{";
        for(auto i : idx) std::cout << i << ", ";
        std::cout << "}, {";
        for(auto x = 0; x < block_size; ++x) {
            if(x % 5 == 0) std::cout << std::endl;
            std::cout << std::fixed << std::setprecision(16);
            std::cout << buffer[x] << ",";
        }
        std::cout << "}}," << std::endl;
    } else { // Adjust the depth-th mode's index
        const auto ntiles = idx_spaces[depth].num_tiles();
        for(size_t blocki = 0; blocki < ntiles; ++blocki) {
            idx[depth] = blocki;
            // Reset indices after depth
            for(size_t dimi = depth + 1; dimi < nmodes; ++dimi) idx[dimi] = 0;
            print_impl(calc, idx, depth + 1);
        }
    }
}

inline void print_integrals(const TensorType& calc) {
    print_impl(calc, IndexType(calc.tiled_index_spaces().size()));
}

inline auto make_molecule() {
    using libchemist::Atom;
    using c_t = typename Atom::coord_type;
    Atom H1{1ul, c_t{1.638033502034240, 1.136556880358410, 0.000000000000000}};
    Atom O{8ul, c_t{0.000000000000000, -0.143222342980786, 0.000000000000000}};
    Atom H2{1ul, c_t{-1.638033502034240, 1.136556880358410, 0.000000000000000}};
    libchemist::Molecule water(O, H1, H2);
    return std::make_tuple(water, libchemist::apply_basis("sto-3g", water));
}
