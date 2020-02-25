#pragma once
#include <libint2.hpp>
#include <tiledarray.h>
#include "integrals/nwx_libint/nwx_libint_factory.hpp"
#include "integrals/nwx_TA/nwx_TA_utils.hpp"

namespace nwx_TA {

    template<typename val_type, libint2::Operator op, std::size_t NBases>
    struct FillNDFunctor {

        using basis_vec = std::vector<libint2::BasisSet>;
        using size_vec = std::vector<std::size_t>;

        // The collected LibInt2 basis sets needed for the integral
        basis_vec LIBasis_sets;

        // The factory that produces the appropriate LibInt2 engines
        nwx_libint::LibintFactory<NBases, op> factory;

        // Number of arrays returned by operator
        std::size_t nopers = libint2::operator_traits<op>::nopers;

        // Initialize and finalize LibInt2
        FillNDFunctor() { libint2::initialize(); }
        ~FillNDFunctor() { libint2::finalize(); }

        // Complies with the TA API for these functions
        float operator()(val_type& tile, const TiledArray::Range& range) {
            return _fill(tile, range);
        }

    private:

        /** @brief The top level function that starts the recursive calls of the other functions.
         *         Gets the tile @p tile and range @p range from TA, then initializes the tile and
         *         declares the vectors for the offsets and shells.
         *
         *  @param[in] tile The tile to be filled
         *  @param[in] range The range of the tile
         *  @returns The norm of the filled tile
         */
        float _fill(val_type& tile, const TiledArray::Range& range) {
            tile = val_type(range);

            // In case something else finalized
            if (not libint2::initialized()) { libint2::initialize(); }

            // Make libint engine
            auto tile_engine = factory();

            // Vector for storing the offsets of the current shells
            size_vec offsets(NBases);

            // Vector for storing the indices of the current shells
            size_vec shells(NBases);

            // Start recurvise process to determine the shells needed for the current tile
            _index_shells(tile, range, tile_engine, offsets, shells, 0);

            // Return norm for new tile
            return tile.norm();
        }

        /** @brief Recursive function that transverses all of the dimensions of the current tile
         *         and finds the shells that need to be computed to fill the tile.
         *
         *  @param tile The tile to be filled
         *  @param range The range of the tile
         *  @param tile_engine The LibInt2 engine that computes integrals
         *  @param offsets Vector containing the coordinate offset based on the current AOs
         *  @param shells Vector tracking the indices of the current shells
         *  @param depth The current dimension of the desired tensor
         */
        void _index_shells(val_type& tile,
                           const TiledArray::Range& range,
                           libint2::Engine& tile_engine,
                           size_vec& offsets,
                           size_vec& shells,
                           int depth) {
            // Deal with the additional dimension of the tile for multipoles
            int tile_depth = (nopers == 1) ? depth : depth + 1;

            // Set initial offset as the dimensions lower bound
            offsets[depth] = range.lobound()[tile_depth];

            // Loop over the shells that contain the AOs spanned by the tile in this dimension
            for (auto s : aos2shells(LIBasis_sets[depth], range.lobound()[tile_depth], range.upbound()[tile_depth])) {
                // Save index of current shell to be passed down the line
                shells[depth] = s;

                if (depth == (NBases - 1)) {
                    // Compute integrals for current shells
                    _call_libint(tile_engine, shells, std::make_index_sequence<NBases>());

                    // Switch for multipole cases
                    if (nopers == 1) {
                        // Initially set indexer to first coordinate of tile
                        size_vec indexer = offsets;

                        // Store the location of the results
                        const auto& int_vals = tile_engine.results()[0];

                        // Index for integrals values; referenced so it gets pasted around
                        int int_i = 0;

                        // Fill in the values of the tile
                        _fill_from_libint(tile, offsets, shells, int_vals, indexer, int_i, 0);
                    } else {
                        // Deal with the additional dimension of the tile for multipoles
                        size_vec indexer = {range.lobound()[0]};
                        indexer.insert(indexer.end(), offsets.begin(), offsets.end());

                        // Loop over multipoles to fill in values
                        for (auto f = 0ul; f != nopers; ++f) {
                            // Store the location of the results
                            const auto& int_vals = tile_engine.results()[f];

                            // Index for integrals values; referenced so it gets pasted around
                            int int_i = 0;

                            // Fill in the values of the multipole
                            _fill_from_libint(tile, offsets, shells, int_vals, indexer, int_i, 0);

                            // Increment the coordinate for the multipoles
                            indexer[0]++;
                        }
                    }
                } else {
                    // Keep going down until all dimensions of the the tensor have been covered
                    _index_shells(tile, range, tile_engine, offsets, shells, depth + 1);
                }
                // Increase the offsets for the current dimension by the size of the completed shell
                offsets[depth] += LIBasis_sets[depth][shells[depth]].size();
            }
        }

        /** @brief Recursive function that transverses all of the dimensions of the current tile
         *         and fills the LibInt2 results into the correct coordinate position in the tile
         *
         *  @param tile The tile to be filled
         *  @param offsets Vector containing the coordinate offset based on the current AOs
         *  @param shells Vector tracking the indices of the current shells
         *  @param int_vals Pointer to the beginning of the integral values
         *  @param indexer Index for placing values into array
         *  @param int_i Index for accessing the current integral value to be filled
         *  @param depth The current dimension of the desired tensor
         */
        void _fill_from_libint(val_type& tile,
                               const size_vec& offsets,
                               const size_vec& shells,
                               const double* int_vals,
                               size_vec& indexer,
                               int& int_i,
                               int depth) {
            // Deal with the additional dimension of the tile for multipoles
            int tile_depth = (nopers == 1) ? depth : depth + 1;

            // Loop over the size of the current shell for this dimension
            for (auto f = 0ul; f != LIBasis_sets[depth][shells[depth]].size(); ++f) {
                if (depth == (NBases - 1)) {
                    // Assign integral value
                    if (int_vals == nullptr) {
                        // Default case of all zeroes
                        tile[indexer] = 0;
                    } else {
                        // Get values from LibInt
                        tile[indexer] = int_vals[int_i];
                        int_i++;
                    }
                } else {
                    // Keep going down until all dimensions of the the tensor have been covered
                    _fill_from_libint(tile, offsets, shells, int_vals, indexer, int_i, depth + 1);
                }
                // Increment the coordinate for this dimension
                indexer[tile_depth]++;
            }
            // Reset the indexer to the initial position for next loop
            indexer[tile_depth] = offsets[depth];
        }

        /** @brief Wrap the call of LibInt2 engine so it can take a variable number of shell inputs.
         *
         * @tparam Is A variadic parameter pack of integers from [0,NBases) to expand.
         * @param tile_engine The LibInt2 engine that computes integrals
         * @param shells The index of the requested shell block
         * @return An std::vector filled with the requested block per operator component
         */
        template<std::size_t... Is>
        void _call_libint(libint2::Engine& tile_engine, size_vec shells, std::index_sequence<Is ...>) {
            tile_engine.compute(LIBasis_sets[Is][shells[Is]] ...);
        }
    };

} // namespace nwx_TA