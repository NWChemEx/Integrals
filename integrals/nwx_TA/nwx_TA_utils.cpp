#include "integrals/nwx_TA/nwx_TA_utils.hpp"

#include <utility>

namespace nwx_TA {

    template<typename T>
    TA::TiledRange1 _make_tiled_range(basis<T> basis_set, size tile_size) {
        std::vector<size> bounds{0};

        auto span = 0;
        for (const auto& atom : basis_set) {
            span += atom.n_aos();

            if (span < tile_size) continue;

            bounds.push_back(bounds.back() + span);
            span = 0;
        }
        if (span != 0) {
            if ((span < (tile_size / 5)) && (bounds.size() != 1)) {
                bounds.back() += span;
            } else {
                bounds.push_back(bounds.back() + span);
            }
        }

        return TA::TiledRange1(bounds.begin(), bounds.end());
    }

    TA::TiledRange1 make_tiled_range(basis<double> basis_set, size tile_size) {
        return _make_tiled_range<double>(std::move(basis_set), tile_size);
    }

    TA::TiledRange1 make_tiled_range(basis<float> basis_set, size tile_size) {
        return _make_tiled_range<float>(std::move(basis_set), tile_size);
    }

    template<typename T>
    TA::TiledRange1 _make_tiled_range(basis<T> basis_set, std::vector<size> tile_sizes) {
        if (tile_sizes.size() == 1) return make_tiled_range(basis_set, tile_sizes[0]);

        std::vector<size> bounds{0};
        auto sizes_index = 0;

        auto span = 0;
        for (const auto& atom : basis_set) {
            span += atom.n_aos();

            if (span < tile_sizes[sizes_index]) continue;

            bounds.push_back(bounds.back() + span);
            span = 0;

            sizes_index++;
            if (sizes_index >= tile_sizes.size()) {
                sizes_index = 0;
            }
        }
        if (span != 0) bounds.push_back(bounds.back() + span);

        return TA::TiledRange1(bounds.begin(), bounds.end());
    }

    TA::TiledRange1 make_tiled_range(basis<double> basis_set, std::vector<size> tile_sizes) {
        return _make_tiled_range<double>(std::move(basis_set), std::move(tile_sizes));
    }

    TA::TiledRange1 make_tiled_range(basis<float> basis_set, std::vector<size> tile_sizes) {
        return _make_tiled_range<float>(std::move(basis_set), std::move(tile_sizes));
    }

    TA::TiledRange1 make_tiled_range(size upper, size tile_size) {
        std::vector<size> bounds{};

        for (auto bound = 0; bound < upper; bound+=tile_size) bounds.push_back(bound);
        if (bounds.back() != upper) bounds.push_back(upper);

        return TA::TiledRange1(bounds.begin(), bounds.end());
    }

    TA::TiledRange1 make_tiled_range(size upper, std::vector<size> tile_sizes) {
        if (tile_sizes.size() == 1) return make_tiled_range(upper, tile_sizes[0]);

        std::vector<size> bounds{};
        auto sizes_index = 0;
        auto bound = 0;

        while (bound < upper) {
            bounds.push_back(bound);

            bound += tile_sizes[sizes_index];

            sizes_index++;
            if (sizes_index >= tile_sizes.size()) {
                sizes_index = 0;
            }
        }
        if (bounds.back() != upper) bounds.push_back(upper);

        return TA::TiledRange1(bounds.begin(), bounds.end());
    }

    template<typename T>
    std::vector<size> _aos2shells(basis<T> basis_set, size lower, size upper) {
        std::vector<size> return_vec;

        for (auto ishell = 0, offset = 0; ishell < basis_set.n_shells(); ++ishell) {
            if (offset >= upper) break;
            if (offset >= lower) return_vec.push_back(ishell);
            offset += basis_set.shell(ishell).size();
        }

        return return_vec;
    }

    std::vector<size> aos2shells(basis<double> basis_set, size lower, size upper) {
        return _aos2shells(std::move(basis_set), lower, upper);
    }

    std::vector<size> aos2shells(basis<float> basis_set, size lower, size upper) {
        return _aos2shells(std::move(basis_set), lower, upper);
    }

    template<typename basis_type>
    TA::TiledRange _make_trange(const std::vector<basis_type>& basis_sets,
                               const std::vector<size>& tile_sizes,
                               std::vector<TA::TiledRange1> ranges) {
        for (const auto& basis_set : basis_sets) ranges.push_back(make_tiled_range(basis_set, tile_sizes));
        TA::TiledRange trange(ranges.begin(), ranges.end());
        return trange;
    }

    TA::TiledRange make_trange(const std::vector<basis<double>>& basis_sets,
                               const std::vector<size>& tile_sizes,
                               std::vector<TA::TiledRange1> ranges) {
        return _make_trange(basis_sets, tile_sizes, std::move(ranges));
    }

    TA::TiledRange make_trange(const std::vector<basis<float>>& basis_sets,
                               const std::vector<size>& tile_sizes,
                               std::vector<TA::TiledRange1> ranges) {
        return _make_trange(basis_sets, tile_sizes, std::move(ranges));
    }

    template<typename basis_type>
    TA::TiledRange _make_trange(const std::vector<basis_type>& basis_sets,
                                const std::vector<std::vector<size>>& atom_ranges,
                                std::vector<TA::TiledRange1> ranges) {

        for (const auto& basis_set : basis_sets) {
            // figure out tilings for each basis set
            std::vector<size> tile_sizes = {};

            // Step through atom ranges [first, second]
            for (const auto& atom_range : atom_ranges) {
                size next_tile_size = 0;

                // Add the total basis set count for current atom range to the last bound
                for (auto atom : atom_range) {
                    next_tile_size += basis_set[atom].n_aos();
                }
                // Push back bound
                tile_sizes.push_back(next_tile_size);
            }

            // Make appropriate range for basis set
            ranges.push_back(make_tiled_range(basis_set, tile_sizes));
        }

        TA::TiledRange trange(ranges.begin(), ranges.end());
        return trange;
    }

    TA::TiledRange make_trange(const std::vector<basis<double>>& basis_sets,
                               const std::vector<std::vector<size>>& atom_ranges,
                               std::vector<TA::TiledRange1> ranges) {
        return _make_trange(basis_sets, atom_ranges, std::move(ranges));
    }

    TA::TiledRange make_trange(const std::vector<basis<float>>& basis_sets,
                               const std::vector<std::vector<size>>& atom_ranges,
                               std::vector<TA::TiledRange1> ranges) {
        return _make_trange(basis_sets, atom_ranges, std::move(ranges));
    }

} // namespace nwx_TA