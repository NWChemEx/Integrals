#include "fill_ND_functor.hpp"
#include "nwx_libint.hpp"

namespace nwx_TA {

template<typename val_type, libint2::Operator op, std::size_t NBases>
void FillNDFunctor<val_type, op, NBases>::initialize(const basis_vec& sets,
                                                     size_type deriv,
                                                     element_type thresh,
                                                     element_type cs) {
    LIBasis_sets = sets;
    cs_thresh    = (NBases == 2) ? 0.0 : cs;

    // Build factory
    factory.max_nprims = nwx_libint::sets_max_nprims(LIBasis_sets);
    factory.max_l      = nwx_libint::sets_max_l(LIBasis_sets);
    factory.thresh     = thresh;
    factory.deriv      = deriv;
}

template<typename val_type, libint2::Operator op, std::size_t NBases>
float FillNDFunctor<val_type, op, NBases>::operator()(
  val_type& tile, const TiledArray::Range& range) {
    // Determine if the entire tile is screened, before anymore work is done
    if((cs_thresh > 0.0) && screen.tile(LIBasis_sets, range, cs_thresh)) {
        return 0.0;
    }
    tile = val_type(range); // Initialize tile

    if(not libint2::initialized()) {
        libint2::initialize();
    }                                       // In case something else finalized
    auto tile_engine = factory(NBases, op); // Make libint engine
    size_vec offsets(
      NBases); // Vector for storing the offsets of the current shells
    size_vec shells(
      NBases); // Vector for storing the indices of the current shells

    std::vector<size_vec> tile_shells; // Shells in the current tile
    for(int depth = 0; depth < NBases; depth++) {
        int tile_depth    = (nopers == 1) ? depth : depth + 1;
        auto depth_shells = nwx_libint::aos2shells(LIBasis_sets[depth],
                                                   range.lobound()[tile_depth],
                                                   range.upbound()[tile_depth]);
        tile_shells.push_back(depth_shells);
    }

    // Start recurvise process to determine the shells needed for the current
    // tile
    _index_shells(tile, range, tile_engine, offsets, shells, tile_shells, 0);

    // Return norm for new tile
    return tile.norm();
}

template<typename val_type, libint2::Operator op, std::size_t NBases>
val_type FillNDFunctor<val_type, op, NBases>::operator()(
  const TiledArray::Range& range) {
    // Determine if the entire tile is screened, before anymore work is done
    if((cs_thresh > 0.0) && screen.tile(LIBasis_sets, range, cs_thresh)) {
        return val_type(range, 0.0);
    }
    val_type tile(range); // Declare and initialize tile

    if(not libint2::initialized()) {
        libint2::initialize();
    }                                       // In case something else finalized
    auto tile_engine = factory(NBases, op); // Make libint engine
    size_vec offsets(
      NBases); // Vector for storing the offsets of the current shells
    size_vec shells(
      NBases); // Vector for storing the indices of the current shells

    std::vector<size_vec> tile_shells; // Shells in the current tile; assumes
                                       // basis sets are tile specific
    for(int depth = 0; depth < NBases; depth++) {
        size_vec depth_shells;
        for(int i = 0; i < LIBasis_sets[depth].size(); ++i) {
            depth_shells.emplace_back(i);
        }
        tile_shells.push_back(depth_shells);
    }

    // Start recurvise process to determine the shells needed for the current
    // tile
    _index_shells(tile, range, tile_engine, offsets, shells, tile_shells, 0);

    // Return norm for new tile
    return tile;
}

template<typename val_type, libint2::Operator op, std::size_t NBases>
void FillNDFunctor<val_type, op, NBases>::_index_shells(
  val_type& tile, const TiledArray::Range& range, libint2::Engine& tile_engine,
  size_vec& offsets, size_vec& shells, std::vector<size_vec>& tile_shells,
  int depth) {
    // Deal with the additional dimension of the tile for multipoles
    int tile_depth = (nopers == 1) ? depth : depth + 1;

    // Set initial offset as the dimensions lower bound
    offsets[depth] = range.lobound()[tile_depth];

    // Loop over the shells that contain the AOs spanned by the tile in this
    // dimension
    for(auto s : tile_shells[depth]) {
        shells[depth] =
          s; // Save index of current shell to be passed down the line

        if(depth == (NBases - 1)) {
            // Check if the shell set can be screened
            auto screened =
              (cs_thresh == 0.0) ? false : screen.shellset(shells, cs_thresh);

            const auto& buf = tile_engine.results();
            if(not screened) { // Compute integrals for current shells
                _call_libint(tile_engine, shells,
                             std::make_index_sequence<NBases>());
            }

            // Switch for multipole cases
            if(nopers == 1) {
                size_vec indexer =
                  offsets; // Initially set indexer to first coordinate of tile
                const auto& int_vals =
                  (screened) ? nullptr :
                               buf[0]; // Store the location of the results
                int int_i = 0; // Index for integrals values; referenced so it
                               // gets pasted around

                // Fill in the values of the tile
                _fill_from_libint(tile, offsets, shells, int_vals, indexer,
                                  int_i, 0);
            } else {
                // Deal with the additional dimension of the tile for multipoles
                size_vec indexer = {range.lobound()[0]};
                indexer.insert(indexer.end(), offsets.begin(), offsets.end());

                // Loop over multipoles to fill in values
                for(auto f = 0ul; f != nopers; ++f) {
                    const auto& int_vals =
                      (screened) ? nullptr :
                                   buf[f]; // Store the location of the results
                    int int_i = 0; // Index for integrals values; referenced so
                                   // it gets pasted around

                    // Fill in the values of the multipole
                    _fill_from_libint(tile, offsets, shells, int_vals, indexer,
                                      int_i, 0);

                    indexer[0]++; // Increment the coordinate for the multipoles
                }
            }
        } else {
            // Keep going down until all dimensions of the the tensor have been
            // covered
            _index_shells(tile, range, tile_engine, offsets, shells,
                          tile_shells, depth + 1);
        }
        // Increase the offsets for the current dimension by the size of the
        // completed shell
        offsets[depth] += LIBasis_sets[depth][shells[depth]].size();
    }
}

template<typename val_type, libint2::Operator op, std::size_t NBases>
void FillNDFunctor<val_type, op, NBases>::_fill_from_libint(
  val_type& tile, const size_vec& offsets, const size_vec& shells,
  const double* int_vals, size_vec& indexer, int& int_i, int depth) {
    // Deal with the additional dimension of the tile for multipoles
    int tile_depth = (nopers == 1) ? depth : depth + 1;

    // Loop over the size of the current shell for this dimension
    for(auto f = 0ul; f != LIBasis_sets[depth][shells[depth]].size(); ++f) {
        if(depth == (NBases - 1)) {
            // Assign integral value
            if(int_vals == nullptr) {
                tile[indexer] = 0; // Default case of all zeroes
            } else {
                tile[indexer] = int_vals[int_i]; // Get values from LibInt
                int_i++;
            }
        } else {
            // Keep going down until all dimensions of the the tensor have been
            // covered
            _fill_from_libint(tile, offsets, shells, int_vals, indexer, int_i,
                              depth + 1);
        }
        indexer[tile_depth]++; // Increment the coordinate for this dimension
    }
    indexer[tile_depth] =
      offsets[depth]; // Reset the indexer to the initial position for next loop
}

template<typename val_type, libint2::Operator op, std::size_t NBases>
template<std::size_t... Is>
void FillNDFunctor<val_type, op, NBases>::_call_libint(
  libint2::Engine& tile_engine, size_vec shells, std::index_sequence<Is...>) {
    tile_engine.compute(LIBasis_sets[Is][shells[Is]]...);
}

template class FillNDFunctor<TA::TensorD, libint2::Operator::overlap, 2>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::kinetic, 2>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::nuclear, 2>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::coulomb, 2>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::coulomb, 3>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::coulomb, 4>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::stg, 2>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::stg, 3>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::stg, 4>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::yukawa, 2>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::yukawa, 3>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::yukawa, 4>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::emultipole1, 2>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::emultipole2, 2>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::emultipole3, 2>;
template class FillNDFunctor<TA::TensorD, libint2::Operator::delta, 4>;
} // namespace nwx_TA