#include "integrals/nwx_libint/cauchy_schwarz_screener.hpp"
#include "integrals/nwx_libint/nwx_libint.hpp"

namespace nwx_libint {

template<std::size_t NBases>
bool CauchySchwarzScreener<NBases>::tile(const basis_vec& basis_sets,
                                         const TiledArray::Range& range,
                                         double cs_thresh) {
    // Shell lists for current tile
    std::vector<size_vec> shell_list;

    // Fill in the above vectors
    for(int i = 0; i < NBases; ++i) {
        shell_list.push_back(nwx_libint::aos2shells(
          basis_sets[i], range.lobound()[i], range.upbound()[i]));
    }

    // Find the largest coefficient in each submatrix and use them for the
    // approximation
    if constexpr(NBases == 2) {
        return false;
    } else if constexpr(NBases == 3) {
        double maxVal = 0.0;

        for(auto i : shell_list[1]) {
            for(auto j : shell_list[2]) {
                maxVal = std::max(maxVal, cs_mat2[i][j]);
            }
        }

        return maxVal < cs_thresh;
    } else if constexpr(NBases == 4) {
        double maxVal1 = 0.0, maxVal2 = 0.0;

        for(auto i : shell_list[0]) {
            for(auto j : shell_list[1]) {
                maxVal1 = std::max(maxVal1, cs_mat1[i][j]);
            }
        }
        for(auto i : shell_list[2]) {
            for(auto j : shell_list[3]) {
                maxVal2 = std::max(maxVal2, cs_mat2[i][j]);
            }
        }

        return (maxVal1 * maxVal2) < cs_thresh;
    }
}

template<std::size_t NBases>
bool CauchySchwarzScreener<NBases>::shellset(size_vec shells,
                                             double cs_thresh) {
    // Check approximation product vs provided threshold value
    if constexpr(NBases == 2) {
        return false;
    } else if constexpr(NBases == 3) {
        return cs_mat2[shells[1]][shells[2]] < cs_thresh;
    } else if constexpr(NBases == 4) {
        return (cs_mat1[shells[0]][shells[1]] * cs_mat2[shells[2]][shells[3]]) <
               cs_thresh;
    }
    return false;
}

template<std::size_t NBases>
void CauchySchwarzScreener<NBases>::set_sub_screen(
  const basis_vec& basis_sets, const TiledArray::Range& range, approx_vec& mat1,
  approx_vec& mat2) const {
    if(cs_mat2.empty()) {
        // Default out if cs_mats aren't set to anything
        mat1 = {}, mat2 = {};
    } else {
        // Shell lists for current tile
        std::vector<size_vec> shell_list;

        // Fill in the above vectors
        for(int i = 0; i < NBases; ++i) {
            shell_list.push_back(nwx_libint::aos2shells(
              basis_sets[i], range.lobound()[i], range.upbound()[i]));
        }

        if constexpr(NBases == 3) {
            mat1 = {};
            mat2 = approx_vec(shell_list[1].size(),
                              double_vec(shell_list[2].size(), 0.0));
            for(int i = 0; i < shell_list[1].size(); ++i) {
                for(int j = 0; j < shell_list[2].size(); ++j) {
                    mat2[i][j] = cs_mat2[shell_list[1][i]][shell_list[2][j]];
                }
            }
        } else if constexpr(NBases == 4) {
            mat1 = approx_vec(shell_list[0].size(),
                              double_vec(shell_list[1].size(), 0.0));
            mat2 = approx_vec(shell_list[2].size(),
                              double_vec(shell_list[3].size(), 0.0));
            for(int i = 0; i < shell_list[0].size(); ++i) {
                for(int j = 0; j < shell_list[1].size(); ++j) {
                    mat1[i][j] = cs_mat1[shell_list[0][i]][shell_list[1][j]];
                }
            }
            for(int i = 0; i < shell_list[2].size(); ++i) {
                for(int j = 0; j < shell_list[3].size(); ++j) {
                    mat2[i][j] = cs_mat2[shell_list[2][i]][shell_list[3][j]];
                }
            }
        }
    }
}

template class CauchySchwarzScreener<2>;
template class CauchySchwarzScreener<3>;
template class CauchySchwarzScreener<4>;

} // namespace nwx_libint