#pragma once
#include <libint2.hpp>

namespace testing {

/// Water STO-3G basis set in Libint Format
inline libint2::BasisSet water_basis_set() {
    libint2::BasisSet bset{};
    bset.push_back(
      libint2::Shell{{130.7093200, 23.8088610, 6.4436083},
                     {{0, true, {0.15432897, 0.53532814, 0.44463454}}},
                     {{0.0, -0.143222342980786, 0.0}}});
    bset.push_back(
      libint2::Shell{{5.0331513, 1.1695961, 0.3803890},
                     {{0, true, {-0.09996723, 0.39951283, 0.70011547}}},
                     {{0.0, -0.143222342980786, 0.0}}});
    bset.push_back(
      libint2::Shell{{5.0331513, 1.1695961, 0.3803890},
                     {{1, true, {0.15591627, 0.60768372, 0.39195739}}},
                     {{0.0, -0.143222342980786, 0.0}}});
    bset.push_back(
      libint2::Shell{{3.425250914, 0.6239137298, 0.1688554040},
                     {{0, true, {0.1543289673, 0.5353281423, 0.4446345422}}},
                     {{1.638033502034240, 1.136556880358410, 0.0}}});
    bset.push_back(
      libint2::Shell{{3.425250914, 0.6239137298, 0.1688554040},
                     {{0, true, {0.1543289673, 0.5353281423, 0.4446345422}}},
                     {{-1.638033502034240, 1.136556880358410, 0.0}}});
    return bset;
}

} // namespace testing