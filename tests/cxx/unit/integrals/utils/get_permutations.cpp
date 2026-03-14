/*
 * Copyright 2026 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "../testing/testing.hpp"
#include <integrals/utils/get_permutations.hpp>

using namespace integrals::utils;

/* Testing Notes:
 *
 * We assume that std::set works and thus won't allow us to insert duplicates.
 * Thus a spot check for an input with all unique values, and one with a
 * repeated index should suffice. This leaves 2^3 = 8 combinations of
 * symmetry to check.
 */

TEST_CASE("get_permutations") {
    using type = std::array<std::size_t, 4>;
    SECTION("No symmetry") {
        SECTION("All unique") {
            auto perms =
              get_permutations(type{0, 1, 2, 3}, false, false, false);
            REQUIRE(perms.size() == 1);
            REQUIRE(perms.count({0, 1, 2, 3}) == 1);
        }
        SECTION("Three repeats") {
            auto perms =
              get_permutations(type{0, 0, 0, 1}, false, false, false);
            REQUIRE(perms.size() == 1);
            REQUIRE(perms.count({0, 0, 0, 1}) == 1);
        }
    }
    SECTION("p0_1 symmetry") {
        SECTION("All unique") {
            auto perms = get_permutations(type{0, 1, 2, 3}, true, false, false);
            REQUIRE(perms.size() == 2);
            REQUIRE(perms.count({0, 1, 2, 3}) == 1);
            REQUIRE(perms.count({1, 0, 2, 3}) == 1);
        }
        SECTION("One repeat (in symmetry)") {
            auto perms = get_permutations(type{0, 0, 1, 2}, true, false, false);
            REQUIRE(perms.size() == 1);
            REQUIRE(perms.count({0, 0, 1, 2}) == 1);
        }
        SECTION("One repeat (not in symmetry)") {
            auto perms = get_permutations(type{0, 1, 1, 2}, true, false, false);
            REQUIRE(perms.size() == 2);
            REQUIRE(perms.count({0, 1, 1, 2}) == 1);
            REQUIRE(perms.count({1, 0, 1, 2}) == 1);
        }
    }

    SECTION("p2_3 symmetry") {
        // Conceptually very similar to p0_1 symmetry so just spot check
        auto perms = get_permutations(type{0, 1, 2, 3}, false, true, false);
        REQUIRE(perms.size() == 2);
        REQUIRE(perms.count({0, 1, 2, 3}) == 1);
        REQUIRE(perms.count({0, 1, 3, 2}) == 1);
    }

    SECTION("p01_23") {
        // Conceptually very similar to p0_1 symmetry so just spot check
        auto perms = get_permutations(type{0, 1, 2, 3}, false, false, true);
        REQUIRE(perms.size() == 2);
        REQUIRE(perms.count({0, 1, 2, 3}) == 1);
        REQUIRE(perms.count({2, 3, 0, 1}) == 1);
    }

    SECTION("P0_1 and p2_3") {
        SECTION("All unique") {
            auto perms = get_permutations(type{0, 1, 2, 3}, true, true, false);
            REQUIRE(perms.size() == 4);
            REQUIRE(perms.count({0, 1, 2, 3}) == 1);
            REQUIRE(perms.count({1, 0, 2, 3}) == 1);
            REQUIRE(perms.count({0, 1, 3, 2}) == 1);
            REQUIRE(perms.count({1, 0, 3, 2}) == 1);
        }
        SECTION("One repeat") {
            auto perms = get_permutations(type{0, 0, 1, 2}, true, true, false);
            REQUIRE(perms.size() == 2);
            REQUIRE(perms.count({0, 0, 1, 2}) == 1);
            REQUIRE(perms.count({0, 0, 2, 1}) == 1);
        }
        SECTION("Two repeats (one in each symmetry group)") {
            auto perms = get_permutations(type{0, 0, 1, 1}, true, true, false);
            REQUIRE(perms.size() == 1);
            REQUIRE(perms.count({0, 0, 1, 1}) == 1);
        }
    }
    SECTION("P0_1 and p01_23") {
        // Conceptually very similar to p0_1 && p2_3 symmetry so just spot check
        auto perms = get_permutations(type{0, 1, 2, 3}, true, false, true);
        REQUIRE(perms.size() == 4);
        REQUIRE(perms.count({0, 1, 2, 3}) == 1);
        REQUIRE(perms.count({1, 0, 2, 3}) == 1);
        REQUIRE(perms.count({2, 3, 0, 1}) == 1);
        REQUIRE(perms.count({2, 3, 1, 0}) == 1);
    }
    SECTION("P2_3 and p01_23") {
        // Conceptually very similar to p0_1 && p2_3 symmetry so just spot check
        auto perms = get_permutations(type{0, 1, 2, 3}, false, true, true);
        REQUIRE(perms.size() == 4);
        REQUIRE(perms.count({0, 1, 2, 3}) == 1);
        REQUIRE(perms.count({0, 1, 3, 2}) == 1);
        REQUIRE(perms.count({2, 3, 0, 1}) == 1);
        REQUIRE(perms.count({3, 2, 0, 1}) == 1);
    }
    SECTION("P0_1 and p2_3 and p01_23") {
        SECTION("All unique") {
            auto perms = get_permutations(type{0, 1, 2, 3}, true, true, true);
            REQUIRE(perms.size() == 8);
            REQUIRE(perms.count({0, 1, 2, 3}) == 1);
            REQUIRE(perms.count({1, 0, 2, 3}) == 1);
            REQUIRE(perms.count({0, 1, 3, 2}) == 1);
            REQUIRE(perms.count({1, 0, 3, 2}) == 1);
            REQUIRE(perms.count({2, 3, 0, 1}) == 1);
            REQUIRE(perms.count({3, 2, 0, 1}) == 1);
            REQUIRE(perms.count({2, 3, 1, 0}) == 1);
            REQUIRE(perms.count({3, 2, 1, 0}) == 1);
        }
        SECTION("One repeat") {
            auto perms = get_permutations(type{0, 0, 1, 2}, true, true, true);
            REQUIRE(perms.size() == 4);
            REQUIRE(perms.count({0, 0, 1, 2}) == 1);
            REQUIRE(perms.count({0, 0, 2, 1}) == 1);
            REQUIRE(perms.count({1, 2, 0, 0}) == 1);
            REQUIRE(perms.count({2, 1, 0, 0}) == 1);
        }
    }
}

/* Testing Notes for get_permutations_with_sigma:
 *
 * For every returned (perm, sigma) pair the invariant perm[d] == abcd[sigma[d]]
 * must hold for all d in {0,1,2,3}.  We verify this invariant explicitly for
 * every pair in every test case, in addition to checking the expected set of
 * permuted arrays and their counts.
 *
 * Deduplication behaviour: when two symmetry operations produce the same
 * permuted array (because abcd has repeated values) only one entry is kept and
 * its sigma is the one that was inserted first (identity wins over later
 * operations).
 */

namespace {
// Helper: verify perm[d] == abcd[sigma[d]] for all d
template<typename Index, typename Sigma>
bool sigma_consistent(const Index& abcd, const Index& perm,
                      const Sigma& sigma) {
    for(std::size_t d = 0; d < 4; ++d)
        if(perm[d] != abcd[sigma[d]]) return false;
    return true;
}

// Helper: find a (perm, sigma) pair by its perm value
template<typename Pairs, typename Index>
auto find_perm(const Pairs& pairs, const Index& target) {
    for(auto& [p, s] : pairs)
        if(p == target) return s;
    // Return identity sigma as sentinel — test will fail on sigma check
    return std::array<std::size_t, 4>{0, 1, 2, 3};
}

// Helper: count how many pairs have a given perm value
template<typename Pairs, typename Index>
std::size_t count_perm(const Pairs& pairs, const Index& target) {
    std::size_t n = 0;
    for(auto& [p, s] : pairs)
        if(p == target) ++n;
    return n;
}
} // namespace

TEST_CASE("get_permutations_with_sigma") {
    using type       = std::array<std::size_t, 4>;
    using sigma_type = std::array<std::size_t, 4>;

    SECTION("No symmetry") {
        SECTION("All unique") {
            type abcd{0, 1, 2, 3};
            auto pairs = get_permutations_with_sigma(abcd, false, false, false);
            REQUIRE(pairs.size() == 1);
            REQUIRE(count_perm(pairs, type{0, 1, 2, 3}) == 1);
            auto s = find_perm(pairs, type{0, 1, 2, 3});
            REQUIRE(s == sigma_type{0, 1, 2, 3});
            for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
        }
        SECTION("Three repeats") {
            type abcd{0, 0, 0, 1};
            auto pairs = get_permutations_with_sigma(abcd, false, false, false);
            REQUIRE(pairs.size() == 1);
            REQUIRE(count_perm(pairs, type{0, 0, 0, 1}) == 1);
            for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
        }
    }

    SECTION("p0_1 symmetry") {
        SECTION("All unique") {
            type abcd{0, 1, 2, 3};
            auto pairs = get_permutations_with_sigma(abcd, true, false, false);
            REQUIRE(pairs.size() == 2);
            REQUIRE(count_perm(pairs, type{0, 1, 2, 3}) == 1);
            REQUIRE(count_perm(pairs, type{1, 0, 2, 3}) == 1);
            REQUIRE(find_perm(pairs, type{0, 1, 2, 3}) ==
                    sigma_type{0, 1, 2, 3});
            REQUIRE(find_perm(pairs, type{1, 0, 2, 3}) ==
                    sigma_type{1, 0, 2, 3});
            for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
        }
        SECTION("One repeat (in symmetry) — deduplicates to 1 entry") {
            type abcd{0, 0, 1, 2};
            auto pairs = get_permutations_with_sigma(abcd, true, false, false);
            REQUIRE(pairs.size() == 1);
            REQUIRE(count_perm(pairs, type{0, 0, 1, 2}) == 1);
            // Identity sigma wins
            REQUIRE(find_perm(pairs, type{0, 0, 1, 2}) ==
                    sigma_type{0, 1, 2, 3});
            for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
        }
        SECTION("One repeat (not in symmetry)") {
            type abcd{0, 1, 1, 2};
            auto pairs = get_permutations_with_sigma(abcd, true, false, false);
            REQUIRE(pairs.size() == 2);
            REQUIRE(count_perm(pairs, type{0, 1, 1, 2}) == 1);
            REQUIRE(count_perm(pairs, type{1, 0, 1, 2}) == 1);
            REQUIRE(find_perm(pairs, type{0, 1, 1, 2}) ==
                    sigma_type{0, 1, 2, 3});
            REQUIRE(find_perm(pairs, type{1, 0, 1, 2}) ==
                    sigma_type{1, 0, 2, 3});
            for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
        }
    }

    SECTION("p2_3 symmetry") {
        type abcd{0, 1, 2, 3};
        auto pairs = get_permutations_with_sigma(abcd, false, true, false);
        REQUIRE(pairs.size() == 2);
        REQUIRE(count_perm(pairs, type{0, 1, 2, 3}) == 1);
        REQUIRE(count_perm(pairs, type{0, 1, 3, 2}) == 1);
        REQUIRE(find_perm(pairs, type{0, 1, 2, 3}) == sigma_type{0, 1, 2, 3});
        REQUIRE(find_perm(pairs, type{0, 1, 3, 2}) == sigma_type{0, 1, 3, 2});
        for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
    }

    SECTION("p01_23 symmetry") {
        type abcd{0, 1, 2, 3};
        auto pairs = get_permutations_with_sigma(abcd, false, false, true);
        REQUIRE(pairs.size() == 2);
        REQUIRE(count_perm(pairs, type{0, 1, 2, 3}) == 1);
        REQUIRE(count_perm(pairs, type{2, 3, 0, 1}) == 1);
        REQUIRE(find_perm(pairs, type{0, 1, 2, 3}) == sigma_type{0, 1, 2, 3});
        REQUIRE(find_perm(pairs, type{2, 3, 0, 1}) == sigma_type{2, 3, 0, 1});
        for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
    }

    SECTION("p0_1 and p2_3") {
        SECTION("All unique") {
            type abcd{0, 1, 2, 3};
            auto pairs = get_permutations_with_sigma(abcd, true, true, false);
            REQUIRE(pairs.size() == 4);
            REQUIRE(count_perm(pairs, type{0, 1, 2, 3}) == 1);
            REQUIRE(count_perm(pairs, type{1, 0, 2, 3}) == 1);
            REQUIRE(count_perm(pairs, type{0, 1, 3, 2}) == 1);
            REQUIRE(count_perm(pairs, type{1, 0, 3, 2}) == 1);
            REQUIRE(find_perm(pairs, type{0, 1, 2, 3}) ==
                    sigma_type{0, 1, 2, 3});
            REQUIRE(find_perm(pairs, type{1, 0, 2, 3}) ==
                    sigma_type{1, 0, 2, 3});
            REQUIRE(find_perm(pairs, type{0, 1, 3, 2}) ==
                    sigma_type{0, 1, 3, 2});
            REQUIRE(find_perm(pairs, type{1, 0, 3, 2}) ==
                    sigma_type{1, 0, 3, 2});
            for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
        }
        SECTION("One repeat — deduplicates") {
            type abcd{0, 0, 1, 2};
            auto pairs = get_permutations_with_sigma(abcd, true, true, false);
            REQUIRE(pairs.size() == 2);
            REQUIRE(count_perm(pairs, type{0, 0, 1, 2}) == 1);
            REQUIRE(count_perm(pairs, type{0, 0, 2, 1}) == 1);
            for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
        }
        SECTION(
          "Two repeats (one in each symmetry group) — deduplicates to 1") {
            type abcd{0, 0, 1, 1};
            auto pairs = get_permutations_with_sigma(abcd, true, true, false);
            REQUIRE(pairs.size() == 1);
            REQUIRE(count_perm(pairs, type{0, 0, 1, 1}) == 1);
            for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
        }
    }

    SECTION("p0_1 and p01_23") {
        type abcd{0, 1, 2, 3};
        auto pairs = get_permutations_with_sigma(abcd, true, false, true);
        REQUIRE(pairs.size() == 4);
        REQUIRE(count_perm(pairs, type{0, 1, 2, 3}) == 1);
        REQUIRE(count_perm(pairs, type{1, 0, 2, 3}) == 1);
        REQUIRE(count_perm(pairs, type{2, 3, 0, 1}) == 1);
        REQUIRE(count_perm(pairs, type{2, 3, 1, 0}) == 1);
        REQUIRE(find_perm(pairs, type{0, 1, 2, 3}) == sigma_type{0, 1, 2, 3});
        REQUIRE(find_perm(pairs, type{1, 0, 2, 3}) == sigma_type{1, 0, 2, 3});
        REQUIRE(find_perm(pairs, type{2, 3, 0, 1}) == sigma_type{2, 3, 0, 1});
        REQUIRE(find_perm(pairs, type{2, 3, 1, 0}) == sigma_type{2, 3, 1, 0});
        for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
    }

    SECTION("p2_3 and p01_23") {
        type abcd{0, 1, 2, 3};
        auto pairs = get_permutations_with_sigma(abcd, false, true, true);
        REQUIRE(pairs.size() == 4);
        REQUIRE(count_perm(pairs, type{0, 1, 2, 3}) == 1);
        REQUIRE(count_perm(pairs, type{0, 1, 3, 2}) == 1);
        REQUIRE(count_perm(pairs, type{2, 3, 0, 1}) == 1);
        REQUIRE(count_perm(pairs, type{3, 2, 0, 1}) == 1);
        REQUIRE(find_perm(pairs, type{0, 1, 2, 3}) == sigma_type{0, 1, 2, 3});
        REQUIRE(find_perm(pairs, type{0, 1, 3, 2}) == sigma_type{0, 1, 3, 2});
        REQUIRE(find_perm(pairs, type{2, 3, 0, 1}) == sigma_type{2, 3, 0, 1});
        REQUIRE(find_perm(pairs, type{3, 2, 0, 1}) == sigma_type{3, 2, 0, 1});
        for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
    }

    SECTION("p0_1 and p2_3 and p01_23") {
        SECTION("All unique — all 8 permutations") {
            type abcd{0, 1, 2, 3};
            auto pairs = get_permutations_with_sigma(abcd, true, true, true);
            REQUIRE(pairs.size() == 8);
            const std::vector<std::pair<type, sigma_type>> expected{
              {{0, 1, 2, 3}, {0, 1, 2, 3}}, {{1, 0, 2, 3}, {1, 0, 2, 3}},
              {{0, 1, 3, 2}, {0, 1, 3, 2}}, {{1, 0, 3, 2}, {1, 0, 3, 2}},
              {{2, 3, 0, 1}, {2, 3, 0, 1}}, {{2, 3, 1, 0}, {2, 3, 1, 0}},
              {{3, 2, 0, 1}, {3, 2, 0, 1}}, {{3, 2, 1, 0}, {3, 2, 1, 0}},
            };
            for(auto& [exp_perm, exp_sigma] : expected) {
                REQUIRE(count_perm(pairs, exp_perm) == 1);
                REQUIRE(find_perm(pairs, exp_perm) == exp_sigma);
            }
            for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
        }
        SECTION("One repeat — deduplicates to 4 entries") {
            type abcd{0, 0, 1, 2};
            auto pairs = get_permutations_with_sigma(abcd, true, true, true);
            REQUIRE(pairs.size() == 4);
            REQUIRE(count_perm(pairs, type{0, 0, 1, 2}) == 1);
            REQUIRE(count_perm(pairs, type{0, 0, 2, 1}) == 1);
            REQUIRE(count_perm(pairs, type{1, 2, 0, 0}) == 1);
            REQUIRE(count_perm(pairs, type{2, 1, 0, 0}) == 1);
            for(auto& [p, sig] : pairs) REQUIRE(sigma_consistent(abcd, p, sig));
        }
    }
}
