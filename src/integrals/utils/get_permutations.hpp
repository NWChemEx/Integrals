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

#pragma once
#include <array>
#include <map>
#include <set>
#include <type_traits>
#include <vector>

namespace integrals::utils {

/** @brief Works out all the permutations of index @p abcd that are symmetry
 *         equivalent.
 *
 *  @return std::set of unique permuted arrays (original API, for backward
 *          compatibility with existing callers and unit tests).
 */
template<typename T>
auto get_permutations(T&& abcd, bool p0_1, bool p2_3, bool p01_23) {
    using element_type = std::decay_t<decltype(abcd[0])>;
    using index_type   = std::array<element_type, 4>;

    std::set<index_type> permutations{abcd};

    if(p0_1) permutations.insert({abcd[1], abcd[0], abcd[2], abcd[3]});
    if(p2_3) permutations.insert({abcd[0], abcd[1], abcd[3], abcd[2]});
    if(p01_23) permutations.insert({abcd[2], abcd[3], abcd[0], abcd[1]});
    if(p0_1 && p01_23)
        permutations.insert({abcd[2], abcd[3], abcd[1], abcd[0]});
    if(p2_3 && p01_23)
        permutations.insert({abcd[3], abcd[2], abcd[0], abcd[1]});
    if(p0_1 && p2_3) permutations.insert({abcd[1], abcd[0], abcd[3], abcd[2]});
    if(p0_1 && p2_3 && p01_23)
        permutations.insert({abcd[3], abcd[2], abcd[1], abcd[0]});

    return permutations;
}

/** @brief Like get_permutations but also returns the index permutation sigma
 *         for each result, where sigma[d] is the canonical dimension that
 *         permuted dimension d came from: perm[d] == abcd[sigma[d]].
 *
 *  Deduplicates by permuted value (first sigma wins), so repeated values in
 *  @p abcd are handled correctly.
 *
 *  @return vector of (permuted_abcd, sigma) pairs.
 */
template<typename T>
auto get_permutations_with_sigma(T&& abcd, bool p0_1, bool p2_3, bool p01_23) {
    using element_type = std::decay_t<decltype(abcd[0])>;
    using index_type   = std::array<element_type, 4>;
    using sigma_type   = std::array<std::size_t, 4>;
    using pair_type    = std::pair<index_type, sigma_type>;

    std::map<index_type, sigma_type> seen;
    auto try_insert = [&](index_type perm, sigma_type sigma) {
        seen.emplace(perm, sigma);
    };

    try_insert({abcd[0], abcd[1], abcd[2], abcd[3]}, {0, 1, 2, 3});
    if(p0_1) try_insert({abcd[1], abcd[0], abcd[2], abcd[3]}, {1, 0, 2, 3});
    if(p2_3) try_insert({abcd[0], abcd[1], abcd[3], abcd[2]}, {0, 1, 3, 2});
    if(p01_23) try_insert({abcd[2], abcd[3], abcd[0], abcd[1]}, {2, 3, 0, 1});
    if(p0_1 && p01_23)
        try_insert({abcd[2], abcd[3], abcd[1], abcd[0]}, {2, 3, 1, 0});
    if(p2_3 && p01_23)
        try_insert({abcd[3], abcd[2], abcd[0], abcd[1]}, {3, 2, 0, 1});
    if(p0_1 && p2_3)
        try_insert({abcd[1], abcd[0], abcd[3], abcd[2]}, {1, 0, 3, 2});
    if(p0_1 && p2_3 && p01_23)
        try_insert({abcd[3], abcd[2], abcd[1], abcd[0]}, {3, 2, 1, 0});

    std::vector<pair_type> result;
    result.reserve(seen.size());
    for(auto& [perm, sigma] : seen) result.emplace_back(perm, sigma);
    return result;
}

} // namespace integrals::utils
