/*
 * Copyright 2022 NWChemEx-Project
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
#include <libint2.hpp>
#include <simde/types.hpp>

namespace integrals::ao_integrals::detail_ {

/** @brief Find the ordinal indices spanned by the shell
 *
 *  @param[in] bases The basis sets.
 *  @param[in] shell The shell indices
 *  @returns A std::vector of the ordinal indices associated with the shell.
 */
inline auto shells2ord(const std::vector<libint2::BasisSet>& bases,
                       std::vector<std::size_t>& shells) {
    using size_vector_t = std::vector<std::size_t>;

    // Calculate the ordinal step of each basis dimension other than the Nth.
    auto N = bases.size();
    size_vector_t step_sizes(N - 1);
    for(decltype(N) i = 1; i <= N - 1; ++i) {
        step_sizes[N - 1 - i] = bases[N - i].nbf();
        if(i > 1) step_sizes[N - 1 - i] *= step_sizes[N - i];
    }

    // Set starting and ending AO indices
    size_vector_t curr_ao(N);
    size_vector_t up_ao(N);
    for(decltype(N) i = 0; i < N; ++i) {
        curr_ao[i] = bases[i].shell2bf()[shells[i]];
        up_ao[i]   = curr_ao[i] + bases[i][shells[i]].size();
    }

    // Increment through the AO indices of the shell and determine the ordinal
    // index for each.
    size_vector_t ords;
    while(curr_ao[0] < up_ao[0]) {
        // ordinal calculation
        auto curr_ord = curr_ao.back();
        for(std::size_t i = 0; i < step_sizes.size(); ++i) {
            curr_ord += curr_ao[i] * step_sizes[i];
        }
        ords.push_back(curr_ord);

        // Increment to the next AO index.
        curr_ao.back() += 1;
        for(decltype(N) i = 1; i < N; ++i) {
            if(curr_ao[N - i] >= up_ao[N - i]) {
                // curr_ao[0] accumalates until it passes up_aos[0]
                // and the loop terminates.
                curr_ao[N - i] = bases[N - i].shell2bf()[shells[N - i]];
                curr_ao[N - i - 1] += 1;
            }
        }
    }
    return ords;
}

} // namespace integrals::ao_integrals::detail_
