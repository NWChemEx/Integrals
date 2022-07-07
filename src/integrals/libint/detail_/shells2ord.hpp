#pragma once
#include <libint2.hpp>

namespace integrals::detail_ {

/** @brief Find the ordinal indices spanned by the shell indices
 *
 *  @param[in] bases A vector of libint basis sets.
 *  @param[in] shell A vector of indices for shells in the basis sets.
 *  @param[in] lo_shell The lower most shell index in the span.
 *  @param[in] up_shell The upper most shell index in the span.
 *  @returns An std::vector of the ordinal indices associated with the shells.
 */
inline auto shells2ord(const std::vector<libint2::BasisSet>& bases,
                       std::vector<std::size_t> shells,
                       std::vector<std::size_t> lo_shells,
                       std::vector<std::size_t> up_shells) {
    using size_vector_t = std::vector<std::size_t>;

    /// Need bases extents and the AO span of the shells
    /// for ordinal position determination
    auto N = bases.size();
    size_vector_t extents;
    size_vector_t lo_ao;
    size_vector_t up_ao;
    for(auto i = 0; i < N; ++i) {
        std::size_t set_extent = 0;
        for(auto j = lo_shells[i]; j <= up_shells[i]; ++j) {
            if(j == shells[i]) {
                lo_ao.push_back(set_extent);
                up_ao.push_back(set_extent + bases[i][j].size() - 1);
            }
            set_extent += bases[i][j].size();
        }
        extents.push_back(set_extent);
    }

    /// Calculate the ordinal step of each basis dimension other than the Nth.
    size_vector_t step_size;
    for(auto i = 0; i < N - 1; ++i) {
        auto step = extents[N - 1];
        for(auto j = i + 1; j < N - 1; ++j) step *= extents[j];
        step_size.push_back(step);
    }

    /// Increment through the AO indices of the shell and determine the ordinal
    /// index for each.
    size_vector_t curr_ao = lo_ao;
    size_vector_t ord_pos;
    while(curr_ao[0] <= up_ao[0]) {
        /// ordinal calculation
        auto curr_ord_pos = curr_ao[N - 1];
        for(auto i = 0; i < N - 1; ++i) {
            curr_ord_pos += curr_ao[i] * step_size[i];
        }
        ord_pos.push_back(curr_ord_pos);

        /// Increment to the next AO index.
        curr_ao[N - 1] += 1;
        for(auto i = 1; i < N; ++i) {
            if(curr_ao[N - i] > up_ao[N - i]) {
                /// curr_ao[0] accumalates until it passes up_aos[0]
                /// and the loop terminates.
                curr_ao[N - i] = lo_ao[N - i];
                curr_ao[N - i - 1] += 1;
            }
        }
    }

    return ord_pos;
}

} // namespace integrals::detail_