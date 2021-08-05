#pragma once
#include <simde/types.hpp>
#include <vector>

namespace integrals::detail_ {

template<std::size_t N, typename ModuleInputs>
auto unpack_bases(const ModuleInputs& inputs) {
    using ao_space_t   = simde::type::ao_space;
    using ao_space_ref = const ao_space_t&;
    std::vector<ao_space_t> aos(N);
    if constexpr(N == 2) {
        aos[0] = inputs.at("bra").template value<ao_space_t>();
        aos[1] = inputs.at("ket").template value<ao_space_t>();
    } else if constexpr(N == 3) {
        aos[0] = inputs.at("bra").template value<ao_space_t>();
        aos[1] = inputs.at("ket 1").template value<ao_space_t>();
        aos[2] = inputs.at("ket 2").template value<ao_space_t>();
    } else if constexpr(N == 4) {
        aos[0] = inputs.at("bra 1").template value<ao_space_t>();
        aos[1] = inputs.at("bra 2").template value<ao_space_t>();
        aos[2] = inputs.at("ket 1").template value<ao_space_t>();
        aos[3] = inputs.at("ket 2").template value<ao_space_t>();
    }
    std::vector<simde::type::ao_basis_set> rv;
    for(auto i = 0u; i < N; ++i) rv.emplace_back(aos[i].basis_set());
    return rv;
}

} // namespace integrals::detail_
