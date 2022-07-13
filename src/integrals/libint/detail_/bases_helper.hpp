#pragma once
#include <array>
#include <chemist/enums.hpp> /// For ShellType
#include <libint2.hpp>
#include <simde/types.hpp>
#include <vector>

namespace integrals::detail_ {

/** @brief Converts an NWX basis set object to a LibInt2 basis set object.
 *
 *  @param[in] bs The NWX basis set to be converted
 *  @returns The basis set as a LibInt2 basis set
 */
inline auto make_libint_basis_set(const simde::type::ao_basis_set& bs) {
    using svec_d_t = libint2::svector<double>;
    using cont_t   = libint2::Shell::Contraction;

    libint2::BasisSet shells;
    for(const auto&& shelli : bs.shells()) {
        const auto nprims     = shelli.n_unique_primitives();
        const auto first_prim = shelli.unique_primitive(0);
        const auto last_prim  = shelli.unique_primitive(nprims - 1);
        const bool pure       = shelli.pure() == chemist::ShellType::pure;
        const int l           = shelli.l();

        svec_d_t alphas(&first_prim.exponent(), &last_prim.exponent() + 1);
        svec_d_t coefs(&first_prim.coefficient(), &last_prim.coefficient() + 1);
        libint2::svector<cont_t> conts{cont_t{l, pure, coefs}};
        std::array<double, 3> center = {shelli.x(), shelli.y(), shelli.z()};
        shells.push_back(libint2::Shell(alphas, conts, center));
    }
    return shells;
}

/** @brief Unpacks the basis sets from the inputs and converts them to Libint.
 *
 *  @param[in] inputs The module inputs containing the basis sets.
 *  @returns A vector of the converted basis sets.
 */
template<std::size_t N, typename ModuleInputs>
auto unpack_bases(const ModuleInputs& inputs) {
    using ao_space_t = simde::type::ao_space;
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
    std::vector<libint2::BasisSet> rv;
    for(auto i = 0u; i < N; ++i)
        rv.emplace_back(make_libint_basis_set(aos[i].basis_set()));
    return rv;
}

} // namespace integrals::detail_
