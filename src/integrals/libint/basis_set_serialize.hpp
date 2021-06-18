#pragma once
#include <libint2.hpp>
#include <tiledarray.h>

namespace nwx_TA {

/** @brief Serialization of a Libint Contraction
 *
 *  @param ar The archive
 *  @param cont The contraction
 */
template<typename Archive>
std::enable_if_t<madness::is_output_archive<Archive>::value, void>
serialize_contraction(Archive& ar, const libint2::Shell::Contraction& cont) {
    int l                           = cont.l;
    bool pure                       = cont.pure;
    libint2::svector<double> coeffs = cont.coeff;
    std::vector<double> std_coeffs(coeffs.begin(), coeffs.end());

    ar& l;
    ar& pure;
    ar& std_coeffs;
}

template<typename Archive>
std::enable_if_t<madness::is_input_archive<Archive>::value, void>
serialize_contraction(Archive& ar, libint2::Shell::Contraction& cont) {
    int l;
    bool pure;
    std::vector<double> std_coeffs;

    ar& l;
    ar& pure;
    ar& std_coeffs;

    libint2::svector<double> coeffs(std_coeffs.begin(), std_coeffs.end());

    cont = libint2::Shell::Contraction({l, pure, coeffs});
}

/** @brief Serialization of a Libint Shell
 *
 *  @param ar The archive
 *  @param shell The shell
 */
template<typename Archive>
std::enable_if_t<madness::is_output_archive<Archive>::value, void>
serialize_basis_shell(Archive& ar, const libint2::Shell& shell) {
    auto alphas                  = shell.alpha;
    auto conts                   = shell.contr;
    std::size_t nconts           = conts.size();
    std::array<double, 3> center = shell.O;
    std::vector<double> std_alphas(alphas.begin(), alphas.end());

    ar& std_alphas;
    ar& center;
    ar& nconts;

    for(auto cont : conts) { serialize_contraction(ar, cont); }
}

template<typename Archive>
std::enable_if_t<madness::is_input_archive<Archive>::value, void>
serialize_basis_shell(Archive& ar, libint2::Shell& shell) {
    std::vector<double> std_alphas;
    libint2::svector<libint2::Shell::Contraction> conts;
    std::size_t nconts = 0;
    std::array<double, 3> center;

    ar& std_alphas;
    ar& center;
    ar& nconts;

    for(int i = 0; i < nconts; ++i) {
        libint2::Shell::Contraction cont;
        serialize_contraction(ar, cont);
        conts.emplace_back(cont);
    }

    libint2::svector<double> alphas(std_alphas.begin(), std_alphas.end());
    shell = libint2::Shell(alphas, conts, center);
}

/** @brief Serialization of a Libint Basis Set
 *
 *  @param ar The archive
 *  @param set The basis set
 */
template<typename Archive>
std::enable_if_t<madness::is_output_archive<Archive>::value, void>
serialize_basis_set(Archive& ar, const libint2::BasisSet& set) {
    auto len = set.size();
    ar& len;

    for(auto shell : set) { serialize_basis_shell(ar, shell); }
}

template<typename Archive>
std::enable_if_t<madness::is_input_archive<Archive>::value, void>
serialize_basis_set(Archive& ar, libint2::BasisSet& set) {
    std::size_t len = 0;
    ar& len;

    for(auto i = 0; i < len; ++i) {
        libint2::Shell shell;
        serialize_basis_shell(ar, shell);
        set.emplace_back(shell);
    }
}

/** @brief Serialization of a std::vector of Libint Basis Sets
 *
 *  @param ar The archive
 *  @param sets The basis set vector
 */
template<typename Archive>
std::enable_if_t<madness::is_output_archive<Archive>::value, void>
serialize_basis_sets(Archive& ar, const std::vector<libint2::BasisSet>& sets) {
    std::size_t len = sets.size();
    ar& len;

    for(auto set : sets) { serialize_basis_set(ar, set); }
}

template<typename Archive>
std::enable_if_t<madness::is_input_archive<Archive>::value, void>
serialize_basis_sets(Archive& ar, std::vector<libint2::BasisSet>& sets) {
    std::size_t len = 0;
    ar& len;

    for(auto i = 0; i < len; ++i) {
        libint2::BasisSet set;
        serialize_basis_set(ar, set);
        sets.emplace_back(set);
    }
}

} // namespace nwx_TA