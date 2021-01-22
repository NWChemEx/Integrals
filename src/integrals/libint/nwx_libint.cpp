#include "nwx_libint.hpp"
#include <array>

namespace nwx_libint {

template<typename T>
LI_basis _make_basis(const NWX_basis<T>& bs) {
    using coord_type  = std::array<double, 3>;
    using Contraction = libint2::Shell::Contraction;
    LI_basis shells;

    for(const auto&& shelli : bs.shells()) {
        const auto nprims     = shelli.n_unique_primitives();
        const auto first_prim = shelli.unique_primitive(0);
        const auto last_prim  = shelli.unique_primitive(nprims - 1);
        const int l           = shelli.l();
        const bool pure       = shelli.pure();

        coord_type center = {shelli.x(), shelli.y(), shelli.z()};

        libint2::svector<double> alphas(&first_prim.exponent(),
                                        &last_prim.exponent() + 1);
        libint2::svector<double> coefs(&first_prim.coefficient(),
                                       &last_prim.coefficient() + 1);

        Contraction cont({l, pure, coefs});
        libint2::svector<Contraction> conts{cont};
        shells.push_back(libint2::Shell(alphas, conts, center));
    }

    return shells;
}

LI_basis make_basis(const NWX_basis<double>& bs) { return _make_basis(bs); }

LI_basis make_basis(const NWX_basis<float>& bs) { return _make_basis(bs); }

template<typename T>
std::vector<LI_basis> _make_basis_sets(const std::vector<NWX_basis<T>>& sets) {
    std::vector<LI_basis> LI_basis_sets{};

    for(auto set : sets) { LI_basis_sets.push_back(make_basis(set)); }

    return LI_basis_sets;
}

std::vector<LI_basis> make_basis_sets(
  const std::vector<NWX_basis<double>>& sets) {
    return _make_basis_sets(sets);
}

std::vector<LI_basis> make_basis_sets(
  const std::vector<NWX_basis<float>>& sets) {
    return _make_basis_sets(sets);
}

size sets_max_nprims(const std::vector<LI_basis>& sets) {
    size max_nprims = 0;

    for(auto set : sets) {
        max_nprims = std::max(max_nprims, libint2::max_nprim(set));
    }

    return max_nprims;
}

int sets_max_l(const std::vector<LI_basis>& sets) {
    int max_l = 0;

    for(auto set : sets) { max_l = std::max(max_l, libint2::max_l(set)); }

    return max_l;
}

std::vector<size> aos2shells(libint2::BasisSet basis_set, size lower,
                             size upper) {
    std::vector<size> return_vec;

    for(auto ishell = 0, offset = 0; ishell < basis_set.size(); ++ishell) {
        if(offset >= upper) break;
        if(offset >= lower) return_vec.push_back(ishell);
        offset += basis_set[ishell].size();
    }

    return return_vec;
}

} // namespace nwx_libint