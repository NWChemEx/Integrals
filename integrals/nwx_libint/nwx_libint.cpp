#include "integrals/nwx_libint/nwx_libint.hpp"

namespace nwx_libint {

template<typename T>
libint2::BasisSet _make_basis(const integrals::type::basis_set<T>& bs) {
    using coord_type = std::array<double, 3>;
    using Contraction = libint2::Shell::Contraction;
    libint2::BasisSet shells;

    for(const auto&& shelli : bs.shells()) {
        const auto nprims = shelli.n_unique_primitives();
        const auto first_prim = shelli.unique_primitive(0);
        const auto last_prim = shelli.unique_primitive(nprims - 1);
        const int l = shelli.l();
        const bool pure = shelli.pure();

        coord_type center = {shelli.x(), shelli.y(), shelli.z()};

        libint2::svector<double> alphas(&first_prim.exponent(), &last_prim.exponent() + 1);
        libint2::svector<double> coefs(&first_prim.coefficient(), &last_prim.coefficient() + 1);

        Contraction cont({l, pure, coefs});
        libint2::svector<Contraction> conts{cont};
        shells.push_back(libint2::Shell(alphas, conts, center));
    }

    return shells;
}

libint2::BasisSet make_basis(const integrals::type::basis_set<double>& bs) {
    return _make_basis(bs);
}

libint2::BasisSet make_basis(const integrals::type::basis_set<float>& bs) {
    return _make_basis(bs);
}

} // namespace nwx_libint
