namespace integrals::libint::detail_ {

using size_type = std::size_t;
using matrix_type =
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template<libint2::Operator op>
matrix_type schwarz_screening(const libchemist::AOBasisSet& bs1,
                              const libchemist::AOBasisSet& bs2) {
    const auto nsh1          = bs1.nshells();
    const auto nsh2          = bs2.nshells();
    const bool bs1_equiv_bs2 = (bs1 == bs2);
    matrix_type rv           = matrix_type::Zero(nsh1, nsh2);
    LibintFunctor<4> fxn;

    fxn.bs[0] = nwx_libint::make_basis(bs1);
    fxn.bs[1] = nwx_libint::make_basis(bs2);
    fxn.bs[2] = fxn.bs[0];
    fxn.bs[3] = fxn.bs[1];
    auto max_prims =
            std::max(fxn.bs[0].max_nprim(fxn.bs[0]), fxn.bs[1].max_nprim(fxn.bs[1]));
    auto max_l =
            std::max(fxn.bs[0].max_l(fxn.bs[0]), fxn.bs[1].max_l(fxn.bs[1]));

    fxn.engine =
            make_engine<op, 4>(libchemist::Molecule{}, max_prims, max_l, 0.0, 0);

    const auto& buf = fxn.engine.results();

    for(size_type s1 = 0, s12 = 0; s1 != nsh1; ++s1) {
        size_type n1 =
                bs1[s1].size(); // number of basis functions in this shell

        size_type s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
        for(size_type s2 = 0; s2 <= s2_max; ++s2, ++s12) {
            size_type n2  = bs2[s2].size();
            size_type n12 = n1 * n2;
            std::array<size_type, 4> shells{s1, s2, s1, s2};
            fxn(shells);

            Eigen::Map<const matrix_type> buf_mat(buf[0], n12, n12);
            auto norm2 = buf_mat.lpNorm<Eigen::Infinity>();
            rv(s1, s2) = std::sqrt(norm2);
            if(bs1_equiv_bs2) rv(s2, s1) = rv(s1, s2);
        }
    }
    return rv;
}

template<size_type NBases, typename element_type>
bool schwarz_estimate(const matrix_type& mat,
                      const std::array<size_type, NBases>& shells,
                      const element_type threshold) {
    if constexpr(NBases == 3) {
        return threshold > mat(shells[1], shells[2]);
    } else if constexpr(NBases == 4) {
        return threshold >
               mat(shells[0], shells[1]) * mat(shells[2], shells[3]);
    } else {
        return false;
    }
}