#include "integrals/tamm_int_functor.hpp"

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

    for(size_type s1 = 0, s12 = 0; s1 != nsh1; ++s1) {
        size_type n1 =
          bs1[s1].size(); // number of basis functions in this shell

        size_type s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
        for(size_type s2 = 0; s2 <= s2_max; ++s2, ++s12) {
            size_type n2  = bs2[s2].size();
            size_type n12 = n1 * n2;
            std::array<size_type, 4> shells{s1, s2, s1, s2};
            auto buf = fxn(shells);

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

template<libint2::Operator op, size_type NBases, typename element_type>
TAMMIntFunctor<op, NBases, element_type>::TAMMIntFunctor(
  const tiled_AO& _tAO,
  const std::array<std::vector<size_type>, NBases>& _atom_blocks,
  const basis_array_type& _bases, fxn_type&& _fxn,
  const element_type _schwarz_thresh) :
  tAO{std::move(_tAO)},
  atom_blocks{std::move(_atom_blocks)},
  bases{std::move(_bases)},
  fxn{std::move(_fxn)},
  schwarz_thresh{_schwarz_thresh} {};

template<>
TAMMIntFunctor<libint2::Operator::coulomb, 3, double>::TAMMIntFunctor(
  const tiled_AO& _tAO,
  const std::array<std::vector<size_type>, 3>& _atom_blocks,
  const basis_array_type& _bases, fxn_type&& _fxn,
  const double _schwarz_thresh) :
  tAO{std::move(_tAO)},
  atom_blocks{std::move(_atom_blocks)},
  bases{std::move(_bases)},
  fxn{std::move(_fxn)},
  schwarz_thresh{_schwarz_thresh} {
    if(schwarz_thresh > 0.0) {
        Scr = schwarz_screening<libint2::Operator::coulomb>(bases[1], bases[2]);
        libint2::initialize();
    }
}

template<>
TAMMIntFunctor<libint2::Operator::coulomb, 4, double>::TAMMIntFunctor(
  const tiled_AO& _tAO,
  const std::array<std::vector<size_type>, 4>& _atom_blocks,
  const basis_array_type& _bases, fxn_type&& _fxn,
  const double _schwarz_thresh) :
  tAO{std::move(_tAO)},
  atom_blocks{std::move(_atom_blocks)},
  bases{std::move(_bases)},
  fxn{std::move(_fxn)},
  schwarz_thresh{_schwarz_thresh} {
    if(schwarz_thresh > 0.0) {
        Scr = schwarz_screening<libint2::Operator::coulomb>(bases[0], bases[1]);
        libint2::initialize();
    }
}

template<libint2::Operator op, size_type NBases, typename element_type>
void TAMMIntFunctor<op, NBases, element_type>::fill(size_type x_tamm,
                                                    size_type depth) {
    const auto first_ao        = ao_ranges[depth].first;
    const auto second_ao       = ao_ranges[depth].second;
    const size_type off_nopers = (nopers > 1) ? 1 : 0;
    const auto tsize           = tAO[depth + off_nopers].tile_size(idx[depth]);

    x_tamm = x_tamm * tsize + (first_ao - ao_off[depth]);

    for(size_type mui = first_ao; mui < second_ao; ++mui) {
        if(depth == NBases - 1) {
            tamm_buf[x_tamm] = libint_buf[x_libint++];
        } else {
            fill(x_tamm, depth + 1);
        }
        x_tamm++;
    }
}

template<libint2::Operator op, size_type NBases, typename element_type>
void TAMMIntFunctor<op, NBases, element_type>::fxn_call(size_type depth) {
    auto first_shell =
      maps[depth].atom_to_shell((atom_blocks[depth])[(idx[depth])]).first;
    auto second_shell =
      maps[depth]
        .atom_to_shell(((atom_blocks[depth])[(idx[depth]) + 1]) - 1)
        .second;

    for(size_type si = first_shell; si < second_shell; ++si) {
        shells[depth]    = si;
        ao_ranges[depth] = maps[depth].shell_to_ao(si);
        if(depth == NBases - 1) {
            size_type nbfs = 1ul;
            for(size_type i = 0; i < NBases; ++i)
                nbfs *= bases[i][shells[i]].size();

            if(schwarz_thresh > 0) {
                if(schwarz_estimate<NBases, element_type>(Scr, shells,
                                                          schwarz_thresh)) {
                    libint_buf = std::vector<double>(nbfs, 0.0);
                } else {
                    auto buffer = fxn(shells);
                    if(buffer[iopers] == nullptr) {
                        libint_buf = std::vector<double>(nbfs, 0.0);
                    } else {
                        libint_buf = std::vector<double>(buffer[iopers],
                                                         buffer[iopers] + nbfs);
                    }
                }
            } else {
                auto buffer = fxn(shells);
                if(buffer[iopers] == nullptr) {
                    libint_buf = std::vector<double>(nbfs, 0.0);
                } else {
                    libint_buf = std::vector<double>(buffer[iopers],
                                                     buffer[iopers] + nbfs);
                }
            }

            x_libint = 0;
            fill(0, 0);
        } else {
            fxn_call(depth + 1);
        }
    }
}

template<libint2::Operator op, size_type NBases, typename element_type>
void TAMMIntFunctor<op, NBases, element_type>::operator()(
  const tamm::IndexVector& blockid, tamm::span<element_type> buff) {
    tamm_buf                   = buff;
    const size_type off_nopers = (nopers > 1) ? 1 : 0;
    iopers                     = (nopers > 1) ? blockid[0] : 0;

    for(size_type i = 0; i < NBases; ++i) {
        idx[i]    = blockid[i + off_nopers];
        maps[i]   = libchemist::BasisSetMap{bases[i]};
        ao_off[i] = maps[i].atom_to_ao(atom_blocks[i][(idx[i])]).first;
    }

    fxn_call(0);
}

template class TAMMIntFunctor<libint2::Operator::overlap, 2, double>;
template class TAMMIntFunctor<libint2::Operator::kinetic, 2, double>;
template class TAMMIntFunctor<libint2::Operator::nuclear, 2, double>;
template class TAMMIntFunctor<libint2::Operator::coulomb, 2, double>;
template class TAMMIntFunctor<libint2::Operator::coulomb, 3, double>;
template class TAMMIntFunctor<libint2::Operator::coulomb, 4, double>;
template class TAMMIntFunctor<libint2::Operator::emultipole1, 2, double>;
template class TAMMIntFunctor<libint2::Operator::emultipole2, 2, double>;
template class TAMMIntFunctor<libint2::Operator::emultipole3, 2, double>;

} // namespace integrals::libint::detail_
