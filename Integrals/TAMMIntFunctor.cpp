#include "Integrals/TAMMIntFunctor.hpp"

namespace Integrals::Libint::detail_ {

using size_type = std::size_t;
using matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template<libint2::Operator op>
matrix schwarz_screening(const LibChemist::AOBasisSet& bs1,
                         const LibChemist::AOBasisSet& bs2) {

    const auto nsh1 = bs1.nshells();
    const auto nsh2 = bs2.nshells();
    const bool bs1_equiv_bs2 = (bs1 == bs2);
    matrix rv = matrix::Zero(nsh1,nsh2);
    LibIntFunctor<4> fxn;

    fxn.bs[0] = nwx_libint::make_basis(bs1);
    fxn.bs[1] = nwx_libint::make_basis(bs2);
    fxn.bs[2] = fxn.bs[0];
    fxn.bs[3] = fxn.bs[1];
    auto max_prims = std::max(fxn.bs[0].max_nprim(fxn.bs[0]), fxn.bs[1].max_nprim(fxn.bs[1]));
    auto max_l = std::max(fxn.bs[0].max_l(fxn.bs[0]), fxn.bs[1].max_l(fxn.bs[1]));

    fxn.engine = make_engine<op, 4>(LibChemist::Molecule{}, max_prims, max_l, 0.0, 0);

    for (size_type s1 = 0, s12 = 0; s1 != nsh1; ++s1) {
        size_type n1 = bs1[s1].size();  // number of basis functions in this shell

        size_type s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
        for (size_type s2 = 0; s2 <= s2_max; ++s2, ++s12) {
            size_type n2 = bs2[s2].size();
            size_type n12 = n1 * n2;
            std::array<size_type,4> shells{s1,s2,s1,s2};
            auto buf = fxn(shells);

            Eigen::Map<const matrix> buf_mat(buf[0], n12, n12);
            auto norm2 = buf_mat.lpNorm<Eigen::Infinity>();
            rv(s1, s2) = std::sqrt(norm2);
            if (bs1_equiv_bs2) rv(s2, s1) = rv(s1, s2);
        }
    }
    return rv;
}

template<libint2::Operator op, size_type NBases, typename element_type>
TAMMIntFunctor<op,NBases,element_type>::TAMMIntFunctor(const tiled_AO &tAO,
               const std::array<std::vector<size_type>, NBases> &atom_blocks,
               const basis_array_type &bases,
               fxn_type &&fxn) : tAO{std::move(tAO)}, atom_blocks{std::move(atom_blocks)},
                                  bases{std::move(bases)}, fxn{std::move(fxn)} {};

template<>
TAMMIntFunctor<libint2::Operator::coulomb,3,double>::TAMMIntFunctor(const tiled_AO& tAO,
                                                                    const std::array<std::vector<size_type>, 3>& atom_blocks,
                                                                    const basis_array_type& bases,
                                                                    fxn_type&& fxn) : tAO{std::move(tAO)}, atom_blocks{std::move(atom_blocks)},
                                                                                      bases{std::move(bases)}, fxn{std::move(fxn)} {
    screen = true;
    Scr = schwarz_screening<libint2::Operator::coulomb>(bases[1],bases[2]);
    libint2::initialize();
}

template<>
TAMMIntFunctor<libint2::Operator::coulomb,4,double>::TAMMIntFunctor(const tiled_AO& tAO,
                                                                    const std::array<std::vector<size_type>, 4>& atom_blocks,
                                                                    const basis_array_type& bases,
                                                                    fxn_type&& fxn) : tAO{std::move(tAO)}, atom_blocks{std::move(atom_blocks)},
                                                                                      bases{std::move(bases)}, fxn{std::move(fxn)} {
    screen = true;
    Scr = schwarz_screening<libint2::Operator::coulomb>(bases[0],bases[1]);
    libint2::initialize();
}

template<libint2::Operator op, size_type NBases, typename element_type>
void TAMMIntFunctor<op,NBases,element_type>::fill(size_type x_tamm,
                                                  size_type depth) {

    const auto first_ao = ao_ranges[depth].first;
    const auto second_ao = ao_ranges[depth].second;
    const size_type off_nopers = (nopers > 1) ? 1 : 0;
    const auto tsize = tAO[depth+off_nopers].tile_size(idx[depth]);

    x_tamm = x_tamm*tsize + (first_ao - ao_off[depth]);

    for (size_type mui = first_ao; mui < second_ao; ++mui) {
        if (depth == NBases - 1) {
            tamm_buf[x_tamm] = libint_buf[x_libint++];
        } else {
            fill(x_tamm,depth+1);
        }
        x_tamm++;
    }
}

template<libint2::Operator op, size_type NBases, typename element_type>
void TAMMIntFunctor<op,NBases,element_type>::fxn_call(size_type depth) {
    auto first_shell = maps[depth].atom_to_shell((atom_blocks[depth])[(idx[depth])]).first;
    auto second_shell = maps[depth].atom_to_shell(((atom_blocks[depth])[(idx[depth]) + 1]) - 1).second;

    for (size_type si = first_shell; si < second_shell; ++si) {
        shells[depth] = si;
        ao_ranges[depth] = maps[depth].shell_to_ao(si);
        if (depth == NBases - 1) {
            size_type nbfs = 1ul;
            for(size_type i=0; i<NBases; ++i)
                nbfs *= bases[i][shells[i]].size();

            const element_type sch_thresh = 1.0E-8;

            if (screen) {
                element_type estimate;
                if (NBases == 3) estimate = Scr(shells[1],shells[2]);
                if (NBases == 4) estimate = Scr(shells[0],shells[1])*Scr(shells[2],shells[3]);
                if (estimate < sch_thresh) {
                    libint_buf = std::vector<double>(nbfs, 0.0);
                } else {
                    auto buffer = fxn(shells);
                    if (buffer[iopers] == nullptr) {
                        libint_buf = std::vector<double>(nbfs, 0.0);
                    } else {
                        libint_buf = std::vector<double>(buffer[iopers], buffer[iopers] + nbfs);
                    }
                }
            } else {
                auto buffer = fxn(shells);
                if (buffer[iopers] == nullptr) {
                    libint_buf = std::vector<double>(nbfs, 0.0);
                } else {
                    libint_buf = std::vector<double>(buffer[iopers], buffer[iopers] + nbfs);
                }
            }

            x_libint = 0;
            fill(0,0);
        } else {
            fxn_call(depth + 1);
        }
    }
}

template<libint2::Operator op, size_type NBases, typename element_type>
void TAMMIntFunctor<op,NBases,element_type>::operator()(const tamm::IndexVector& blockid, tamm::span<element_type> buff) {

    tamm_buf = buff;
    const size_type off_nopers = (nopers > 1) ? 1 : 0;
    iopers = (nopers > 1) ? blockid[0] : 0;

    for (size_type i = 0; i < NBases; ++i) {
        idx[i] = blockid[i+off_nopers];
        maps[i] = LibChemist::BasisSetMap{bases[i]};
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

}