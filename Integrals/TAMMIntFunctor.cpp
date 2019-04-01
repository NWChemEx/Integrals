#include "Integrals/TAMMIntFunctor.hpp"

namespace Integrals::Libint::detail_ {

template<libint2::Operator op, std::size_t NBases, typename element_type>
TAMMIntFunctor<op,NBases,element_type>::TAMMIntFunctor(const tiled_AO &tAO,
               const std::array<std::vector<size_type>, NBases> &atom_blocks,
               const basis_array_type &bases,
               fxn_type &&fxn) : tAO{std::move(tAO)}, atom_blocks{std::move(atom_blocks)},
                                  bases{std::move(bases)}, fxn{std::move(fxn)} {};

template<>
TAMMIntFunctor<libint2::Operator::coulomb,4,double>::TAMMIntFunctor(const tiled_AO& tAO,
                                                                    const std::array<std::vector<size_type>, 4>& atom_blocks,
                                                                    const basis_array_type& bases,
                                                                    fxn_type&& fxn) : tAO{std::move(tAO)}, atom_blocks{std::move(atom_blocks)},
                                                                                      bases{std::move(bases)}, fxn{std::move(fxn)} {

    screen = true;
    const auto nsh1 = bases[0].nshells();
    const auto nsh2 = bases[2].nshells();
    const bool bs1_equiv_bs2 = (bases[0] == bases[2]);
    Scr = matrix::Zero(nsh1,nsh2);
    fxn.engine.set_precision(0);

    for (size_type s1 = 0, s12 = 0; s1 != nsh1; ++s1) {
        size_type n1 = bases[0][s1].size();  // number of basis functions in this shell

        size_type s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
        for (size_type s2 = 0; s2 <= s2_max; ++s2, ++s12) {
            size_type n2 = bases[2][s2].size();
            size_type n12 = n1 * n2;
            std::array<size_type,4> shells{s1,s2,s1,s2};
            auto buf = fxn(shells);

            Eigen::Map<const matrix> buf_mat(buf[0], n12, n12);
            auto norm2 = buf_mat.lpNorm<Eigen::Infinity>();
            Scr(s1, s2) = std::sqrt(norm2);
            if (bs1_equiv_bs2) Scr(s2, s1) = Scr(s1, s2);
        }
    }
}

template<libint2::Operator op, std::size_t NBases, typename element_type>
void TAMMIntFunctor<op,NBases,element_type>::fill(size_type x_tamm,
                                                  size_type x_libint,
                                                  size_type depth,
                                                  size_type off_nopers,
                                                  std::array<size_type, NBases>& idx,
                                                  range_array& ao_ranges,
                                                  std::array<size_type, NBases>& ao_off,
                                                  tamm::span<element_type> tamm_buf,
                                                  std::vector<double>& libint_buf) {

    const auto first_ao = ao_ranges[depth].first;
    const auto second_ao = ao_ranges[depth].second;
    const auto tsize = tAO[depth+off_nopers].tile_size(idx[depth]);

    x_tamm = x_tamm*tsize + (first_ao - ao_off[depth]);

    for (size_type mui = first_ao; mui < second_ao; ++mui) {
        if (depth == NBases - 1) {
            tamm_buf[x_tamm] = libint_buf[x_libint++];
        } else {
            fill(x_tamm,x_libint,depth+1,off_nopers,idx,ao_ranges,ao_off,tamm_buf,libint_buf);
        }
        x_tamm++;
    }
}

template<libint2::Operator op, std::size_t NBases, typename element_type>
void TAMMIntFunctor<op,NBases,element_type>::fxn_call(typename fxn_type::shell_index& shells,
                                                      range_array& ao_ranges,
                                                      size_type depth,
                                                      std::array<LibChemist::BasisSetMap, NBases>& maps,
                                                      std::array<size_type, NBases>& idx,
                                                      size_type off_nopers,
                                                      size_type iopers,
                                                      std::array<size_type, NBases>& ao_off,
                                                      tamm::span<element_type> tamm_buf) {
    auto first_shell = maps[depth].atom_to_shell((atom_blocks[depth])[(idx[depth])]).first;
    auto second_shell = maps[depth].atom_to_shell(((atom_blocks[depth])[(idx[depth]) + 1]) - 1).second;

    for (size_type si = first_shell; si < second_shell; ++si) {
        shells[depth] = si;
        ao_ranges[depth] = maps[depth].shell_to_ao(si);
        if (depth == NBases - 1) {
            std::size_t nbfs = 1ul;
            for(std::size_t i=0; i<NBases; ++i)
                nbfs *= bases[i][shells[i]].size();

            std::vector<double> libint_buf;

            const element_type sch_thresh = 1.0E-5;

            if (screen) {
                auto estimate = Scr(shells[0],shells[1])*Scr(shells[2],shells[3]);
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

            size_type x_libint = 0;
            fill(0,0,0,off_nopers,idx,ao_ranges,ao_off,tamm_buf,libint_buf);
        } else {
            fxn_call(shells, ao_ranges, depth + 1, maps, idx, off_nopers, iopers, ao_off, tamm_buf);
        }
    }
}

template<libint2::Operator op, std::size_t NBases, typename element_type>
void TAMMIntFunctor<op,NBases,element_type>::operator()(const tamm::IndexVector& blockid, tamm::span<element_type> tamm_buf) {

    std::array<size_type, NBases> idx;
    std::array<LibChemist::BasisSetMap, NBases> maps;
    std::array<size_type, NBases> ao_off;
    const size_type off_nopers = (nopers > 1) ? 1 : 0;
    for (size_type i = 0; i < NBases; ++i) {
        idx[i] = blockid[i+off_nopers];
        maps[i] = LibChemist::BasisSetMap{bases[i]};
        ao_off[i] = maps[i].atom_to_ao(atom_blocks[i][(idx[i])]).first;
    }
    const size_type iopers = (nopers > 1) ? blockid[0] : 0;

    typename fxn_type::shell_index shells;
    using range_array = std::array<typename LibChemist::BasisSetMap::range,NBases>;
    range_array ao_ranges;

    fxn_call(shells, ao_ranges, 0, maps, idx, off_nopers, iopers, ao_off, tamm_buf);
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