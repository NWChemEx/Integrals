#include "integrals/tamm_int_functor.hpp"
#include "schwarz_screening.hpp"

namespace integrals::libint::detail_ {

using size_type = std::size_t;

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
    schwarz_thresh{_schwarz_thresh} {
  if(schwarz_thresh > 0.0) {
    if constexpr (NBases == 3)
      Scr = schwarz_screening<op>(bases[1], bases[2]);
    else if constexpr (NBases == 4)
      Scr = schwarz_screening<op>(bases[0], bases[1]);
    libint2::initialize();
  }
};


template<libint2::Operator op, size_type NBases, typename element_type>
void TAMMIntFunctor<op, NBases, element_type>::operator()(
        const tamm::IndexVector& blockid, tamm::span<element_type> buff) {
    tamm_buf                   = buff;
    const size_type off_nopers = (nopers > 1) ? 1 : 0;
    iopers                     = (nopers > 1) ? blockid[0] : 0;

    for(size_type i = 0; i < NBases; ++i) {
        idx[i]    = blockid[i + off_nopers];
        maps[i]   = chemist::BasisSetMap{bases[i]};
        ao_off[i] = maps[i].atom_to_ao(atom_blocks[i][(idx[i])]).first;
    }

    fxn_call(0);
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
            x_libint = 0;
            if(schwarz_thresh > 0) {
                if(schwarz_estimate<NBases, element_type>(Scr, shells,
                                                          schwarz_thresh)) {
                    fill<true>(0, 0);
                } else {
                    fxn(shells);
                    if(fxn.engine.results()[iopers] == nullptr) {
                        fill<true>(0, 0);
                    } else {
                        fill<false>(0, 0);
                    }
                }
            } else {
                fxn(shells);
                if(fxn.engine.results()[iopers] == nullptr) {
                    fill<true>(0, 0);
                } else {
                    fill<false>(0, 0);
                }
            }
        } else {
            fxn_call(depth + 1);
        }
    }
}

template<libint2::Operator op, size_type NBases, typename element_type>
template<bool zero>
void TAMMIntFunctor<op, NBases, element_type>::fill(size_type x_tamm,
                                                    size_type depth) {
    const auto first_ao        = ao_ranges[depth].first;
    const auto second_ao       = ao_ranges[depth].second;
    const size_type off_nopers = (nopers > 1) ? 1 : 0;
    const auto tsize           = tAO[depth + off_nopers].tile_size(idx[depth]);

    x_tamm = x_tamm * tsize + (first_ao - ao_off[depth]);

    for(size_type mui = first_ao; mui < second_ao; ++mui) {
        if(depth == NBases - 1) {
            if (zero) {
                tamm_buf[x_tamm] = 0.0;
            } else {
                tamm_buf[x_tamm] = fxn.engine.results()[iopers][x_libint++];
            }
        } else {
            fill<zero>(x_tamm, depth + 1);
        }
        ++x_tamm;
    }
}

template class TAMMIntFunctor<libint2::Operator::overlap, 2, double>;
template class TAMMIntFunctor<libint2::Operator::kinetic, 2, double>;
template class TAMMIntFunctor<libint2::Operator::nuclear, 2, double>;
template class TAMMIntFunctor<libint2::Operator::coulomb, 2, double>;
template class TAMMIntFunctor<libint2::Operator::coulomb, 3, double>;
template class TAMMIntFunctor<libint2::Operator::coulomb, 4, double>;
template class TAMMIntFunctor<libint2::Operator::stg, 2, double>;
template class TAMMIntFunctor<libint2::Operator::stg, 3, double>;
template class TAMMIntFunctor<libint2::Operator::stg, 4, double>;
template class TAMMIntFunctor<libint2::Operator::yukawa, 2, double>;
template class TAMMIntFunctor<libint2::Operator::yukawa, 3, double>;
template class TAMMIntFunctor<libint2::Operator::yukawa, 4, double>;
template class TAMMIntFunctor<libint2::Operator::emultipole1, 2, double>;
template class TAMMIntFunctor<libint2::Operator::emultipole2, 2, double>;
template class TAMMIntFunctor<libint2::Operator::emultipole3, 2, double>;

} // namespace integrals::libint::detail_
