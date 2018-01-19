#pragma once
#include <memory>
#include <array>
#include <limits>
#include "IntegralsEx/nwx_libint/nwx_libint.hpp"
#include "IntegralsEx/FourCenterIntegral.hpp"
#include <libint2.hpp>


namespace nwx_libint {

inline size_t get_max(size_t v1,size_t v2, size_t v3, size_t v4){
    return std::max(std::max(std::max(v1,v2),v3),v4);
}

/** @brief Class with implementation of libint four-center integrals.
 *
 *  @tparam Op The libint2 operator to be used in the integral calculation.
 */
template<libint2::Operator Op>
class LibInt4C : public IntegralsEx::FourCenterIntegral {
private:
    std::array<libint2::BasisSet,4> bs_;
    libint2::Engine engine_;
public:
     LibInt4C(unsigned int deriv,
              const LibChemist::SetOfAtoms &atoms,
              const LibChemist::BasisSet &bs1,
              const LibChemist::BasisSet &bs2,
              const LibChemist::BasisSet &bs3,
              const LibChemist::BasisSet &bs4,
              double thresh = std::numeric_limits<double>::epsilon())
              :FourCenterIntegral(deriv,atoms,bs1,bs2,bs3,bs4,thresh)
     {
         libint2::initialize();

         bs_=std::array<libint2::BasisSet,4>({make_basis(bs1),
                                              make_basis(bs2),
                                              make_basis(bs3),
                                              make_basis(bs4)});
         const size_t max_prims=get_max(bs_[0].max_nprim(bs_[0]),
                                        bs_[1].max_nprim(bs_[1]),
                                        bs_[2].max_nprim(bs_[2]),
                                        bs_[3].max_nprim(bs_[3]));
         const size_t max_l=get_max(bs_[0].max_l(bs_[0]),
                                    bs_[1].max_l(bs_[1]),
                                    bs_[2].max_l(bs_[2]),
                                    bs_[3].max_l(bs_[3]));
         engine_=libint2::Engine(Op,max_prims,max_l,deriv,thresh);
     }

    ~LibInt4C()
    {
        libint2::finalize();
    }

     std::vector<const double*> calculate(ShellIndex shell1, 
                                          ShellIndex shell2,
                                          ShellIndex shell3, 
                                          ShellIndex shell4) override
     {
         const auto& buf_vec=engine_.results();
         engine_.compute(bs_[0][shell1],bs_[1][shell2],
                            bs_[2][shell3],bs_[3][shell4]);
         std::vector<const double*> rv;
         for (size_t i=0; i<n_components(); i++)
             rv.insert(rv.end(),{buf_vec[i]});
         return rv;
     }

};

using ERI=LibInt4C<libint2::Operator::coulomb>;
using Delta=LibInt4C<libint2::Operator::delta>;

}//End namespace
