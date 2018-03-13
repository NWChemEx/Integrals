#pragma once
#include <memory>
#include <array>
#include "Integrals/nwx_libint/nwx_libint.hpp"
#include "Integrals/ThreeCenterIntegral.hpp"
#include <libint2.hpp>


namespace nwx_libint {

/** @brief Class with implementation of libint three-center integrals.
 *
 *  @tparam Op The libint2 operator to be used in the integral calculation.
 */
template<libint2::Operator Op>
class LibInt3C : public Integrals::ThreeCenterIntegral {
protected:
    std::array<libint2::BasisSet,3> bs_;
    libint2::Engine engine_;
public:
     LibInt3C(unsigned int deriv,
              const LibChemist::Molecule &molecule,
              const LibChemist::BasisSet &bs1,
              const LibChemist::BasisSet &bs2,
              const LibChemist::BasisSet &bs3):ThreeCenterIntegral(deriv,molecule,bs1,bs2,bs3)
     {
         libint2::initialize();
         bs_=std::array<libint2::BasisSet,3>({make_basis(bs1),
                                              make_basis(bs2),
                                              make_basis(bs3)});
         const size_t max_prims=std::max(
                      std::max(bs_[0].max_nprim(bs_[0]),bs_[1].max_nprim(bs_[1])),
                               bs_[2].max_nprim(bs_[2]));
         const size_t max_l=std::max(
                     std::max(bs_[0].max_l(bs_[0]),bs_[1].max_l(bs_[1])),
                               bs_[2].max_l(bs_[2]));

         engine_=libint2::Engine(Op,max_prims,max_l,deriv);
         engine_.set_braket(libint2::BraKet::xs_xx);
     }

    ~LibInt3C()
    {
        libint2::finalize();
    }

     std::vector<const double*> calculate(ShellIndex shell1,
                                          ShellIndex shell2,
                                          ShellIndex shell3) override
     {
         const auto& buf_vec=engine_.results();
         engine_.compute(bs_[0][shell1],bs_[1][shell2],bs_[2][shell3]);

         std::vector<const double*> rv(buf_vec.begin(),buf_vec.end());
         return rv;
     }

};

using DF3C2E=LibInt3C<libint2::Operator::coulomb>;

}//End namespace
