#pragma once
#include <memory>
#include <array>
#include "IntegralsEx/nwx_libint/nwx_libint.hpp"
#include "IntegralsEx/TwoCenterIntegral.hpp"
#include <libint2.hpp>


namespace nwx_libint {
/** @brief Class with implementation of libint two-center integrals.
 *
 *  @tparam Op The libint2 operator to be used in the integral calculation.
 */
template<libint2::Operator Op>
class LibInt2C : public IntegralsEx::TwoCenterIntegral {
protected:
    std::array<libint2::BasisSet,2> bs_;
    libint2::Engine engine_;
public:
     LibInt2C(unsigned int deriv,
              const LibChemist::SetOfAtoms &atoms,
              const LibChemist::BasisSet &bs1,
              const LibChemist::BasisSet &bs2):TwoCenterIntegral(deriv,atoms,bs1,bs2)
     {
       libint2::initialize();

       bs_=std::array<libint2::BasisSet,2>({make_basis(bs1),make_basis(bs2)});
       const size_t max_prims=std::max(bs_[0].max_nprim(bs_[0]),
                                       bs_[1].max_nprim(bs_[1]));
       const size_t max_l=std::max(bs_[0].max_l(bs_[0]),
                                   bs_[1].max_l(bs_[1]));

       engine_=libint2::Engine(Op,max_prims,max_l,deriv);
     }

     const double* calculate_(uint64_t shell1, uint64_t shell2)
     {
         const auto& buf_vec=engine_.results();
         engine_.compute(bs_[0][shell1],bs_[1][shell2]);
         return buf_vec[0];
     }
  
    ~LibInt2C()
     {
        libint2::finalize();
     }

};

using Kinetic=LibInt2C<libint2::Operator::kinetic>;
using Overlap=LibInt2C<libint2::Operator::overlap>;

struct NuclearElectron: public LibInt2C<libint2::Operator::nuclear> {
	using base_type=LibInt2C<libint2::Operator::nuclear>;

	NuclearElectron(unsigned int deriv,
                    const LibChemist::SetOfAtoms &atoms,
                    const LibChemist::BasisSet &bs1,
                    const LibChemist::BasisSet &bs2) : base_type(deriv, atoms, bs1, bs2)
    {
        std::vector<std::pair<double,std::array<double,3>>> qs;
        for(const auto& ai: atoms)
            qs.push_back({ai.charge==0? ai.Z: ai.charge ,{ai.coord[0],ai.coord[1],ai.coord[2]}});
        engine_.set_params(qs);

    }
};

struct Metric: public LibInt2C<libint2::Operator::coulomb> {
	using base_type=LibInt2C<libint2::Operator::coulomb>;

    Metric(unsigned int deriv,
           const LibChemist::SetOfAtoms &atoms,
           const LibChemist::BasisSet &bs1,
           const LibChemist::BasisSet &bs2) : base_type(deriv, atoms, bs1, bs2)
    {
         engine_.set_braket(libint2::BraKet::xs_xs);
    }
};

}//End namespace
