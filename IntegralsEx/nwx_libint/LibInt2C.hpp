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

     std::vector<const double*> calculate(ShellIndex shell1, ShellIndex shell2) override
     {
         const auto& buf_vec=engine_.results();
         engine_.compute(bs_[0][shell1],bs_[1][shell2]);

         std::vector<const double*> rv(buf_vec.begin(),buf_vec.end());
         return rv;
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
		qs.push_back({ai.property(LibChemist::AtomProperty::charge)==0? ai.property(LibChemist::AtomProperty::Z): ai.property(LibChemist::AtomProperty::charge) ,{ai.coord[0],ai.coord[1],ai.coord[2]}});
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

struct EDipole: public LibInt2C<libint2::Operator::emultipole1> {
	using base_type=LibInt2C<libint2::Operator::emultipole1>;

    EDipole(unsigned int deriv,
           const LibChemist::SetOfAtoms &atoms,
           const LibChemist::BasisSet &bs1,
           const LibChemist::BasisSet &bs2) : base_type(deriv, atoms, bs1, bs2){}    

    unsigned int n_components(void) const override {return 4;}
};

struct EQuadrupole: public LibInt2C<libint2::Operator::emultipole2> {
	using base_type=LibInt2C<libint2::Operator::emultipole2>;

    EQuadrupole(unsigned int deriv,
           const LibChemist::SetOfAtoms &atoms,
           const LibChemist::BasisSet &bs1,
           const LibChemist::BasisSet &bs2) : base_type(deriv, atoms, bs1, bs2){}    
   
    unsigned int n_components(void) const override {return 10;}
};

struct EOctopole: public LibInt2C<libint2::Operator::emultipole3> {
	using base_type=LibInt2C<libint2::Operator::emultipole3>;

    EOctopole(unsigned int deriv,
           const LibChemist::SetOfAtoms &atoms,
           const LibChemist::BasisSet &bs1,
           const LibChemist::BasisSet &bs2) : base_type(deriv, atoms, bs1, bs2){}    
   
    unsigned int n_components(void) const override {return 20;}
};

}//End namespace
