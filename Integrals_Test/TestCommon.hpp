#pragma once
#include <SDE/NWXDefaults.hpp>
#include <SDE/MoleculeFileParser.hpp>
#include <catch/catch.hpp>
#include <fstream>

using namespace LibChemist;
using namespace SDE;

struct ptr_wrapper{
    const double* ptr_;
    size_t n_;
    const double* begin()const{return ptr_;}
    const double* end()const{return ptr_+n_;}
};

    const double eps = 1000*std::numeric_limits<double>::epsilon();
    const double marg = 10*std::numeric_limits<double>::epsilon();

inline void compare_integrals(const ptr_wrapper& calc, 
                              const std::vector<double>& corr)
{
    size_t i = 0;
    for (const double& x : calc) {
        INFO("epsilon = " << eps); 
        REQUIRE(x == Approx(corr[i++]).epsilon(eps).margin(marg));
    }
}

Molecule make_molecule()
{
    std::map<size_t,std::vector<BasisShell>> bs;
    bs[1].push_back(
	  BasisShell{ShellType::CartesianGaussian,0,1,
	      {3.42525091,0.62391373,0.16885540},
	      {0.15432897,0.53532814,0.44463454}});
    bs[8].push_back(
	  BasisShell{ShellType::CartesianGaussian,0,1,
	      {130.7093200,23.8088610,6.4436083},
	      {0.15432897,0.53532814,0.44463454}});
    bs[8].push_back(
	  BasisShell{ShellType::CartesianGaussian,-1,2,
	      {5.0331513,1.1695961,0.3803890},
	      {-0.09996723,0.39951283,0.70011547,
	        0.15591627,0.60768372,0.39195739}});

    auto crt = default_runtime();
    auto xyzfile = std::ifstream("water.xyz");
    auto water = parse_molecule_file(xyzfile, XYZParser(), crt);

    auto charge = Atom::Property::charge;
    for(auto& atomi : water.atoms) {
	const size_t Z = std::lround(atomi.properties.at(charge));
	atomi.bases["sto-3g_cart"] = bs.at(Z);
    }

    return water;
}

