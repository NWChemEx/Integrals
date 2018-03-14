#pragma once
#include <SDE/NWXDefaults.hpp>
#include <SDE/BasisSetFileParser.hpp>
#include <SDE/MoleculeFileParser.hpp>
#include <catch/catch.hpp>
#include <string>
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

Molecule apply_basis_file(const std::string& key, std::istream& is, 
                                Molecule mol, ChemistryRuntime& crt) {
    auto charge = LibChemist::Atom::Property::charge;
    auto bs = parse_basis_set_file(is, G94(), crt);
    for(auto& atomi : mol.atoms) {
        const size_t Z = std::lround(atomi.properties.at(charge));
        atomi.bases[key] = bs.at(Z);
    }
    return mol;
}

Molecule make_molecule()
{
    auto crt = default_runtime();
    auto xyzfile = std::ifstream("water.xyz");
    auto water = parse_molecule_file(xyzfile, XYZParser(), crt);
    auto fs = std::ifstream("sto-3g.gbs");    
    auto water_w_basis = apply_basis_file("sto-3gfile", fs, water, crt);

    return water_w_basis;
}

