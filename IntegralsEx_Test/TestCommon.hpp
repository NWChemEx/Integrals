#pragma once
#include <LibChemist/SetOfAtomsParser.hpp>
#include <LibChemist/BasisSetParser.hpp>

using namespace LibChemist;

struct ptr_wrapper{
    const double* ptr_;
    size_t n_;
    const double* begin()const{return ptr_;}
    const double* end()const{return ptr_+n_;}
};

inline SetOfAtoms make_atoms()
{
const double angstrom_to_bohr = 1 / 0.52917721092;

const double Oy=-0.07579*angstrom_to_bohr;
const double Hx=0.86681*angstrom_to_bohr;
const double Hy=0.60144*angstrom_to_bohr;

SetOfAtoms atoms;
atoms.insert(create_atom({0.0000,Oy,0.0000},8));
atoms.insert(create_atom({Hx,Hy,0.000000},1));
atoms.insert(create_atom({-Hx,Hy,0.00000},1));

// Build STO-3G basis
std::map<size_t,std::vector<BasisShell>> bs;
bs[1].push_back(
      BasisShell(ShellType::CartesianGaussian,0,1,
          std::vector<double>({3.42525091,0.62391373,0.16885540}),
          std::vector<double>({0.15432897,0.53532814,0.44463454})));
bs[8].push_back(
      BasisShell(ShellType::CartesianGaussian,0,1,
          std::vector<double>({130.7093200,23.8088610,6.4436083}),
          std::vector<double>({0.15432897,0.53532814,0.44463454})));
bs[8].push_back(
      BasisShell(ShellType::CartesianGaussian,-1,2,
          std::vector<double>({5.0331513,1.1695961,0.3803890}),
          std::vector<double>({-0.09996723,0.39951283,0.70011547,
                                0.15591627,0.60768372,0.39195739})));

auto atoms_with_basis=apply_basis_set("PRIMARY",bs,atoms);
return atoms_with_basis;
}
