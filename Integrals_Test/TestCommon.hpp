#pragma once
#include <SDE/NWXDefaults.hpp>
#include <catch/catch.hpp>

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

inline Molecule make_molecule()
{
    auto crt = default_runtime();
    auto water = crt.pubchem.at("water");
    auto water_w_basis = crt.apply_basis("sto-3g", water);

    return water_w_basis;
}
