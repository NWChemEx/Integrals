#pragma once
#include <LibChemist/LibChemist.hpp>
#include <catch/catch.hpp>

struct ptr_wrapper{
    const double* ptr_;
    size_t n_;
    const double* begin()const{return ptr_;}
    const double* end()const{return ptr_+n_;}
};

inline void compare_integrals(const ptr_wrapper& calc,
                              const std::vector<double>& corr)
{
    const double eps = 1000*std::numeric_limits<double>::epsilon();
    const double marg = 10*std::numeric_limits<double>::epsilon();
    size_t i = 0;
    for (const double& x : calc) {
        INFO("epsilon = " << eps); 
        REQUIRE(x == Approx(corr[i++]).epsilon(eps).margin(marg));
    }
}

inline auto make_molecule() {
   using LibChemist::Atom;
   using c_t = typename Atom::coord_type;
   Atom H1{1ul, c_t{1.638033502034240, 1.136556880358410, 0.000000000000000}};
   Atom O{8ul, c_t{0.000000000000000, -0.143222342980786, 0.000000000000000}};
   Atom H2{1ul, c_t{-1.638033502034240, 1.136556880358410, 0.000000000000000}};
   LibChemist::Molecule water(O, H1, H2);
   return std::make_tuple(water, LibChemist::apply_basis("sto-3g", water));
}

