#include "libint/detail_/nwx_libint.hpp"
#include <integrals/factory.hpp>

size_t integrals::Factory::instances;

integrals::Factory::Factory(integrals::property prop, integrals::NWX_basis bs1,
                            integrals::NWX_basis bs2, integrals::NWX_basis bs3,
                            integrals::NWX_basis bs4) {
    if(instances == 0 && not libint2::initialized()) { libint2::initialize(); };
    instances++;
    libint2::Operator op;
    switch(prop) {
        case integrals::property::overlap:
            op = libint2::Operator::overlap;
            break;
        case integrals::property::kinetic:
            op = libint2::Operator::kinetic;
            break;
        case integrals::property::nuclear:
            op = libint2::Operator::nuclear;
            break;
        case integrals::property::eri: op = libint2::Operator::coulomb; break;
        default: assert(false && "invalid integral case in switch");
    };
    lbs1            = nwx_libint::make_basis(bs1);
    lbs2            = nwx_libint::make_basis(bs2);
    lbs3            = nwx_libint::make_basis(bs3);
    lbs4            = nwx_libint::make_basis(bs4);
    auto max_nprims = nwx_libint::sets_max_nprims({lbs1, lbs2, lbs3, lbs4});
    auto max_l      = nwx_libint::sets_max_l({lbs1, lbs2, lbs3, lbs4});
    double thresh   = 1.0e-16;
    size_t deriv    = 1;
    engine.reset(new libint2::Engine(op, max_nprims, max_l, thresh, deriv));
};

integrals::Factory::Factory(integrals::property prop, integrals::NWX_basis bs1,
                            integrals::NWX_basis bs2,
                            integrals::NWX_basis bs3) :
  Factory(prop, bs1, bs2, bs3, bs3){};

integrals::Factory::Factory(integrals::property prop, integrals::NWX_basis bs1,
                            integrals::NWX_basis bs2) :
  Factory(prop, bs1, bs2, bs2){};

integrals::Factory::~Factory() {
    // delete engine; Apparently std::unique_ptr self destructs.
    instances--;
    if(instances == 0 && libint2::initialized()) { libint2::finalize(); };
};

const double* const integrals::Factory::compute(size_t s1, size_t s2) {
    engine->compute(lbs1[s1], lbs2[s2]);
    return engine->results()[0];
};

const double* const integrals::Factory::compute(size_t s1, size_t s2, size_t s3,
                                                size_t s4) {
    engine->compute(lbs1[s1], lbs2[s2], lbs3[s3], lbs4[s4]);
    return engine->results()[0];
};
