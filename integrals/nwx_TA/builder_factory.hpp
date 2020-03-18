#include <tiledarray.h>
#include "integrals/nwx_TA/fill_ND_functor.hpp"

namespace nwx_TA {

    // Factory that produces tile-specific builders based on a full array master functor
    template<typename val_type, libint2::Operator op, std::size_t NBases>
    struct BuilderFactory {

        using Builder = FillNDFunctor<val_type, op, NBases>;

        Builder master;

        BuilderFactory(Builder master) : master(master) {};
        ~BuilderFactory() = default;

        /** @brief Produce a builder for the provided tile range.
         *
         *  @param range The TiledArray Range
         *  @returns The builder functor for the integral tile
         */
        Builder operator()(const TA::Range& range) const;

    }; // struct BuilderFactory

    extern template class BuilderFactory<TA::TensorD, libint2::Operator::overlap, 2>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::kinetic, 2>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::nuclear, 2>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::coulomb, 2>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::coulomb, 3>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::coulomb, 4>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::stg, 2>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::stg, 3>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::stg, 4>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::yukawa, 2>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::yukawa, 3>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::yukawa, 4>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::emultipole1, 2>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::emultipole2, 2>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::emultipole3, 2>;
    extern template class BuilderFactory<TA::TensorD, libint2::Operator::delta, 4>;
    
} // namespace nwx_TA