#include "builder_factory.hpp"
#include <integrals/nwx_TA/nwx_TA_utils.hpp>

namespace nwx_TA {

    template<typename val_type, libint2::Operator op, std::size_t NBases>
    using Builder = FillNDFunctor<val_type, op, NBases>;

    template<typename val_type, libint2::Operator op, std::size_t NBases>
    Builder<val_type, op, NBases> BuilderFactory<val_type, op, NBases>::operator()(const TA::Range& range) const {
        auto builder = Builder(); // The builder to be returned

        std::vector<libint2::BasisSet> basis_subsets = {}; // Subsets of full basis sets that are needed for tile

        for (int depth = 0; depth < NBases; depth++) {
            int tile_depth = (master.nopers == 1) ? depth : depth + 1; // Deals with multipole cases
            auto depth_shells = nwx_TA::aos2shells(master.LIBasis_sets[depth],
                    range.lobound()[tile_depth], range.upbound()[tile_depth]); // Get indices for shells at this depth

            libint2::BasisSet current_subset; // Shells subset for this depth

            for (auto shell : depth_shells) {
                current_subset.emplace_back(master.LIBasis_sets[depth][shell]);
            }

            basis_subsets.emplace_back(current_subset);
        }

        // Initialize the new builder with the sub-basis sets and master's parameters.
        builder.initialize(basis_subsets, master.factory.deriv, master.factory.thresh, 0.0);

        // New builder factory special parameters to master (STG exponent, origin, charges...)
        builder.factory.stg_exponent = master.factory.stg_exponent;
        builder.factory.origin = master.factory.origin;
        builder.factory.qs = master.factory.qs;

        // Get Cauchy-Schwarz information from master, it's only calculated once.
        master.screen.set_sub_screen(master.LIBasis_sets, range, builder.screen.cs_mat1, builder.screen.cs_mat2);
        builder.cs_thresh = master.cs_thresh;

        return builder;
    }

    template class BuilderFactory<TA::TensorD, libint2::Operator::overlap, 2>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::kinetic, 2>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::nuclear, 2>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::coulomb, 2>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::coulomb, 3>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::coulomb, 4>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::stg, 2>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::stg, 3>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::stg, 4>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::yukawa, 2>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::yukawa, 3>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::yukawa, 4>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::emultipole1, 2>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::emultipole2, 2>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::emultipole3, 2>;
    template class BuilderFactory<TA::TensorD, libint2::Operator::delta, 4>;

}