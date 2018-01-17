#pragma once
#include <LibChemist/SetOfAtoms.hpp>

namespace IntegralsEx{

/*! @brief One-electron integral implementation
 */
class TwoCenterIntegral
{
public:
        typedef uint64_t ShellIndex;

        /*! @brief initialize the integral computation
         *
         * @param [in] deriv Derivative to calculate
         * @param [in] atoms SetofAtoms to use to calculate integrals
         * @param [in] bs1 BasisSet to use on the first center
         * @param [in] bs2 BasisSet to use on the second center
         */
        TwoCenterIntegral(unsigned int deriv,
                          const LibChemist::SetOfAtoms &atoms,
                          const LibChemist::BasisSet &bs1,
                          const LibChemist::BasisSet &bs2){}

        virtual ~TwoCenterIntegral(){}

        /*! Return the number of components calculated by this module
         *
         * For example, something that calculates x,y,z component would return 3
         */
        virtual unsigned int n_components(void) const {return 1;}

        /*! @brief calculate an integral
         *
         * @param [in] shell1 Shell index on the first center
         * @param [in] shell2 Shell index on the second center
         * @returns A pointer to the beginning of the integral buffer
         */
        virtual const double* calculate(ShellIndex shell1, ShellIndex shell2) = 0;

        /*! @brief calculate multiple integrals
         *
         * @param [in] shells1 Shell indices on the first center
         * @param [in] shells2 Shell indices on the second center
         * @returns A pointer to the beginning of the integral buffer
         */
        virtual const double* calculate_multi(const std::vector<ShellIndex> & shells1,
                                 const std::vector<ShellIndex> & shells2)
        {
            return nullptr;
        }
};

} // close namespace IntegralsEx

