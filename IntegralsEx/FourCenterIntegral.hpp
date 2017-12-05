#pragma once
#include <LibChemist/SetOfAtoms.hpp>

namespace IntegralsEx{

/*! @brief Two-electron integral implementation
 *
 */
class FourCenterIntegral
{

public:

        /*! @brief initialize the integral computation
         *
         * @param [in] deriv Derivative to calculate
         * @param [in] atoms SetOfAtoms to use to calculate integrals
         * @param [in] bs1 BasisSet to use for the first center
         * @param [in] bs2 BasisSet to use for the second center
         * @param [in] bs3 BasisSet to use for the third center
         * @param [in] bs4 BasisSet to use for the fourth center
         */
        FourCenterIntegral(unsigned int deriv,
                           const LibChemist::SetOfAtoms &atoms,
                           const LibChemist::BasisSet & bs1,
                           const LibChemist::BasisSet & bs2,
                           const LibChemist::BasisSet & bs3,
                           const LibChemist::BasisSet & bs4,
                           double thresh){}

        virtual ~FourCenterIntegral(){}

        /*! Return the number of components calculated by this module
         *
         * For example, something that calculates x,y,z component would return 3
         */
        unsigned int n_components(void) const
        {
            return n_components_();
        }

        /*! @brief calculate an integral
         *
         * @param [in] shell1 Shell index on the first center
         * @param [in] shell2 Shell index on the second center
         * @param [in] shell3 Shell index on the third center
         * @param [in] shell4 Shell index on the fourth center
         * @param [in] outbuffer Where to place the completed integrals
         * @param [in] bufsize Size of \p outbuffer (as the number of doubles)
         * @returns A pointer to the beginning of the integral buffer
         */
        const double* calculate(uint64_t shell1, uint64_t shell2,
                           uint64_t shell3, uint64_t shell4)
        {
            return calculate_(shell1,shell2,shell3,shell4);
        }

        /*! @brief calculate multiple integrals
         *
         * @param [in] shells1 Shell indicies on the first center
         * @param [in] shells2 Shell indicies on the second center
         * @param [in] shells3 Shell indicies on the third center
         * @param [in] shells4 Shell indicies on the fourth center
         * @param [in] outbuffer Where to place the completed integrals
         * @param [in] bufsize Size of \p outbuffer (as the number of doubles)
         * @returns Number of integrals calculated
         */
        const double* calculate_multi(const std::vector<uint64_t> & shells1,
                                 const std::vector<uint64_t> & shells2,
                                 const std::vector<uint64_t> & shells3,
                                 const std::vector<uint64_t> & shells4)
        {
            return calculate_multi_(shells1, shells2, shells3, shells4);
        }

        /////////////////////////////////////////
        // To be implemented by derived classes
        /////////////////////////////////////////

        //! @copydoc n_components
        virtual unsigned int n_components_(void) const
        {
            return 1;
        }


        //! @copydoc calculate
        virtual const double* calculate_(uint64_t shell1, uint64_t shell2,
                                    uint64_t shell3, uint64_t shell4) = 0;

        //! \copydoc calculate_multi
        virtual const double* calculate_multi_(const std::vector<uint64_t> & /*shells1*/,
                                          const std::vector<uint64_t> & /*shells2*/,
                                          const std::vector<uint64_t> & /*shells3*/,
                                          const std::vector<uint64_t> & /*shells4*/)
        {
            return nullptr;
        }

};

} // close namespace pulsar

