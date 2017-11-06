#pragma once
#include <LibChemist/SetOfAtoms.hpp>


namespace IntegralsEx{

/*! @brief Two-electron, three center integral implementation
 *
 */
class ThreeCenterIntegral
{
    public:

        /*! @brief initialize the integral computation
         *
         * @param [in] deriv Derivative to calculate
         * @param [in] atoms SetOfAtoms to use to calculate integrals
         * @param [in] bs1 BasisSet to use for the first center
         * @param [in] bs2 BasisSet to use for the second center
         * @param [in] bs3 BasisSet to use for the third center
         */
        ThreeCenterIntegral(unsigned int deriv,
                            const LibChemist::SetOfAtoms &atoms,
                            const LibChemist::BasisSet &bs1,
                            const LibChemist::BasisSet &bs2,
                            const LibChemist::BasisSet &bs3){}
        
        virtual ~ThreeCenterIntegral(){}

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
         * @return A pointer to the beginning of the integral buffer
         */
        const double* calculate(uint64_t shell1, uint64_t shell2,
                           uint64_t shell3)
        {
            return calculate_(shell1,shell2,shell3); 
        }

        /*! @brief calculate multiple integrals
         *
         * @param [in] shells1 Shell indicies on the first center
         * @param [in] shells2 Shell indicies on the second center
         * @param [in] shells3 Shell indicies on the third center
         * @return A pointer to the beginning of the integral buffer
         */
        const double* calculate_multi(const std::vector<uint64_t> &shells1,
                                 const std::vector<uint64_t> &shells2,
                                 const std::vector<uint64_t> &shells3)
        {
            return calculate_multi_(shells1,shells2,shells3);
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
                                    uint64_t shell3) = 0;

        //! @copydoc calculate_multi
        virtual const double* calculate_multi_(const std::vector<uint64_t> & /*shells1*/,
                                          const std::vector<uint64_t> & /*shells2*/,
                                          const std::vector<uint64_t> & /*shells3*/)
        {
//            //////////////////////////////////////////////////////////
//            // default implementation - just loop over and do them all
//            //////////////////////////////////////////////////////////

//            uint64_t ntotal = 0;

//            for(uint64_t s1 : shells1)
//            for(uint64_t s2 : shells2)
//            for(uint64_t s3 : shells3)
//            for(uint64_t s4 : shells4)
//            {
//                uint64_t nbatch = calculate_(s1, s2, s3, s4);
//                ntotal += nbatch;
//                outbuffer += nbatch;

//                // be safe with unsigned types
//                bufsize = ( (nbatch >= bufsize) ? 0 : (bufsize - nbatch) );
//            }

            return nullptr;
        }

};

} // close namespace IntegralsEx
