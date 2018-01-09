#include "TestCommon.hpp"
#include "nwx_simint.hpp"
#include <iostream>

//Computes the ERI integrals for water in STO-3G
int main(){
    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);

    simint_init();
   
    std::vector<simint_shell> shells = nwx_simint::make_basis(bs);

    for (int x=0; x<5; x++) {
      std::cout << "memsize = " << shells[x].memsize << std::endl;
      //std::cout << "X = " << shells[x].x << std::endl;
      for (int y=0; y<3; y++) {
        std::cout << "alpha = " << shells[x].alpha[y];
        std::cout << "   coef = " << shells[x].coef[y] << std::endl;
      }
    }

    // we need the starting indices of the individual basis functions within the
    // shells, as well as how many basis functions are in each shell.

    //                     s  s  p  s  s
    int shell_start[5] = { 0, 1, 2, 5, 6 };
    int shell_nbf[5]   = { 1, 1, 3, 1, 1 };


       ////////////////////////////////////////
    // 4. Allocate some memory
    // Final (packed) storage
    // This example has 7 basis functions. That would
    // be 7^4 = 2401 unpacked integrals, or 406 packed integrals
    double * integrals = (double *)SIMINT_ALLOC(406 * sizeof(double));

    // This is where individual calls to simint_compute_eri will place
    // its integrals.
    // This must be big enough to hold an integral computation for the highest
    // angular momentum in the basis set (multiplied by number of batches, however
    // we aren't doing batches here).
    // For water/sto-3g, highest am = p, so we have to hold
    // a (p p | p p) cartesian integral, which is 3*3*3*3 = 81 elements
    double * target = (double *)SIMINT_ALLOC(81 * sizeof(double));


    // Shared workspace. This prevents needless allocation/deallocation
    // for each call
    //
    // No need to multiply by sizeof(double) here. It is included
    // in the return value from simint_ostei_workmem
    //
    // Arguments to simint_ostei_workmem: Derivative order, max am
    //
    // (for your information, OSTEI = Obara-Saika Two-Electron Integral)
    const size_t simint_workmem = simint_ostei_workmem(0,1);
    double * work = (double *)SIMINT_ALLOC(simint_workmem);


    ////////////////////////////////////////
    // 5. Loop over the shells
    // The following is somewhat naive due
    // to only calculating the unique shell
    // quartets, but is just a demonstration
    ////////////////////////////////////////

    // number of shell quartets we calculated
    int ncomputed_shell = 0;

    // number of integrals calcualated
    int ncomputed_integrals = 0;

    // we only have to initialize these once
    simint_multi_shellpair left_pair;
    simint_multi_shellpair right_pair;
    simint_initialize_multi_shellpair(&left_pair);
    simint_initialize_multi_shellpair(&right_pair);

    for(int i = 0; i < 5; i++)
    for(int j = 0; j <= i; j++)
    {
        int ij = (i*(i+1))/2 + j;

        // form the left shell pair
        simint_create_multi_shellpair(1, &shells[i],
                                      1, &shells[j], &left_pair, 0);
        
        std::cout << "am1 = " << left_pair.am1 << "   am2 = " << left_pair.am2
                  << "   nshell = " << left_pair.nshell12 << std::endl;

        for(int k = 0; k < 5; k++)
        for(int l = 0; l <= k; l++)
        {
            int kl = (k*(k+1))/2 + l;
            if(ij < kl) // skip due to permutational symmetry?
                continue;

            simint_create_multi_shellpair(1, &shells[k],
                                          1, &shells[l], &right_pair, 0);


            /////////////////////////////////////
            // Integrals actually computed here
            /////////////////////////////////////
            
            // Should always return one in this case (since we aren't batching)
            simint_compute_eri(&left_pair, &right_pair, 0.0, work, target);

            // keep track of how many shells and integrals we calculated
            ncomputed_shell++;
            ncomputed_integrals += shell_nbf[i] * shell_nbf[j] * shell_nbf[k] * shell_nbf[l];


            // index of where we are in the intermediate "target" buffer
            int target_idx = 0;

            // Now, determine where the integrals go in the final integral storage and put them there
            for(int m = 0; m < shell_nbf[i]; m++)
            for(int n = 0; n < shell_nbf[j]; n++)
            for(int o = 0; o < shell_nbf[k]; o++)
            for(int p = 0; p < shell_nbf[l]; p++)
            {
                int m_idx = shell_start[i] + m;
                int n_idx = shell_start[j] + n;
                int o_idx = shell_start[k] + o;
                int p_idx = shell_start[l] + p;

                int mn_idx = m_idx < n_idx ? (n_idx*(n_idx+1))/2 + m_idx : (m_idx*(m_idx+1))/2 + n_idx;
                int op_idx = o_idx < p_idx ? (p_idx*(p_idx+1))/2 + o_idx : (o_idx*(o_idx+1))/2 + p_idx;
                int mnop_idx = mn_idx < op_idx ? (op_idx*(op_idx+1))/2 + mn_idx : (mn_idx*(mn_idx+1))/2 + op_idx;

                integrals[mnop_idx] = target[target_idx];
                std::cout << target[target_idx] << " ";
                target_idx++;

            }
        std::cout << std::endl;
        }
    }

    //////////////////////////
    // Print out some info
    //////////////////////////
    printf("\n");
    printf("    Number of shell quartets (unpacked): %d\n", 5*5*5*5);
    printf("  Number of shell quartets (calculated): %d\n", ncomputed_shell);
    printf("         Number of integrals (unpacked): %d\n", 7*7*7*7);
    printf("       Number of integrals (calculated): %d\n", ncomputed_integrals);
    printf("\n");


    // cleanup
    simint_free_multi_shellpair(&left_pair);
    simint_free_multi_shellpair(&right_pair);

    for(int i = 0; i < 5; i++)
        simint_free_shell(&shells[i]);

    SIMINT_FREE(target);
    SIMINT_FREE(work);
    SIMINT_FREE(integrals);

    simint_finalize(); 
}
