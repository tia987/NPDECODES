/**
 * @file radauthreetimesteppingode.cc
 * @brief NPDE homework RadauThreeTimestepping
 * @author Erick Schulz
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "radauthreetimesteppingode.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace RadauThreeTimestepping {

/* SAM_LISTING_BEGIN_1 */
std::vector<double> twoStageRadauTimesteppingLinScalODE(unsigned int m) {
    std::vector<double> sol_vec;
    //====================
    double tau = 5./m;

    auto evo = [](double step){
        double k = (step*(1.+step/6.))/((1.+step*5./12.)*(1.+step/4.)+(step*step/16.));
        return 1-k;
    };
    sol_vec.push_back(1.);
    for(unsigned i = 0; i < m; i++){
        sol_vec.push_back(evo(tau)*sol_vec.back());
    }
    //====================
    return sol_vec;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void testConvergenceTwoStageRadauLinScalODE() {
    constexpr int nIter = 10;       // total number of iterations
    double max_norm_errors[nIter];  // errors vector for all approx. sols
    double rates[nIter - 1];        // The rates of convergence
    double avg_rate = 0.0;  // The average rate of convergence over all iterations

    //====================
    
    std::vector<double> err, conv;
    for(unsigned i = 0; i < nIter; i++){
        unsigned m = 10*std::pow(2,i);
        double step = 5.0/i;
        std::vector<double> exact_sol;
        for(unsigned i = 0; i < m; i++){
            exact_sol.push_back(std::exp(-i*step));
        }
        auto sol = twoStageRadauTimesteppingLinScalODE(m);
        
    }
    //====================
    /* SAM_LISTING_END_2 */

    // Printing results
    std::cout << "\n" << std::endl;
    std::cout << "*********************************************************"
              << std::endl;
    std::cout << "         Convergence of two-stage Radau Method           "
              << std::endl;
    std::cout << "*********************************************************"
              << std::endl;
    std::cout << "--------------------- RESULTS ---------------------------"
              << std::endl;
    std::cout << "Iteration"
              << "\t| Nsteps"
              << "\t| error"
              << "\t\t| rates" << std::endl;
    std::cout << "---------------------------------------------------------"
              << std::endl;
    for (int k = 0; k < nIter; k++) {
      std::cout << k << "\t"
                << "\t|" << 10 * std::pow(2, k) << "\t\t|" << max_norm_errors[k];
      if (k > 0) {
        std::cout << "\t|" << rates[k - 1];
      }
      std::cout << "\n";
    }
    std::cout << "---------------------------------------------------------"
              << std::endl;
    std::cout << "Average rate of convergence: " << avg_rate << "\n" << std::endl;
}

}  // namespace RadauThreeTimestepping
