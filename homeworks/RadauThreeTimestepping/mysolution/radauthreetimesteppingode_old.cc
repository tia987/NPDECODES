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
  sol_vec.push_back(1.);
  for(int i = 0; i < m; i++){
    sol_vec.push_back(1-tau*(1+tau/6.)/((1+tau*5/12.)*(1+tau/4.)+tau*tau/16.)*sol_vec.at(i));
  }
  //====================
  return sol_vec;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void testConvergenceTwoStageRadauLinScalODE() {
  constexpr int nIter = 10;       // total number of iterations
  std::vector<double>  max_norm_errors(nIter);  // errors vector for all approx. sols
  std::vector<double>  rates(nIter-1);        // The rates of convergence
  double avg_rate = 0.0;  // The average rate of convergence over all iterations

  //====================
  for(unsigned int i = 0; i < nIter; i++){
    // We use the exact solution of the ODE
    unsigned m = 10*std::pow(2,i);
    std::vector<double> real_vec;
    std::vector<double> sol_vec;
    for(unsigned int k = 0; k < m+1; k++){
      real_vec.push_back(std::exp(-k*5./m));
    }
    // To extract the approximation from the sol_vec with the 2 stage Radau
    sol_vec = twoStageRadauTimesteppingLinScalODE(m);
    for(unsigned j = 0; j < m; j++){
      double diff = std::abs(real_vec.at(j)-sol_vec.at(j));
      if(max_norm_errors.at(i) < diff) max_norm_errors.at(i) = diff;
    }    
  }
  // And calculate the error
  for(unsigned i = 1; i < nIter; i++){
    rates.at(i-1) = std::abs(std::abs(max_norm_errors.at(i))-std::abs(max_norm_errors.at(i-1)));
    avg_rate += rates.at(i-1);
  }
  avg_rate /= nIter-1;
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
