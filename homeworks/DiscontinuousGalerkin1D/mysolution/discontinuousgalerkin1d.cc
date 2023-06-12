/**
 * @file discontinuousgalerkin1d.cc
 * @brief NPDE homework "DiscontinuousGalerkin1D" code
 * @author Oliver Rietmann
 * @date 22.05.2019
 * @copyright Developed at ETH Zurich
 */

#include "discontinuousgalerkin1d.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace DiscontinuousGalerkin1D {

/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<double> compBmat(int Ml, int Mr, double h) {
    const int N = 2 * (Ml + Mr + 1);
    Eigen::SparseMatrix<double> A(N, N);
    //====================
    A.reserve(Eigen::VectorXi::Ones(N));
    for (int i = 0; i < N; i += 2) {
        A.insert(i, i) = h;
        A.insert(i + 1, i + 1) = h * h * h / 12.0;
    }
    //====================
    return A;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double Feo(double v, double w) {
    double result;
    //====================
    
    //====================
    return result;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
Solution solveTrafficFlow() {
  int Ml = 40;
  int Mr = 40;
  int N_half = Mr + Ml + 1;
  int N = 2 * N_half;

  double h = 0.05;
  double tau = h / 3;
  double T = 1.0;
  unsigned int m = (unsigned int)(T / tau);


  //====================
  // Your code goes here
  // Fill the following vectors
  Eigen::VectorXd x;
  Eigen::VectorXd u;
  //====================

  return Solution(std::move(x), std::move(u));
}
/* SAM_LISTING_END_3 */

}  // namespace DiscontinuousGalerkin1D
