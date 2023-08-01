/**
 * @file solvecauchyproblem.cc
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 19.07.2019
 * @copyright Developed at ETH Zurich
 */

#include "solvecauchyproblem.h"

#include <Eigen/Core>
#include <cmath>

#include "uniformcubicspline.h"

namespace CLEmpiricFlux {

/* SAM_LISTING_BEGIN_1 */
Eigen::Vector2d findSupport(const UniformCubicSpline &f,
                            Eigen::Vector2d initsupp,
                            double t) {
  Eigen::Vector2d result;
  //====================
  Eigen::Vector2d speed = {f.derivative(-1.),f.derivative(1.)};
  result = initsupp+t*speed;  
  //====================
  return result;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
Eigen::VectorXd semiDiscreteRhs(const Eigen::VectorXd &mu0, double h,
                                FUNCTOR &&numFlux) {
  int m = mu0.size();
  Eigen::VectorXd mu1(m);
  //====================
  // Your code goes here
  //====================
  return mu1;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR>
Eigen::VectorXd RalstonODESolver(FUNCTOR &&rhs, Eigen::VectorXd mu0, double tau,
                                 int n) {
  
  //====================
  for(unsigned i = 0; i < n; i++){
    auto k1 = rhs(mu0);
    auto k2 = rhs(mu0+tau*k1*(2/3.));
    mu0 = mu0+tau*(k1*1/4.+k2*3/4.);
  }
  //====================
  
  return mu0;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd solveCauchyProblem(const UniformCubicSpline &f,
                                   const Eigen::VectorXd &mu0, double h,
                                   double T) {
    Eigen::VectorXd muT(mu0.size());

    //====================
    double tau = std::min(h / std::abs(f.derivative(-1.0)), h / std::abs(f.derivative(1.0)));
    GodunovFlux GF(f);
    double n = (int)std::floor(T/tau);
    auto rhs = [h, &GF](Eigen::VectorXd mu){
        return semiDiscreteRhs(mu, h, GF);
    };
    muT = RalstonODESolver(rhs, mu0, tau, n);
    //====================

    return muT;
}
/* SAM_LISTING_END_4 */

}  // namespace CLEmpiricFlux
