/**
 * @file engquistoshernumericalflux.cc
 * @brief NPDE homework "EngquistOsherNumericalFlux" code
 * @author Oliver Rietmann
 * @date 23.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "engquistoshernumericalflux.h"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>

namespace EngquistOsherNumericalFlux {

/* SAM_LISTING_BEGIN_1 */
double EngquistOsherNumFlux(double v, double w) {
    double result;

    //====================
    auto f = [](double x){return std::cosh(x);};
    if(v >= 0 && w >= 0) result = f(v);
    else if(v < 0 && w >= 0) result = f(0);
    else if(v >= 0 && w < 0) result = f(v)+f(w)-f(0);
    else result = f(w);
    //====================

    return result;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solveCP(double a, double b, Eigen::VectorXd u0, double T) {
    // Find the maximal speed of propagation
    double A = u0.minCoeff();
    double B = u0.maxCoeff();
    double K = std::max(std::abs(std::sinh(A)), std::abs(std::sinh(B)));
    // Set uniform timestep according to CFL condition
    int N = u0.size();
    double h = (b - a) / N;
    double tau_max = h / K;
    double timesteps = std::ceil(T / tau_max);
    double tau = T / timesteps;

    // Main timestepping loop
    //====================
    Eigen::VectorXd mu(N);
    for(unsigned j = 0; j < timesteps; j++){
        mu.swap(u0);
        double Fm;
        double Fp = EngquistOsherNumFlux(mu(0),mu(0));
        for(unsigned i = 0; i < N-1; i++){
            Fm = Fp;
            Fp = EngquistOsherNumFlux(mu(i),mu(i+1));
            u0(i) = mu(i)-tau/h*(Fp-Fm);
        }
        Fm = Fp;
        Fp = EngquistOsherNumFlux(mu(N-1),mu(N-1));
        u0(N-1) = mu(N-1)-tau/h*(Fp-Fm);
    }
    //====================

    return u0;
}
/* SAM_LISTING_END_2 */

}  // namespace EngquistOsherNumericalFlux
