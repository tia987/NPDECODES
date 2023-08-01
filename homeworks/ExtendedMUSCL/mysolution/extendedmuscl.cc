/**
 * @file extendedmuscl.cc
 * @brief NPDE homework ExtendedMUSCL code
 * @author Oliver Rietmann
 * @date 04.08.2019
 * @copyright Developed at ETH Zurich
 */

#include "extendedmuscl.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>

namespace ExtendedMUSCL {

/* SAM_LISTING_BEGIN_1 */
double logGodunovFlux(double v, double w) {
    double godunov_numerical_flux;
    auto f = [](double u) { return u * (std::log(u) - 1.0); };
    auto df = [](double u) { return std::log(u); };
    //====================
    // calculate sdot
    // auto sdot = [&f](double v, double w){
    //     return (f(v)-f(w))/(v-w);
    // };
    // if(v > w){
    //     if(v/w < sdot(v,w)){
    //         return v;
    //     } else {
    //         return w;
    //     }
    // } else {
    //     if(v/w < df(v)){
    //         return v;
    //     } else if (v/w > df(w)) {
    //         return w;
    //     } else {
    //         return std::exp(v/w);
    //     }
    // }
    if(v > w) return std::max(f(v),f(w));
    else {
        if(df(w) < 0){
            return f(w);
        } else if(df(v) > 0){
            return f(v);
        } else {
            return f(1.);
        }
    }
    //====================
    return godunov_numerical_flux;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_4 */
double limiterMC(double mu_left, double mu_center, double mu_right) {
    double scaled_slope = 0;

    //====================
    double rho_left, rho_center, rho_right;
    rho_left = (mu_right-mu_left)/2.;
    rho_center = (mu_center-mu_left)*2.;
    rho_right = (mu_right-mu_center)*2.;
    if(rho_left > 0 && rho_center > 0 && rho_right > 0) return std::min(rho_left,std::min(rho_center,rho_right));
    else if(rho_left < 0 && rho_center < 0 && rho_right < 0) return std::max(rho_left,std::max(rho_center,rho_right));
    //====================

    return scaled_slope;
}
/* SAM_LISTING_END_4 */

}  // namespace ExtendedMUSCL
