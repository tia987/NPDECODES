/**
 * @file conslawwithsource.cc
 * @brief NPDE exam TEMPLATE CODE FILE
 * @author Oliver Rietmann
 * @date 20.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include "conslawwithsource.h"

#include <cmath>

namespace ConsLawWithSource {

/* SAM_LISTING_BEGIN_1 */
double godnfn(double v, double w) {
  auto f = [](double u) { return std::exp(u) - u; };

  //====================
  // Your code goes here
  // Replace the dummy return value below
  // Use 11.3.4.31
  auto df = [](double u){return std::log(std::exp(u)-1);};
  auto s = [f](double v, double w){return (f(v)-f(w))/(v-w);};
  if(v == w || v > w && s(v,w) < 0 || v < w && df(w) < 0){
      return w;
  } else if(v > w && s(v,w) > 0 || v < w && df(v) > 0){
      return v;
  } else if(v < w && df(v) <= 0 && 0 <= df(w)){
      return 1/df(0.);
  }
  return 0.;
  //====================
}
/* SAM_LISTING_END_1 */

}  // namespace ConsLawWithSource
