/**
 * @file laxwendroffscheme.cc
 * @brief NPDE homework "LaxWendroffScheme" code
 * @author Oliver Rietmann
 * @date 29.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "laxwendroffscheme.h"

#include <Eigen/Core>
#include <cmath>

namespace LaxWendroffScheme {

namespace Constant {
constexpr double e = 2.71828182845904523536;
constexpr double pi = 3.14159265358979323846;
}  // namespace Constant

constexpr double Square(double x) { return x * x; }

/**
 * @brief Computes the right-hand side according to the Lax-Wendroff scheme.
 * @param mu mu^(k-1) (i.e. mu at timestep k-1)
 * @param gamma tau / h, where tau = timestep size and h = spatial meshwidth
 * @return mu^(k) (i.e. mu at timestep k)
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::VectorXd LaxWendroffRhs(const Eigen::VectorXd &mu, double gamma) {
    int N = mu.size();
    Eigen::VectorXd result(N);

    //====================
    auto f = [](double x){return std::exp(x);};
    double gamma2 = gamma*gamma;
    result(0) = mu(0)-gamma*(f(mu(1))-f(mu(0)))/2.+
                    gamma2*(f((mu(1)+mu(0))/2.)*f((mu(1)+mu(0))/2.)*(mu(1)-mu(0))-
                            f((mu(1)+mu(0))/2.)*f((mu(1)+mu(0))/2.)*(mu(1)-mu(0)))/2.;
    for(unsigned i = 1; i < N-1; i++){
        result(i) = mu(i)-gamma*(f(mu(i+1))-f(mu(i-1)))/2.+
                    gamma2*(f((mu(i+1)+mu(i))/2.)*f((mu(i+1)+mu(i))/2.)*(mu(i+1)-mu(i))-
                            f((mu(i)+mu(i-1))/2.)*f((mu(i)+mu(i-1))/2.)*(mu(i)-mu(i-1)))/2.;
    }
    result(N-1) = mu(N-1)-gamma*(f(mu(N-1))-f(mu(N-2)))/2.+
                gamma2*(f((mu(N-1)+mu(N-2))/2.)*f((mu(N-1)+mu(N-2))/2.)*(mu(N-1)-mu(N-2))-
                        f((mu(N-1)+mu(N-2))/2.)*f((mu(N-1)+mu(N-2))/2.)*(mu(N-1)-mu(N-2)))/2.;
    //====================

    return result;
}

Eigen::VectorXd solveLaxWendroff(const Eigen::VectorXd &u0, 
                                 double T,
                                 unsigned int M) {
    double gamma = 1.0 / Constant::e;
    Eigen::VectorXd mu = u0;
    // Main timestepping loop
    for (int j = 0; j < M; ++j) mu = LaxWendroffRhs(mu, gamma);
    return mu;
}

/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_2 */
// Build spatial grid
Eigen::VectorXd getXValues(double T, unsigned int M) {
  double tau = T / M;
  double h = Constant::e * tau;
  int j_max = (int)(std::ceil((3.0 * T + 1.0) / h) + 0.5);
  int j_min = (int)(std::floor(-3.0 * T / h) - 0.5);
  unsigned int N = j_max - j_min + 1;
  return Eigen::VectorXd::LinSpaced(N, j_min * h, j_max * h);
}
Eigen::VectorXd numexpLaxWendroffRP(const Eigen::VectorXi &M) {
    const double T = 1.0;
    const int M_size = M.size();
    Eigen::VectorXd error(M_size);
    // Initial values for the Riemann problem
    auto u_initial = [](double x) { return 0.0 <= x ? 1.0 : 0.0; };
    // Exact solution \prbeqref{eq:solrp} at time $T = 1.0$
    auto u_exact = [](double x) {
      return (x <= 1.0) ? 0.0 : ((Constant::e <= x) ? 1.0 : std::log(x));
    };
  
    //====================    
    for(unsigned i = 0; i < M_size; i++){
        Eigen::VectorXd x = getXValues(T, M(i));
        Eigen::VectorXd u0 = x.unaryExpr(u_initial);
        auto u_approx = LaxWendroffScheme::solveLaxWendroff(u0, T, M(i));
        double tau = T/M(i);
        double h = Constant::e*tau;
        error(i) = h*(x.unaryExpr(u_exact)-u_approx).lpNorm<1>();
    }
    //====================
    return error;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Evaluates the discrete function u at position x by linear
 * interpolation
 * @param u descrete function values at spatial positions y
 * @param y vector of same length as u, representing the nodes of u
 * @return best linear interpolation of u at spacial position x
 */
double eval(const Eigen::VectorXd &u, const Eigen::VectorXd &y, double x) {
  int N = y.size();
  double a = y(0);
  double b = y(N - 1);

  if (x <= a) return u(0);
  if (b <= x) return u(N - 1);

  double lambda = (x - a) / (b - a);
  int k0 = (int)(lambda * (N - 1));
  int k1 = k0 + 1;

  lambda = (x - y(k0)) / (y(k1) - y(k0));
  return lambda * u(k1) + (1.0 - lambda) * u(k0);
}

/* SAM_LISTING_BEGIN_9 */
double smoothU0(double x) {
  return (x < 0.0)
             ? 0.0
             : ((1.0 < x) ? 1.0 : Square(std::sin(0.5 * Constant::pi * x)));
}
Eigen::VectorXd referenceSolution(const Eigen::VectorXd &x) {
    double T = 1.0;
    // Reference solution on a very fine mesh
    unsigned int M = 3200;

    Eigen::VectorXd y = getXValues(T, M);
    Eigen::VectorXd u0 = y.unaryExpr(&smoothU0);
    Eigen::VectorXd u = solveLaxWendroff(u0, T, M);
    int N = x.size();
    Eigen::VectorXd u_ref(N);
    // The vector u is larger than u_ref. Use eval() from above the "evaluate" u
    // at the positions x(i) and thus obtain the reference solution u_ref.
    //====================
    for(unsigned i = 0; i < N; i++){
        u_ref(i) = eval(u,y,x(i));
    }      
    //===================
    return u_ref;
}
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd numexpLaxWendroffSmoothU0(const Eigen::VectorXi &M) {
    const double T = 1.0;
    const int M_size = M.size();
    Eigen::VectorXd error(M_size);

     // Initial values for the Riemann problem
    auto u_initial = [](double x) { 
        if(x < 0) return 0.;
        else if(x > 1) return 1.;
        else return std::sin(M_PI/2.*x)*std::sin(M_PI/2.*x);
    };
    // Exact solution \prbeqref{eq:solrp} at time $T = 1.0$
    
  
    //====================
    //auto u_ref = referenceSolution(u_approx);
    for(unsigned i = 0; i < M_size; i++){
        Eigen::VectorXd x = getXValues(T, M(i));
        Eigen::VectorXd u0 = x.unaryExpr(&smoothU0);
        auto u_approx = solveLaxWendroff(u0, T, M(i));
        double tau = T/M(i);
        double h = Constant::e*tau;
        error(i) = h*(referenceSolution(x)-u_approx).lpNorm<1>();
    }
    //====================
    
    return error;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_7 */
Eigen::VectorXd solveGodunov(const Eigen::VectorXd &u0, 
                             double T,
                             unsigned int M) {
    double tau = T / M;
    double h = Constant::e * tau;
    unsigned int N = u0.size();
    Eigen::VectorXd mu = u0;

    //====================
    for(unsigned i = 0; i < M; i++){
        for(unsigned j = N-1; j > 0; j--){
            mu(j) -= tau/h*(std::exp(mu(j))-std::exp(mu(j-1)));;
        }
    }
    //====================
    return mu;
}

/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_8 */
Eigen::VectorXd numexpGodunovSmoothU0(const Eigen::VectorXi &M) {
    const double T = 1.0;
    const int M_size = M.size();
    Eigen::VectorXd error(M_size);

    //====================
    for(unsigned i = 0; i < M_size; i++){
        Eigen::VectorXd x = getXValues(T, M(i));
        Eigen::VectorXd u0 = x.unaryExpr(&smoothU0);
        auto u_approx = solveGodunov(u0, T, M(i));
        double tau = T/M(i);
        double h = Constant::e*tau;
        error(i) = h*(referenceSolution(x)-u_approx).lpNorm<1>();
    }
    //====================

    return error;
}
/* SAM_LISTING_END_8 */

}  // namespace LaxWendroffScheme
