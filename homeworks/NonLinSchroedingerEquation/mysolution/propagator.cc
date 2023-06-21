/**
 * @file propagator.h
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 04.05.2020
 * @copyright Developed at ETH Zurich
 */

#include "propagator.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <complex>

namespace NonLinSchroedingerEquation {

// KineticPropagator
/* SAM_LISTING_BEGIN_1 */
KineticPropagator::KineticPropagator(const SparseMatrixXd &A,
                                     const SparseMatrixXcd &M, double tau) {
    //====================
    B_ = M+0.5*tau*A.cast<std::complex<double>>();
    SparseMatrixXcd M_ = M-0.5*tau*A.cast<std::complex<double>>();
    solver_.compute(M_); 
    //====================
}

Eigen::VectorXcd KineticPropagator::operator()
  (const Eigen::VectorXcd &mu) const {
    //====================
    return solver_.solve(B_*mu);
    //====================
    
}
/* SAM_LISTING_END_1 */

// InteractionPropagator
/* SAM_LISTING_BEGIN_2 */
InteractionPropagator::InteractionPropagator(double tau) {
    //====================
    phase_ = [&tau](std::complex<double> z){
        const std::complex<double> i(0,1);
        return std::exp(-i*std::norm(z)*tau)*z;
    };
    //====================
}

Eigen::VectorXcd InteractionPropagator::operator()
  (const Eigen::VectorXcd &mu) const {
    //====================
    // Your code goes here
    // Replace mu by its value after a timestep tau    
    return mu.unaryExpr(phase_);
    //====================
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
//====================
// Your code goes here
// Change this dummy implementation of the constructor:
SplitStepPropagator::SplitStepPropagator(const SparseMatrixXd &A,
                                         const SparseMatrixXcd &M, double tau) :
                                         kineticPropagator_ (A, M, 0.5 * tau) , interactionPropagator_ (tau) {
}
//====================

Eigen::VectorXcd SplitStepPropagator::operator()
  (const Eigen::VectorXcd &mu) const {
    Eigen::VectorXcd nu(mu.size());
    //====================
    nu = kineticPropagator_(mu);
    nu = interactionPropagator_(nu);
    nu = kineticPropagator_(nu);
    //====================
    return nu;
}
/* SAM_LISTING_END_3 */

}  // namespace NonLinSchroedingerEquation
