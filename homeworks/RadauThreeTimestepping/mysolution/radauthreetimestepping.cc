/**
 * @file radauthreetimestepping.cc
 * @brief NPDE homework RadauThreeTimestepping
 * @author Erick Schulz, edited by Oliver Rietmann
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "radauthreetimestepping.h"

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <cmath>
#include <iostream>
#include <unsupported/Eigen/KroneckerProduct>

namespace RadauThreeTimestepping {

/**
 * @brief Implementation of the right hand side (time dependent) source vector
 * for the parabolic heat equation
 * @param dofh A reference to the DOFHandler
 * @param time The time at which to evaluate the source vector
 * @returns The source vector at time `time`
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd rhsVectorheatSource(const lf::assemble::DofHandler &dofh,
                                    double time) {
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  // Right-hand side vector has to be set to zero initially
  Eigen::VectorXd phi(N_dofs);
  //====================
  // Build the rhs function
  phi.setZero();
  auto f = [time](Eigen::Vector2d x)->double {
    Eigen::Vector2d v(std::cos(M_PI*time),
                      std::sin(M_PI*time));
    if((x-v*0.5).norm() < 0.5) return 1.;
    else return 0.;
  };
  // Build the Entity vector and assemble it
  TrapRuleLinFEElemVecProvider<decltype(f)> Entity_Vector(f);
  AssembleVectorLocally(0, dofh, Entity_Vector, phi);
  // Now we need to assemble the vector for the boundary
  auto mesh_p = dofh.Mesh();
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  for(auto vertex : mesh_p->Entities(1)){
    if(bd_flags(*vertex)){
      auto index = dofh.GlobalDofIndices(*vertex);
      phi(index[0]) = 0.;
    }
  }
  //====================
  return phi;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Heat evolution solver: the solver obtains the
 * discrete evolution operator from the Radau3MOLTimestepper class and
 * repeatedly iterates its applicaiton starting from the initial condition
 * @param dofh The DOFHandler object
 * @param m is total number of steps until final time final_time (double)
 * @param final_time The duration for which to solve the PDE
 * @returns The solution at the final timestep
 */
/* SAM_LISTING_BEGIN_6 */
Eigen::VectorXd solveHeatEvolution(const lf::assemble::DofHandler &dofh,
                                   unsigned int m, double final_time) {
  Eigen::VectorXd discrete_heat_sol(dofh.NumDofs());
  //====================
  double tau = final_time/m;
  Radau3MOLTimestepper sol(dofh);
  Eigen::VectorXd curr = sol.discreteEvolutionOperator(0.0, tau, Eigen::VectorXd::Zero(dofh.NumDofs()));
  for(unsigned i = 1; i < m; i++){
    curr = sol.discreteEvolutionOperator(i*tau, tau, curr);
  }
  discrete_heat_sol = curr;
  //====================
  return discrete_heat_sol;
}
/* SAM_LISTING_END_6 */

/* Implementing member function Eval of class LinFEMassMatrixProvider*/
Eigen::Matrix<double, 3, 3> LinFEMassMatrixProvider::Eval(
    const lf::mesh::Entity &tria) {
  Eigen::Matrix<double, 3, 3> elMat;
  //====================
  // Throw error in case no triangular cell
  LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
		  "Unsupported cell type " << tria.RefEl());
  // Obtain vertex coordinates of the triangle in a 2x3 matrix
  const double area = lf::geometry::Volume(*(tria.Geometry()));		  
  elMat << 2.0, 1.0, 1.0,
           1.0, 2.0, 1.0,
           1.0, 1.0, 2.0;
  // clang-format on
  elMat *= area / 12.0;;
  //====================
  return elMat;  // return the local mass element matrix
}

/* Implementing constructor of class Radau3MOLTimestepper */
/* SAM_LISTING_BEGIN_4 */
Radau3MOLTimestepper::Radau3MOLTimestepper(const lf::assemble::DofHandler &dofh)
    : dofh_(dofh) {
  //====================
  auto N_dofs = dofh.NumDofs();
  lf::assemble::COOMatrix<double> A_COO(N_dofs,N_dofs);
  lf::assemble::COOMatrix<double> M_COO(N_dofs,N_dofs);
  Eigen::Matrix2d I = Eigen::Matrix2d::Identity();
  Eigen::Matrix2d Coeff;
  Coeff << 5/12., -1/12., 3/4.,  1/4.;      
  lf::uscalfe::LinearFELaplaceElementMatrix EL_Mat;
  LinFEMassMatrixProvider EL_Mass;
  // Let's assemble it locally with the two FEM
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, EL_Mat, A_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, EL_Mass, M_COO);
  // Set the boundary conditions
  auto mesh_p = dofh.Mesh();
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  auto bd_vertices = [&bd_flags,&dofh](unsigned int id) -> bool { 
    return bd_flags(dofh.Entity(id));
  };
  // And set the zeros in the matrices
  dropMatrixRowsColumns(bd_vertices,M_COO);
  dropMatrixRowsColumns(bd_vertices,A_COO);
  // Then compute the Kronecker product as in 9.1-e)
  A = A_COO.makeSparse();
  Eigen::SparseMatrix<double> M = M_COO.makeSparse();
  A_ = Eigen::kroneckerProduct(Coeff,A);
  M_ = Eigen::kroneckerProduct(I,M);
  //====================
}
/* SAM_LISTING_END_4 */

/* Implementation of Radau3MOLTimestepper member functions */
// The function discreteEvolutionOperator() returns the discretized evolution
// operator as obtained from the Runge-Kutta Radau IIA 2-stages method using the
// Butcher table as stored in the Radau3MOLTimestepper class
/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd Radau3MOLTimestepper::discreteEvolutionOperator(
    double time, double tau, const Eigen::VectorXd &mu) const {
  Eigen::VectorXd discrete_evolution_operator(dofh_.NumDofs());
  //====================
  Eigen::VectorXd RHS(2*dofh_.NumDofs());
  RHS << rhsVectorheatSource(dofh_,time+1/3*tau)-A*mu,rhsVectorheatSource(dofh_,time+tau)-A*mu;
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(tau*A_+M_);
  auto sol = solver.solve(RHS);
  discrete_evolution_operator = mu+tau*(3/4.*sol.topRows(dofh_.NumDofs())+1/4.*sol.bottomRows(dofh_.NumDofs())); // (9.2.7.8)
  //====================
  return discrete_evolution_operator;
}
/* SAM_LISTING_END_5 */

}  // namespace RadauThreeTimestepping
