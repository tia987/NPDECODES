/**
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss, edited Amélie Loher
 * @date   18.03.2019, 03.03.20
 * @copyright Developed at ETH Zurich
 */

#ifndef NUMPDE_SOLVE_CR_DIRICHLET_BVP_H
#define NUMPDE_SOLVE_CR_DIRICHLET_BVP_H

#include <lf/assemble/assemble.h>
#include <lf/uscalfe/uscalfe.h>

#include "crfespace.h"

namespace NonConformingCrouzeixRaviartFiniteElements {

template <typename GAMMA_COEFF, typename F_FUNCTOR>
Eigen::VectorXd solveCRDirichletBVP(std::shared_ptr<CRFeSpace> fe_space,
                                    GAMMA_COEFF &&gamma, F_FUNCTOR &&f) {
  Eigen::VectorXd sol;
// TODO: task 2-14.v)
  //====================
  // No need to correct
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // Dimension of finite element space`
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // identity mesh function for very simple problem
  lf::mesh::utils::MeshFunctionConstant mf_identity(1.0);

  auto zero = [](const Eigen::Vector2d&) -> double { return 1.; };
  lf::mesh::utils::MeshFunctionGlobal mf_zero{zero};

  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);

  // Obtain an object that computes the element matrix for the
  // volumne part of the bilinear form
  lf::uscalfe::ReactionDiffusionElementMatrixProvider elmat_builder(
      fe_space, mf_identity, mf_identity);

  // Invoke assembly on cells (co-dimension = 0)
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

   // Right-hand side vector; has to be set to zero initially
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  lf::mesh::utils::CodimMeshDataSet<bool> edges{ lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(),1)};
  lf::assemble::FixFlaggedSolutionComponents<double> ([&edges, &dofh](lf::assemble::glb_idx_t gdof_idx)->std::pair<bool,double>{
    const lf::mesh::Entity &edge{dofh.Entity(gdof_idx)};
    return {edges(edge),0.0};
  }, A, phi);

  Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  sol = solver.solve(phi);
  //====================
  return sol;
}

}  // namespace NonConformingCrouzeixRaviartFiniteElements

#endif  // NUMPDE_SOLVE_CR_DIRICHLET_BVP_H
