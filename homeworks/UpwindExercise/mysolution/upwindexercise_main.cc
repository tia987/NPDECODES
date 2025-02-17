/**
* @file upwindquadrature_main.cc
* @brief NPDE homework template main
* @author Philippe Peter
* @date June 2020
* @copyright Developed at SAM, ETH Zurich
*/
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <fstream>
#include <memory>

#include "../../../lecturecodes/ConvectionDiffusion/convection_emp.h"
#include "upwindexercise.h"

int main() {
   // PARAMETERS
   // mesh specification (number of cells in both sides of the tensor-product
   // triangular mesh)
   int M = 49;

   // coefficient functions:
   // Dirichlet functor
   const auto g = [](const Eigen::Vector2d &x) {
       return x(1) == 0 ? 0.5 - std::abs(x(0) - 0.5) : 0.0;
   };
   lf::mesh::utils::MeshFunctionGlobal mf_g{g};

   // velocity field
   const auto v = [](const Eigen::Vector2d &x) {
       return (Eigen::Vector2d() << -x(1), x(0)).finished();
   };

   // diffusion coefficient
   const double eps = 1e-4;
   lf::mesh::utils::MeshFunctionConstant mf_eps{eps};

   // MESH CONSTRUCTION
   // construct a triangular tensor product mesh on the unit square
   std::unique_ptr<lf::mesh::MeshFactory> mesh_factory_ptr =
       std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
   lf::mesh::utils::TPTriagMeshBuilder builder(std::move(mesh_factory_ptr));
   builder.setBottomLeftCorner(Eigen::Vector2d{0.0, 0.0})
       .setTopRightCorner(Eigen::Vector2d{1.0, 1.0})
       .setNumXCells(M)
       .setNumYCells(M);
   std::shared_ptr<lf::mesh::Mesh> mesh_p = builder.Build();

   // DOF HANDLER & FINITE ELEMENT SPACE
   // Construct dofhanlder for linear finite elements on the mesh.
   // TODO: 1.1 get dofhandler
   auto fe_space =
       std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
   const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

   // PREPARING DATA TO IMPOSE DIRICHLET CONDITIONS
   // Create a dataset of boolean flags indicating edges on the boundary of the
   // mesh
   auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

   // Fetch flags and values for degrees of freedom located on Dirichlet
   // boundary.
   // TODO: 1.2 make an edge flagger (see InitEssentialConditionFromFunction in lehrfem++ Documentation)
   auto ess_bdc_flags_values =  InitEssentialConditionFromFunction(*fe_space, lf::base::PredicateTrue{}, mf_boundary);

   //============================================================================
   // SOLVE LAPLACIAN WITH NON-HOMOGENEOUS DIRICHLET BC (STANDARD: UNSTABLE)
   //============================================================================
   // Matrix in triplet format holding Galerkin matrix, zero initially.
   lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());

   // ASSEMBLE GALERKIN MATRIX
   // First the part corresponding to the laplacian
   lf::uscalfe::ReactionDiffusionElementMatrixProvider laplacian_provider(
       fe_space, mf_eps, lf::mesh::utils::MeshFunctionConstant(0.0));

   // TODO: 1.3 assemble the laplace Galerkin matrix part
   lf::assemble::AssembleMatrixLocally(0, dofh, dofh, laplacian_provider, A);

   // Next part corresponding to the convection term:
   ConvectionDiffusion::ConvectionElementMatrixProvider convection_provider(v);

   //TODO: 1.4 assemble the convection Galerkin matrix part   
   lf::assemble::AssembleMatrixLocally(0, dofh, dofh, convection_provider, A);

   // RIGHT-HAND SIDE VECTOR
   Eigen::VectorXd phi(dofh.NumDofs());
   phi.setZero();

   // IMPOSE DIRICHLET CONDITIONS:
   // Eliminate Dirichlet dofs from linear system
   // TODO: 1.5 pass ess_bdc_flags_values to FixFlaggedSolutionComponents
   // make sure it has a call operator (Hint: make a lambda)
   lf::assemble::FixFlaggedSolutionComponents<double>(
       [&](int gdof_idx) {
           return ess_bdc_flags_values[gdof_idx];
       },
       A, phi);

   // SOLVE LINEAR SYSTEM
   // TODO: 1.6 convert A into an Eigen::SparseMatrix and solve the linear system of equations
   Eigen::SparseMatrix A_crs = A.makeSparse();
   Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
   solver.compute(A_crs);
   Eigen::VectorXd sol_vec = solver.solve(phi);

   // output solution
   std::ofstream solution_file_unstable(CURRENT_SOURCE_DIR "/upwind_quadrature_solution_unstable.txt");
   if (solution_file_unstable.is_open()) {
       solution_file_unstable << sol_vec << std::endl;
       std::cout << "Unstable solution saved!" << std::endl;
   }
   //============================================================================
   // SOLVE LAPLACIAN WITH NON-HOMOGENEOUS DIRICHLET BC (UPWIND: STABLE)
   //============================================================================
   /* SAM_LISTING_BEGIN_7 */
   // Matrix in triplet format holding Galerkin matrix, zero initially.
   // TODO: 2 now the whole thing again but with upwind quadrature
   lf::assemble::COOMatrix<double> A_stable(dofh.NumDofs(), dofh.NumDofs());

   // ASSEMBLE GALERKIN MATRIX
   // First the part corresponding to the laplacian, computed using standard
   // Galerkin approach
   // TODO: 2.1 assemble the laplace Galerkin matrix (see 1.3)
   lf::assemble::AssembleMatrixLocally(0, dofh, dofh, laplacian_provider,
                                       A_stable);

   // Next part corresponding to the convection term, computed using upwind
   // quadrature:
   // TODO: 2.2 go into the files upwindexercise.cc/.h and complete the code there
   UpwindQuadrature::UpwindConvectionElementMatrixProvider
       convection_provider_stable(v, UpwindQuadrature::initializeMasses(mesh_p));
   // TODO: 2.3 assemble the advection Gelerkin matrix (see 1.4)
   lf::assemble::AssembleMatrixLocally(0, dofh, dofh, convection_provider_stable,
                                       A_stable);   

   // RIGHT-HAND SIDE VECTOR
   Eigen::VectorXd phi_stable(dofh.NumDofs());
   phi_stable.setZero();

   // IMPOSE DIRICHLET CONDITIONS:
   // Eliminate Dirichlet dofs from linear system
   // TODO: 2.4 pass ess_bdc_flags_values to FixFlaggedSolutionComponents (see 1.5)
   lf::assemble::FixFlaggedSolutionComponents<double>(
       [&](int gdof_idx) {
           return ess_bdc_flags_values[gdof_idx];
       },
       A_stable, phi_stable);

   // SOLVE LINEAR SYSTEM
   // TODO: 2.5 convert A into an Eigen::SparseMatrix and solve the linear system of equations (see 1.6)
   Eigen::SparseMatrix A_stable_crs = A.makeSparse();
   Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_stable;
   solver_stable.compute(A_stable_crs);
   Eigen::VectorXd sol_vec_stable = solver_stable.solve(phi_stable);

   std::ofstream solution_file(CURRENT_SOURCE_DIR
                               "/upwind_quadrature_solution_stable.txt");
   if (solution_file.is_open()) {
       solution_file << sol_vec_stable << std::endl;
       std::cout << "Upwind solution saved!" << std::endl;
   }
   /* SAM_LISTING_END_7 */
   return 0;
}
