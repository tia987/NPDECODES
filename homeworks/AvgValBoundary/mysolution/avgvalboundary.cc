/**
 * @ file avgvalboundary.cc
 * @ brief NPDE homework AvgValBoundary code
 * @ author Simon Meierhans, edited by Oliver Rietmann
 * @ date 11.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "avgvalboundary.h"

#include <lf/assemble/assemble.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace AvgValBoundary {

/**
 * @brief computes H1 seminorm over the computational domain
 * @param dofh DofHandler of FEspace.
 *        u coefficient vector
 */
/* SAM_LISTING_BEGIN_1 */
double compH1seminorm(const lf::assemble::DofHandler &dofh,
                      const Eigen::VectorXd & u) {
  double result = 0.0;
  //====================
  // From 3.1.2 we can see that for computing the gradient (H1 norm) we only need
  // to have alpha = 1 ^ beta = 0 ^ gamma = 0
  auto a = [](Eigen::VectorXd)->double{return 1.0;};
  auto b = [](Eigen::VectorXd)->double{return 0.0;};
  auto g = [](Eigen::VectorXd)->double{return 0.0;};
  // After compute the matrix and apply the norm
  Eigen::SparseMatrix<double> A_crs = compGalerkinMatrix(dofh,a,b,g);
  result = std::sqrt(u.transpose()*A_crs*u);
  //====================
  return result;
}
/* SAM_LISTING_END_1 */

/**
 * @brief solves pde for some simple test problem with
 *        alpha = beta = gamma := 1.0 and load f := 1.0
 * @param dofh DofHandler of FEspace.
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solveTestProblem(const lf::assemble::DofHandler &dofh) {
  // Obtain Galerkin matrix for alpha = beta = gamma := 1.0
  auto const_one = [](Eigen::Vector2d x) -> double { return 1.0; };
  auto A = compGalerkinMatrix(dofh, const_one, const_one, const_one);

  // Set up load vector
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(dofh.Mesh());
  lf::mesh::utils::MeshFunctionConstant mf_identity{1.0};
  lf::uscalfe::ScalarLoadElementVectorProvider elvec_builder(fe_space,
                                                             mf_identity);
  Eigen::VectorXd phi(dofh.NumDofs());
  phi.setZero();
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // Solve system of linear equations
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  Eigen::VectorXd mu = solver.solve(phi);

  return mu;
}
/* SAM_LISTING_END_2 */

/** @brief generate sequence of nested triangular meshes with L+1 levels */
/* SAM_LISTING_BEGIN_3 */
std::shared_ptr<lf::refinement::MeshHierarchy> generateTestMeshSequence(
    unsigned int L) {
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  std::shared_ptr<lf::refinement::MeshHierarchy> meshes =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh, L);
  return meshes;
}
/* SAM_LISTING_END_3 */

/**
 * @brief Compute the boundary functional values for a sequence of L meshes
 * @param L The number of meshes in the sequence
 * @returns A vector of pairs containing the number of DOFs and value of the
 *	    boundary functional for each level
 */
/* SAM_LISTING_BEGIN_5 */
/*std::vector<std::pair<unsigned int, double>> approxBoundaryFunctionalValues(
    unsigned int L) {
  std::vector<std::pair<unsigned int, double>> result;
  //====================
  // First we generate the meshes from generateTestMeshSequence
  std::shared_ptr<lf::refinement::MeshHierarchy> gTMS = generateTestMeshSequence(L);
  // Afterwards we need to compute the values at F(u_l)
  for(unsigned i = 0; i < L; i++){
      double f;
      gTMS->getMesh(i)
      result.push_back(i,f);
    }
  }
  //====================
  return result;*/
std::vector<std::pair<unsigned int, double>> approxBoundaryFunctionalValues(
    unsigned int L) {
  std::vector<std::pair<unsigned int, double>> result;
  auto meshes = generateTestMeshSequence(L - 1);
  int num_meshes = meshes->NumLevels();
  for (int level = 0; level < num_meshes; ++level) {
    auto mesh_p = meshes->getMesh(level);

    // Set up global FE space; lowest order Lagrangian finite elements
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
    const lf::base::size_type N_dofs(dofh.NumDofs());

    // compute galerkin matrix with alpha = 1.0, gamma = 1.0, beta = 0.0
    auto const_one = [](Eigen::Vector2d x) -> double { return 1.0; };
    auto const_zero = [](Eigen::Vector2d x) -> double { return 0.0; };
    auto A = compGalerkinMatrix(dofh, const_one, const_one, const_zero);

    // compute load vector for f(x) = x.norm()
    auto f = [](Eigen::Vector2d x) -> double { return x.norm(); };
    lf::mesh::utils::MeshFunctionGlobal mf_f{f};
    lf::uscalfe::ScalarLoadElementVectorProvider elvec_builder(fe_space, mf_f);
    Eigen::VectorXd phi(N_dofs);
    phi.setZero();
    AssembleVectorLocally(0, dofh, elvec_builder, phi);

    // solve system of equations
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd u = solver.solve(phi);

    // set up weight function
    auto w = [](Eigen::Vector2d x) -> double { return 1.0; };
    double functional_value = compBoundaryFunctional(dofh, u, w);

    result.push_back({N_dofs, functional_value});
  }
  return result;

}
/* SAM_LISTING_END_5 */

}  // namespace AvgValBoundary
