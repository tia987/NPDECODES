/** @file
 * @brief NPDE OutputImpedanceBVP
 * @author Erick Schulz
 * @date 12/07/2019
 * @copyright Developed at ETH Zurich
 */

#include "outputimpedancebvp.h"

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>

namespace OutputImpedanceBVP {

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd solveImpedanceBVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p,
    Eigen::Vector2d g) {
    // Related implementations:
    // Homework problem ErrorEstimatesForTraces:
    // https://gitlab.math.ethz.ch/ralfh/npdecodes/tree/master/homeworks/ErrorEstimatesForTraces

    // Pointer to current mesh
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
    // Dimension of finite element space
    const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
    // Obtain specification for shape functions on edges
    const lf::fe::ScalarReferenceFiniteElement<double> *rsf_edge_p =
        fe_space_p->ShapeFunctionLayout(lf::base::RefEl::kSegment());

    Eigen::VectorXd discrete_solution(N_dofs);

    // I : ASSEMBLY
    // Matrix in triplet format holding Galerkin matrix, zero initially.
    lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
    // Right hand side vector, must be initialized with 0!
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
    phi.setZero();

    // I.i : Computing volume matrix for negative Laplace operator
    //====================
    auto mf_one = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d){return 1.;});
    auto mf_zero = [](Eigen::Vector2d){return 0.;};
    lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(mf_one), decltype(mf_one)> built_mat(fe_space_p, mf_one, mf_one);
    lf::assemble::AssembleMatrixLocally(0,dofh,dofh,built_mat,A_COO);
    //====================
    /* SAM_LISTING_END_1 */

    /* SAM_LISTING_BEGIN_2 */
    // I.ii : Computing mass edge matrix resulting from Robin B.C.
    // Obtain an array of boolean flags for the edges of the mesh, 'true'
    // indicates that the edge lies on the boundary
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

    //====================
    auto edges_predicate_RobinBC = [&bd_flags](const lf::mesh::Entity &edge) -> bool {
        if (bd_flags(edge)) {
            auto endpoints = lf ::geometry::Corners(*(edge.Geometry()));
            if (endpoints(0, 0) < 0.01 || 0.9 < endpoints(0, 0) || endpoints(1, 0) < 0.1 || 0.9 < endpoints(1, 0)) {
                return false;
            }
            return true;
        }
        return false;
    };
    //====================
    /* SAM_LISTING_END_2 */

    /* SAM_LISTING_BEGIN_9 */
    // I.iii : Computing right-hand side vector
    // Right-hand side source function f
    auto mf_f = lf::mesh::utils::MeshFunctionGlobal(
        [](Eigen::Vector2d x) -> double { return 0.0; });
    lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
        elvec_builder(fe_space_p, mf_f);
    // Invoke assembly on cells (codim == 0)
    AssembleVectorLocally(0, dofh, elvec_builder, phi);

    // I.iv : Imposing essential boundary conditions
    // Dirichlet data
    auto mf_g = lf::mesh::utils::MeshFunctionGlobal(
        [&g](Eigen::Vector2d x) -> double { return g.dot(x); });

    //====================
    auto edges_predicate_Dirichlet = [&bd_flags](const lf::mesh::Entity &edge) -> bool {
        if (bd_flags(edge)) {
            auto endpoints = lf ::geometry::Corners(*(edge.Geometry()));
            if (endpoints(0, 0) < 0.01 || 0.9 < endpoints(0, 0) || endpoints(1, 0) < 0.1 || 0.9 < endpoints(1, 0)) {
                return true;
            }
        }
        return false;
    };

    auto edges_flag_values_Dirichlet = lf::fe::InitEssentialConditionFromFunction(fe_space_p, edges_predicate_Dirichlet , mf_g); // Eliminate Dirichlet dofs from the linear system
lf::assemble::FixFlaggedSolutionCompAlt<double>( [&edges_flag_values_Dirichlet](lf::assemble::glb_idx_t gdof_idx) {
return edges_flag_values_Dirichlet [ gdof_idx ]; },A, phi);
    //====================

    // Assembly completed! Convert COO matrix A into CRS format using Eigen's
    // internal conversion routines.

    // II : SOLVING  THE LINEAR SYSTEM
    //====================
    Eigen::SparseMatrix<double> A = A_COO.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.factorize(A);
    discrete_solution = solver.solve(phi);
    //====================

    return discrete_solution;
};
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_3 */
double computeBoundaryOutputFunctional(
    const Eigen::VectorXd eta,
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p,
    Eigen::Vector2d d) {
    double func_val = 0.0;
    // Pointer to current mesh
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
  
    // Obtain an array of boolean flags for the edges of the mesh, 'true'
    // indicates that the edge lies on the boundary
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  
    //====================
    // Your code goes here
    //====================
  
    // Computing value of the functional
    for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
        //====================
        // Your code goes here
        //====================
    }
    return func_val;
};
/* SAM_LISTING_END_3 */

}  // namespace OutputImpedanceBVP
