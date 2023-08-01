/** @file
 * @brief NPDE BoundaryWave
 * @author Erick Schulz
 * @date 24/07/2019
 * @copyright Developed at ETH Zurich
 */

#include "boundarywave.h"

namespace BoundaryWave {

/* SAM_LISTING_BEGIN_1 */
lf::assemble::COOMatrix<double> buildM(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p) {
    // I. TOOLS AND DATA
    // Pointer to current fe_space and mesh
    std::shared_ptr<const lf::mesh::Mesh> mesh_p(fe_space_p->Mesh());
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
    // Dimension of finite element space
    const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

    // II : ASSEMBLY
    // Matrix in triplet format holding Galerkin matrix, zero initially.
    lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);
    //====================
    // Set up gamma function
    auto gamma = lf::mesh::utils::MeshFunctionGlobal([](Eigen::VectorXd){return 1.0;});

    // Set Selector
    auto bd_flags {lf::mesh::utils::flagEntitiesOnBoundary(mesh_p,1)};

    auto edge_sel = [bd_flags](const lf::mesh::Entity &edge){
        return bd_flags(edge);
    };

    // Build Matrix M
    lf::uscalfe::MassEdgeMatrixProvider<double, decltype(gamma), decltype(edge_sel)> M_builder(fe_space_p,gamma,edge_sel);
    lf::assemble::AssembleMatrixLocally(1,dofh,dofh,M_builder,M);
    //====================
    return M;
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
lf::assemble::COOMatrix<double> buildA(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p) {
    // I. TOOLS AND DATA
    // Pointer to current fe_space and mesh
    std::shared_ptr<const lf::mesh::Mesh> mesh_p(fe_space_p->Mesh());
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
    // Dimension of finite element space
    const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

    // II : ASSEMBLY
    // Matrix in triplet format holding Galerkin matrix, zero initially.
    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);

    //====================
    // Set up alpha function
    auto alpha = lf::mesh::utils::MeshFunctionGlobal([](Eigen::VectorXd x){return 1+x.dot(x);});

    // Set gamma  function
    auto zero = lf::mesh::utils::MeshFunctionGlobal([](Eigen::VectorXd x){return 0.;});

    // Build Matrix M
    lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(alpha), decltype(zero)> A_builder(fe_space_p,alpha,zero);
    lf::assemble::AssembleMatrixLocally(0,dofh,dofh,A_builder,A);
    //====================

    return A;
};

/* SAM_LISTING_END_2 */

}  // namespace BoundaryWave
