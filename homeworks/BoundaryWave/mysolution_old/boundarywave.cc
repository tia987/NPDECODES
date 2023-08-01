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
    lf ::mesh::utils::CodimMeshDataSet<bool> bd_flags{
        lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)
    };
    //lf::uscalfe::MassEdgeMatrixProvider m_builder(lf::assemble::COOMatrix<double>, mesh_p, bd_flags);
    auto edges_predicate = [&bd_flags](const lf::mesh::Entity &edge) -> bool {
        return bd_flags(edge);
    };
    auto one_fun = [](Eigen::Vector2d x) -> double {
        return 1.0;
    };
    auto eta = lf::mesh::utils::MeshFunctionGlobal(one_fun);
    lf::uscalfe::MassEdgeMatrixProvider<double,decltype(eta),decltype(edges_predicate)>
        edgemat_builder(fe_space_p, eta, edges_predicate);
    lf::assemble::AssembleMatrixLocally(1, dofh, dofh, edgemat_builder, M);
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
    auto one_fun = [](Eigen::Vector2d x) -> double {
        return 1.0+x.norm()*x.norm();
    };

    auto zero_fun = [](Eigen::Vector2d x) -> double {
        return 0.;
    };

    auto alpha = lf::mesh::utils::MeshFunctionGlobal(one_fun);
    auto gamma = lf::mesh::utils::MeshFunctionGlobal(zero_fun);

    lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(alpha), decltype(gamma)> 
        mat_builder(fe_space_p, alpha, gamma);
    
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, mat_builder, A);


    /*lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
        lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)
    };
    //lf::uscalfe::MassEdgeMatrixProvider m_builder(lf::assemble::COOMatrix<double>, mesh_p, bd_flags);
    auto edges_predicate = [&bd_flags](const lf::mesh::Entity &edge) −> bool { return bd_flags (edge);};
    auto eta = lf ::mesh:: utils ::MeshFunctionGlobal( [](Eigen::Vector2d x) -> double { return 1.0; });
    lf::uscalfe::MassEdgeMatrixProvider<double,decltype(eta),decltype(edges_predicate)>
        edgemat_builder(fe_space_p, eta, edges_predicate);
    A = lf::assemble::AssembleMatrixLocally(1, dofh, dofh, edgemat_builder, m_builder);

    
    Eigen::Sparse<double> A_crs = A.makeSparse();*/
    //====================
    return A;
};

/* SAM_LISTING_END_2 */

}  // namespace BoundaryWave
