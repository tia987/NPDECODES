/**
 * @file minimalgraphsurface.cc
 * @brief NPDE homework 5-3 Minimal Graph Surface code
 * @author R. Hiptmair & W. Tonnonw
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "minimalgraphsurface.h"

#include <lf/base/lf_assert.h>
#include <lf/fe/fe_tools.h>
#include <lf/fe/mesh_function_fe.h>

#include <Eigen/Core>


namespace MinimalGraphSurface {

/* SAM_LISTING_BEGIN_1 */
double computeGraphArea(
  std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
  const Eigen::VectorXd& mu_vec) {
    double area;
    //====================
    lf::fe::MeshFunctionGradFE<double,double> graduh(fes_p,mu_vec);
    auto norm = [&graduh](const lf::mesh::Entity& e, const Eigen::MatrixXd& refc) -> std::vector<double> {
        const std::vector<Eigen::VectorXd> gradvals{graduh(e, refc)};
        std::vector<double> ret(gradvals.size()); for (int i = 0; i < gradvals.size(); ++i) {
        ret[i] = std::sqrt(1.0 + gradvals[i].squaredNorm()); }
        return ret;
    };
    const unsigned order = 2;
    area = lf::fe::IntegrateMeshFunction(*fes_p->Mesh(),norm,order);
    //====================
    return area;
}
/* SAM_LISTING_END_1 */

// Implementation of the constructor
/* SAM_LISTING_BEGIN_2 */
CoeffTensorA::CoeffTensorA(
  std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
  const Eigen::VectorXd& mu) : graduh_(fes_p, mu) {
    //====================
    //====================
}
/* SAM_LISTING_END_2 */

// Implementation of evaluation operator for class CoeffTensorA
/* SAM_LISTING_BEGIN_3 */
std::vector<Eigen::Matrix2d> CoeffTensorA::operator()(
  const lf::mesh::Entity& e, const Eigen::MatrixXd& refc) const {
    // Number of points for which evaluation is requested
    const int nvals = refc.cols();
    // For returning values
    std::vector<Eigen::Matrix2d> Avals(nvals);
    //====================
    const std::vector<Eigen::VectorXd> gradvals{graduh_(e,refc)};
    for(unsigned i = 0; i < nvals; i++){
        const Eigen::Vector2d g{gradvals[i]};
        const double norms_g = g.squaredNorm();
        Avals[i] = 1. / std::sqrt(1+norms_g)*(Eigen::Matrix2d::Identity()-g*g.transpose()/(1.+norms_g));
    }
    //====================
    return Avals;
}
/* SAM_LISTING_END_3 */

// Implementation of constructor for CoeffScalarc
/* SAM_LISTING_BEGIN_5 */
CoeffScalarc::CoeffScalarc(
  std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
  const Eigen::VectorXd& mu) : graduh_(fes_p, mu) {
    //====================
    //====================
}
/* SAM_LISTING_END_5 */

// Implementation of evaluation operator for class CoeffScalarc
/* SAM_LISTING_BEGIN_4 */
std::vector<double> CoeffScalarc::operator()(
  const lf::mesh::Entity& e, const Eigen::MatrixXd& refc) const {
    // Number of points for which evaluation is requested
    const int nvals = refc.cols();
    // For returning values
    std::vector<double> cvals(nvals);
    //====================
    const std::vector<Eigen::VectorXd> gradvals{graduh_(e,refc)};
    for(unsigned i = 0; i < nvals; i++){
        const Eigen::Vector2d g{gradvals[i]};        
        cvals[i] = -1. / std::sqrt(1+gradvals[i].squaredNorm());
    }
    //====================
    return cvals;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_6 */
Eigen::VectorXd computeNewtonCorrection(
  std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fes_p,
  const Eigen::VectorXd& mu_vec) {
    // Obtain reference to the underlying finite element mesh
    const lf::mesh::Mesh& mesh{*fes_p->Mesh()};
    // The local-to-global index mapping
    const lf::assemble::DofHandler& dofh{fes_p->LocGlobMap()};
    // Get the number of degrees of freedom = dimension of FE space
    const lf::base::size_type N_dofs(dofh.NumDofs());
    LF_ASSERT_MSG(mu_vec.size() == N_dofs, "Vector length mismatch!");
    // Solution vector = return value
    Eigen::VectorXd sol_vec(N_dofs);
    //====================
    //auto c = CoeffScalarc(fes_p,mu_vec);
    // Assemble matrix A
    Eigen::VectorXd phi{N_dofs};
    lf::assemble::COOMatrix<double> A(N_dofs,N_dofs); 
    auto mf_A = CoeffTensorA(fes_p,mu_vec);
    lf::mesh::utils::MeshFunctionConstant<double> mf_zero(0.0);
    lf::uscalfe::ReactionDiffusionElementMatrixProvider mf_eval(fes_p,std::move(mf_A),mf_zero);
    lf::assemble::AssembleMatrixLocally(0,dofh,dofh,mf_eval,A);

    // Assemble vector c through AssembleMatrixLocally
    CoeffScalarc mf_c(fes_p,mu_vec);
    lf::assemble::COOMatrix<double> T(N_dofs,N_dofs);
    lf::uscalfe::ReactionDiffusionElementMatrixProvider Tmat_builder(fes_p,std::move(mf_c),mf_zero);
    lf::assemble::AssembleMatrixLocally(0,dofh,dofh,Tmat_builder,T);
    phi = T.MatVecMult(1.0,mu_vec);

    // Set boundary conditions
    auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(fes_p->Mesh(),2);
    auto selector = [&bd_flags, &dofh](lf::assemble::glb_idx_t gdof_idx) -> std::pair<bool,double>{
        const lf::mesh::Entity& node{dofh.Entity(gdof_idx)};
        return (bd_flags(node) ? std::make_pair(true,0.) : std::make_pair(false,0.));
    };
    lf::assemble::FixFlaggedSolutionComponents<double>(selector,A,phi);

    // compute Galerking inverse with sparse matrix solver
    Eigen::SparseMatrix<double> A_crs = A.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_crs);
    sol_vec = solver.solve(phi);
    //====================
    return sol_vec;
}
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_7 */
void graphMinSurfVis(std::string meshfile, std::string vtkfile) {
    // Read mesh for unit square from file
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    const lf::io::GmshReader reader(std::move(mesh_factory), meshfile.c_str());
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
    // Finite element space
    auto fe_space_p =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Solve non-linear BVP by means of Newton's method
    std::vector<Eigen::VectorXd> iterates{};
    int itcnt = 0;
    Eigen::VectorXd mu = graphMinimalSurface(
        //====================
        
        //====================
        fe_space_p, [](Eigen::Vector2d x) { return 0.0; }, 0.0, 0.0, 0
    );
    // Tabulate progress of iteration
    unsigned int N = mu.size();
    for (int i = 0; i < itcnt; ++i) {
      double mu_norm = iterates[i].norm() / std::sqrt(N);
      std::cout << "k = " << i << ": |mu(" << i << ")| = " << mu_norm;
      if (i > 0) {
        std::cout << ", |correction| = "
                  << (iterates[i] - iterates[i - 1]).norm() / std::sqrt(N);
      }
      std::cout << std::endl;
    }
    // Output solution for visualization
    //====================
    lf::io::VtkWriter vtk_writer (mesh_p, vtkfile ) ;
    // Create a MeshFunction encoding the FE solution
    lf::fe::MeshFunctionFE<double, double> mf_uh(fe_space_p, mu); // Write nodal data to file
    vtk_writer.WritePointData( "mPoint" , mf_uh) ;
    std::cout << "\n The finite element solution was written to:" << std::endl; std::cout << ">> " << vtkfile << std::endl;
    // Supplement: output for visualization with MATLAB
    //PrepTriMesh::prepTriMesh(fe_space_p,mu," minsurf ");
    //====================
}
/* SAM_LISTING_END_7 */

}  // namespace MinimalGraphSurface
