/** @file
 * @brief NPDE ZienkiewiczZhuEstimator
 * @author Erick Schulz
 * @date 25/07/2019
 * @copyright Developed at ETH Zurich
 */

#include "zienkiewiczzhuestimator.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
// Eigen includes
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

namespace ZienkiewiczZhuEstimator {

/* Implementing member function Eval of class VectorProjectionMatrixProvider*/
/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd VectorProjectionMatrixProvider::Eval(
  const lf::mesh::Entity &entity) {
    Eigen::MatrixXd elMat_vec;  // element matrix to be returned
    // Throw error in case cell is not Tria nor Quad
    LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria() ||
                      entity.RefEl() == lf::base::RefEl::kQuad(),
                  "Unsupported cell type " << entity.RefEl());

    if (entity.RefEl() == lf::base::RefEl::kTria()) {
        elMat_vec = Eigen::MatrixXd::Zero(6, 6);
        // For TRIANGULAR CELLS
        // Compute the area of the triangle cell
        const double area = lf::geometry::Volume(*(entity.Geometry()));
        //====================
        bool chessboard = true;
        for(unsigned i = 0; i < 6; i++){
            for(unsigned j = 0; j < 6; j++){
                if(chessboard == true){
                    elMat_vec(i,j) = 1;
                    chessboard = false;
                } else{
                    chessboard = true;
                }
                if(i == j){
                   elMat_vec(i,j) += 1;
                }
            }
        }
        elMat_vec = elMat_vec*area/12.;
        //====================
    } else {
        // for QUADRILATERAL CELLS
        elMat_vec = Eigen::MatrixXd::Zero(8, 8);
        Eigen::MatrixXd elMat_scal =
            Eigen::MatrixXd::Zero(4, 4);  // element matrix for scalar FEM
        //====================
        const double area = lf::geometry::Volume(*(entity.Geometry()));
        bool chessboard = true;
        for(unsigned i = 0; i < 4; i++){
            for(unsigned j = 0; j < 4; j++){
                if(chessboard == true){
                    elMat_vec(i,j) = 1;
                    chessboard = false;
                } else{
                    chessboard = true;
                }
            }
        }
        elMat_vec.block(0,0,4,4) = elMat_scal;
        elMat_vec.block(4,0,4,4) = elMat_scal;
        elMat_vec.block(0,4,4,4) = elMat_scal;
        elMat_vec.block(4,4,4,4) = elMat_scal;
        for(unsigned i = 0; i < 8; i++){
            elMat_vec(i,i) += 1;
        }
        elMat_vec = elMat_vec*area/12.;
        //====================
    }
    return elMat_vec;  // return the local mass element matrix
}  //
/* SAM_LISTING_END_1 */

/* Implementing member function Eval of class GradientProjectionVectorProvider*/
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd GradientProjectionVectorProvider::Eval(
  const lf::mesh::Entity &entity) {
    Eigen::VectorXd elVec(6);  // for returning the element vector
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh{_fe_space_p->LocGlobMap()};
    // Obtain global indices of the vertices of the triangle entity
    auto dof_idx_vec = dofh.GlobalDofIndices(entity);
    LF_ASSERT_MSG(dofh.NumLocalDofs(entity) == 3,
                  "Too many global indices were returned for a triangle entity!");

    // Obtain the gradients of the barycentric coordinate functions
    Eigen::Matrix<double, 2, 3> elgrad_Mat = gradbarycoordinates(entity);
    // Compute the local constant gradient of the finite element solution
    Eigen::Vector2d grad_vec(0.0, 0.0);
    
    //====================
    ;
    //====================
    
    // Assemble local element vector
    // Compute the area of the triangle cell
    const double area = lf::geometry::Volume(*(entity.Geometry()));
    
    //====================
    ;
    //====================
    
    return elVec;
}  // GradientProjectionVectorProvider::Eval
/* SAM_LISTING_END_2 */

Eigen::Matrix<double, 2, 3> gradbarycoordinates(
    const lf::mesh::Entity &entity) {
  LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << entity.RefEl());

  // Get vertices of the triangle
  auto endpoints = lf::geometry::Corners(*(entity.Geometry()));

  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = endpoints.transpose();

  return X.inverse().block<2, 3>(1, 0);
}  // gradbarycoordinates

/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd computeLumpedProjection(
  const lf::assemble::DofHandler &scal_dofh, 
  const Eigen::VectorXd &mu,
  const lf::assemble::DofHandler &vec_dofh) {
    // Obtain shared_ptr to mesh
    std::shared_ptr<const lf::mesh::Mesh> mesh_p = scal_dofh.Mesh();
    // Dimension of vector-valued finite element space
    const lf::uscalfe::size_type N_vec_dofs(vec_dofh.NumDofs());
    // Initialize vector FE basis expansion coefficient vector with zeros
    Eigen::VectorXd proj_vec(N_vec_dofs);
    proj_vec.setZero();
    // Initialize temporary helper nodal DataSet (codim 2)
    auto nodal_sum_of_areas =
        lf::mesh::utils::CodimMeshDataSet<double>(mesh_p, 2, 0.0);

    // Loop over the triangular cells of the mesh in the spirit of
    // cell oriented assembly
    for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
        LF_VERIFY_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                      "Unsupported cell type " << cell->RefEl());
        // Obtain global scalar-FE indices of the vertices
        const auto scal_dof_idx_vec = scal_dofh.GlobalDofIndices(*cell);
        // Obtain the gradients of the barycentric coordinate functions
        const Eigen::Matrix<double, 2, 3> elgrad_Mat = gradbarycoordinates(*cell);
        // Obtain area of the triangular cell
        const double area = lf::geometry::Volume(*(cell->Geometry()));

        // Compute the gradient of the passed coefficient vector
        //====================
        //====================

    }

    // Scaling of components of vector of dofs
    for (const lf::mesh::Entity *node : mesh_p->Entities(2)) {
        LF_VERIFY_MSG(node->RefEl() == lf::base::RefEl::kPoint(),
                      "Expected kPoint type!" << node->RefEl());
        //====================
        //====================
    }
    return proj_vec;
};  // computeLumpedProjection

/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
double computeL2Deviation(const lf::assemble::DofHandler &scal_dofh,
                          const Eigen::VectorXd &eta,
                          const lf::assemble::DofHandler &vec_dofh,
                          const Eigen::VectorXd &gamma) {
  double deviation_norm_value = 0.0;  // For retrurning the result
  // Obtain shared_ptr to mesh
  auto mesh_p = scal_dofh.Mesh();
  // Cell-oriented computation of deviation norm (squared)
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
    LF_VERIFY_MSG(cell->RefEl() == lf::base::RefEl::kTria(),
                  "Unsupported cell type " << cell->RefEl());
    // Obtain area of the triangular cell
    const double area = lf::geometry::Volume(*(cell->Geometry()));

    // Obtain global scalar-FE indices of the vertices
    auto scal_dof_idx_vec = scal_dofh.GlobalDofIndices(*cell);
    // Obtain the gradients of the barycentric coordinates functions
    Eigen::Matrix<double, 2, 3> elgrad_Mat = gradbarycoordinates(*cell);
// Compute the gradient of the passed coefficient vector eta
//====================
// Your code goes here
//====================
  }
  return std::sqrt(deviation_norm_value);
};  // computeL2Deviation
/* SAM_LISTING_END_4 */

Eigen::VectorXd solveBVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p) {
  Eigen::VectorXd discrete_solution;

  // TOOLS AND DATA
  // Pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  // Obtain specification for shape functions on edges
  const lf::fe::ScalarReferenceFiniteElement<double> *rsf_edge_p =
      fe_space_p->ShapeFunctionLayout(lf::base::RefEl::kSegment());

  // Dirichlet data
  auto mf_g = lf::mesh::utils::MeshFunctionGlobal(
      [](coord_t x) -> double { return 0.0; });
  // Right-hand side source function f
  auto mf_f = lf::mesh::utils::MeshFunctionGlobal([](coord_t x) -> double {
    return sin(M_PI * x[0]) * sin(2 * M_PI * x[1]);
  });

  // I : ASSEMBLY
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Right hand side vector, must be initialized with 0!
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  // I.i : Computing volume matrix for negative Laplace operator
  // Initialize object taking care of local mass (volume) computations.
  lf::uscalfe::LinearFELaplaceElementMatrix elmat_builder{};
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix A.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // I.ii : Computing right-hand side vector
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder(fe_space_p, mf_f);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // I.iii : Imposing essential boundary conditions
  // Obtain an array of boolean flags for the edges of the mesh, 'true'
  // indicates that the edge lies on the boundary (codim = 1)
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // Inspired by the example in the documentation of
  // InitEssentialConditionFromFunction()
  // https://craffael.github.io/lehrfempp/namespacelf_1_1uscalfe.html#a5afbd94919f0382cf3fb200c452797ac
  // Creating a predicate that will guarantee that the computations are carried
  // only on the exterior boundary edges of the mesh using the boundary flags
  auto edges_predicate_Dirichlet =
      [&bd_flags](const lf::mesh::Entity &edge) -> bool {
    return bd_flags(edge);
  };
  // Determine the fixed dofs on the boundary and their values
  // Alternative: See lecturedemoDirichlet() in
  // https://github.com/craffael/lehrfempp/blob/master/examples/lecturedemos/lecturedemoassemble.cc
  auto edges_flag_values_Dirichlet{lf::fe::InitEssentialConditionFromFunction(
      *fe_space_p, edges_predicate_Dirichlet, mf_g)};
  // Eliminate Dirichlet dofs from the linear system
  lf::assemble::FixFlaggedSolutionCompAlt<double>(
      [&edges_flag_values_Dirichlet](lf::assemble::glb_idx_t gdof_idx) {
        return edges_flag_values_Dirichlet[gdof_idx];
      },
      A, phi);

  // Assembly completed! Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_sparse = A.makeSparse();

  // II : SOLVING  THE LINEAR SYSTEM
  // II.i : Setting up Eigen's sparse direct elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_sparse);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  // II.ii : Solving
  discrete_solution = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

  return discrete_solution;
};  // solveBVP

Eigen::VectorXd solveGradVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p,
    const Eigen::VectorXd &mu, const lf::assemble::DofHandler &vec_dofh) {
  Eigen::VectorXd approx_grad;

  // TOOLS AND DATA
  // Pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  // Dimension of finite vector elements space
  const lf::uscalfe::size_type N_vec_dofs(vec_dofh.NumDofs());
  // Matrix in triplet format holding the Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> M_COO(N_vec_dofs, N_vec_dofs);
  // Right-hand side vector has to be set to zero initially
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_vec_dofs);
  phi.setZero();

  // Initialize classes containing the information required for the
  // local computations of the Galerkin matrix and the load vector
  VectorProjectionMatrixProvider elMat_builder;
  GradientProjectionVectorProvider elVec_builder(fe_space_p, mu);
  // Compute the Galerkin matrix and load vector
  lf::assemble::AssembleMatrixLocally(0, vec_dofh, vec_dofh, elMat_builder,
                                      M_COO);
  lf::assemble::AssembleVectorLocally(0, vec_dofh, elVec_builder, phi);

  // Solve the linear problem
  Eigen::SparseMatrix<double> M = M_COO.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(M);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  approx_grad = solver.solve(phi);

  return approx_grad;
};  // solveGradVP

double getMeshSize(const std::shared_ptr<const lf::mesh::Mesh> &mesh_p) {
  double mesh_size = 0.0;

  // Find maximal edge length
  double edge_length;
  for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
    // Compute the length of the edge
    auto endpoints = lf::geometry::Corners(*(edge->Geometry()));
    edge_length = (endpoints.col(0) - endpoints.col(1)).norm();
    if (mesh_size < edge_length) {
      mesh_size = edge_length;
    }
  }

  return mesh_size;
};  // getMeshSize

void progress_bar::write(double fraction) {
  // clamp fraction to valid range [0,1]
  if (fraction < 0)
    fraction = 0;
  else if (fraction > 1)
    fraction = 1;

  auto width = bar_width - message.size();
  auto offset = bar_width - static_cast<unsigned>(width * fraction);

  std::string sign = "%";
  os << '\r' << message;
  os.write(full_bar.data() + offset, width);
  os << " [completed: " << std::setw(3) << static_cast<int>(100 * fraction)
     << sign + " of meshes]" << std::flush;
};

}  // namespace ZienkiewiczZhuEstimator
