/**
 * @file heatevolutionanimation_main.cc
 * @brief NPDE homework HeatEvolutionAnimation code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <memory>

#include "heatevolutionanimation.h"

#define N 20
#define T 2.
#define steps 100

#define PI 3.1415926535897

std::shared_ptr<lf::mesh::Mesh> build_mesh() {
    // Obtain a triangular mesh of the unit square from the collection of
    // test meshes
    lf::mesh::utils::TPTriagMeshBuilder builder(
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
    // Set mesh parameters following the Builder pattern
    // Domain is the unit square
    builder.setBottomLeftCorner(Eigen::Vector2d{-1., -1.})
        .setTopRightCorner(Eigen::Vector2d{1, 1})
        .setNumXCells(N)
        .setNumYCells(N);
    return builder.Build();
}

int main() {
    const double tau = T / steps;

    // Source function
    auto source = [](const double t) {
        return [t](const Eigen::Vector2d &x) {
            Eigen::Vector2d c{sin(PI*t), cos(PI*t)};
            return (x - 0.5*c).norm() < 0.5 ? 1. : 0.;
        };
    };

    std::shared_ptr<lf::mesh::Mesh> mesh_p = build_mesh();
    const lf::mesh::Mesh &mesh{*mesh_p};

    // TODO 1. Setup of Linear Finite Element Simplicial Lagrangian Space.
    // TODO 1.1: Create finite element simplicial space from mesh_p
    auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Obtain local->global index mapping for current finite element space
    // TODO 1.2: Get the dofhandler from the fe_space
    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

    // Dimension of finite element space`
    const lf::base::size_type N_dofs(dofh.NumDofs());

    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
    lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);


    // TODO 2. Build lhs matrix
    // TODO 2.1 Build Stiffness matrix A.
    auto alpha_A = [](const Eigen::Vector2d&){        
        return 1.;
    };
    auto gamma_A = [](const Eigen::Vector2d&){
        return 0.;
    };
    
    lf::mesh::utils::MeshFunctionGlobal mf_alpha_A{alpha_A};
    lf::mesh::utils::MeshFunctionGlobal mf_gamma_A{gamma_A};

    //  Use the ReactionDiffusionElementMatrixProvider and the
    //  AssembleMatrixLocally functions to build the matrix.

    lf::uscalfe::ReactionDiffusionElementMatrixProvider elmat_A{fe_space,mf_alpha_A, mf_gamma_A};
    lf::assemble::AssembleMatrixLocally(0,dofh,dofh,elmat_A,A);
    Eigen::SparseMatrix<double> A_crs = A.makeSparse();

    // TODO 2.2 Build Mass Matrix M.
    auto alpha_M = [](const Eigen::Vector2d&){
        return 0.;
    };
    auto gamma_M = [](const Eigen::Vector2d&){
        return 1.;
    };
    lf::mesh::utils::MeshFunctionGlobal mf_alpha_M{alpha_M};
    lf::mesh::utils::MeshFunctionGlobal mf_gamma_M{gamma_M};

    //  Use the ReactionDiffusionElementMatrixProvider and the
    //  AssembleMatrixLocally functions to build the matrix.

    lf::uscalfe::ReactionDiffusionElementMatrixProvider elmat_M{fe_space,mf_alpha_M, mf_gamma_M};
    lf::assemble::AssembleMatrixLocally(0,dofh,dofh,elmat_M,M);
    Eigen::SparseMatrix<double> M_crs = M.makeSparse();

    // TODO 2.3 Create the final lhs matrix
    Eigen::SparseMatrix<double> lhs = M_crs+tau*A_crs;


    // TODO 4. Boundary Conditions
    // TODO 4.1 Create a CodimMeshDataSet<bool> from the mesh pointer for the
    //  vertices that belong to the boundary.
    //  Hint: flagEntitiesOnBoundary
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
    // TODO 4.2 Create a predicate (lambda that returns a boolean) that takes
    //  a vertex (of type lf::mesh::Entity) as argument and returns if it belongs
    //  to the boundary or not
    auto bd_predicate = [bd_flags](auto vertex){
         if (bd_flags(*vertex)) {
             return true;
         }
         else {
             return false;
         }
    };


    // TODO 4.3: Uncomment the following section. You also have to store in
    //  the variable fill_me the correct value.

    // Ensure BCs: note: this is not efficient, could be better done if we
    // could have the matrix in RowMajor Format, but Lehrfem++ only allow for
    // ColumnMajor as this is the default of Eigen

    const double fill_me = 1.;
    for (const lf::mesh::Entity* vertex : mesh.Entities(2)) {
        if(bd_predicate(*vertex)) {
            int global_idx = dofh.GlobalDofIndices(*vertex)[0];
            lhs.row(global_idx) *= 0.;
            lhs.coeffRef(global_idx, global_idx) = fill_me;
        }
    }

    // Initial Conditions
    Eigen::VectorXd mu = Eigen::VectorXd::Zero(N_dofs);

    auto ICs = [](const Eigen::Vector2d &x) {
        if(x.norm() < 1./2.) return 1.;
        else return 0.;
    };
    lf::mesh::utils::MeshFunctionGlobal mf_ICs{ICs};

    // Evaluate ICs on vertices of the mesh
    const Eigen::Matrix<double, 0, 1> dummy;
    for (const lf::mesh::Entity* vertex : mesh.Entities(2)) {
        mu[dofh.GlobalDofIndices(*vertex)[0]] = mf_ICs(*vertex, dummy)[0];
    }


    Eigen::MatrixXd solution(steps, N_dofs);
    solution.row(0) = mu;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(lhs);

    // TODO 3. Timestepping
    double t = 0.;
    for(int i = 1; i < steps; ++i) {
        t += tau;

        Eigen::VectorXd phi = Eigen::VectorXd::Zero(N_dofs);
        lf::mesh::utils::MeshFunctionGlobal mf_source{source(t)};

        // TODO 3.3: Evaluate mf_source on every vertex and store it.
        //  in the corresponding global_idx entry
        //  Use the dofh.GlobalDofIndices(*vertex) function
        //  Hint: Look at the documentation of MeshFunctionGlobal and how
        //  ICs initialization looks like.


        // TODO 3.1 Create the rhs vector
       Eigen::VectorXd rhs = M_crs*mu+tau*phi;

        // TODO 4.4: Uncomment and fill the global_idx with the correct
        //  expression.

        //  Ensure BCs by setting every entry of the rhs corresponding to a
        //  boundary vertex to 0.
        for (const lf::mesh::Entity* vertex : mesh.Entities(2)) {
            if(bd_predicate(*vertex)) {
                int global_idx = dofh.GlobalDofIndices(*vertex)[0];
                rhs[global_idx] = 0.;
            }
        }

        // TODO 3.2 Iterate. After solving this, you should be able to build and
        //  run the code.
        mu = solver.solve(rhs);
        

        solution.row(i) = mu;
    }

    std::ofstream solution_file(CURRENT_SOURCE_DIR "/solution.txt");
    solution_file << solution << std::endl;

    return 0;
}
