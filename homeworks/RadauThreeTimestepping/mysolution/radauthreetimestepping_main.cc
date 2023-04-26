/**
 * @file radauthreetimestepping_main.cc
 * @brief NPDE homework RadauThreeTimestepping
 * @author Erick Schulz
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <iostream>
#include <memory>

#include "radauthreetimestepping.h"
#include "radauthreetimesteppingode.h"

using namespace RadauThreeTimestepping;

int main(int /*argc*/, char ** /*argv*/) {
  /* Solving the ODE problem */
  // This function prints to the terminal the convergence rates and average rate
  // of a convergence study performed for the ODE (d/dt)y = -y.
  testConvergenceTwoStageRadauLinScalODE();

  /* Solving the parabolic heat equation */
  // Create a Lehrfem++ square tensor product mesh
  lf::mesh::utils::TPTriagMeshBuilder builder(
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
  // Set mesh parameters following the Builder pattern
  // Domain is the unit square
  builder.setBottomLeftCorner(Eigen::Vector2d{-1.0, -1.0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNumXCells(50)
      .setNumYCells(50);
  auto mesh_p = builder.Build();

  /* SAM_LISTING_BEGIN_1 */
  //====================
  //Eigen::VectorXd solveHeatEvolution(const lf::assemble::DofHandler &dofh,
//                                   unsigned int m, double final_time);
  auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
// Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler &dofh {fe_space->LocGlobMap()};
  Eigen::VectorXd output = solveHeatEvolution(dofh, 1., 10);
  lf::io::VtkWriter vtk_writer (mesh_p , CURRENT_BINARY_DIR "/discrete_heat_solution.vtk ");
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < dofh.NumDofs(); global_idx++) {
    nodal_data->operator()(dofh.Entity(global_idx)) =
        output[global_idx];
  };
  vtk_writer.WritePointData("discrete_heat_solution", *nodal_data);
  std::cout << "\n The discrete_heat_solution was written to:" << std::endl;
  std::cout << ">> discrete_heat_solution.vtk\n" << std::endl;
  //====================

  return 0;
}
