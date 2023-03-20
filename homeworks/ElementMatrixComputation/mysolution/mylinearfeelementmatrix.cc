/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "mylinearfeelementmatrix.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>

namespace ElementMatrixComputation {

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 4, 4> MyLinearFEElementMatrix::Eval(
    const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());
  // Matrix for returning element matrix
  Eigen::Matrix<double, 4, 4> elem_mat;

  //====================
  // First we need to specify which case we have to compute
  // through the kind of ref_el      
  switch(ref_el){
    //case ref_el.kTria():{
    case lf::base::RefEl::kTria():{
      // Triangle
      // Calculate area for triangle
      double K = ((vertices(0,1)-vertices(0,0))*(vertices(1,2)-vertices(1,0))-(vertices(1,1)-vertices(1,0))*(vertices(0,2)-vertices(0,0)))/2.; 
      // Setup matrix (Need to be 4x4 given the laplace_elmat_builder)
      elem_mat << 2,1,1,0,
                  1,2,1,0,
                  1,1,2,0,
                  0,0,0,0;
      elem_mat = elem_mat*K/12.;   
      break; 
    }

    case lf::base::RefEl::kQuad():{
      // square
      // Calculate area for square
      double K = (vertices(0,1)-vertices(0,0))*(vertices(1,3)-vertices(1,0)); 
      elem_mat << 4,2,1,2,
                  2,4,2,1,
                  1,2,4,2,
                  2,1,2,4;
      elem_mat = elem_mat*K/36.;
      break;    
    }
  }
  lf::uscalfe::LinearFELaplaceElementMatrix laplace_elmat_builder;
  auto lap_elem_mat = laplace_elmat_builder.Eval(cell);
  elem_mat += lap_elem_mat;
  //====================

  return elem_mat;
}
/* SAM_LISTING_END_1 */
}  // namespace ElementMatrixComputation
