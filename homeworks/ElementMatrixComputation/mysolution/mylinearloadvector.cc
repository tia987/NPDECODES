/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "mylinearloadvector.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <functional>

namespace ElementMatrixComputation {

namespace {

/* SAM_LISTING_BEGIN_1 */
Eigen::Vector4d computeLoadVector(
    const Eigen::MatrixXd &vertices,
    std::function<double(const Eigen::Vector2d &)> f) {
  // Number of nodes of the element: triangles = 3, rectangles = 4
  const int num_nodes = vertices.cols();
  // Vector for returning element vector
  Eigen::Vector4d elem_vec = Eigen::Vector4d::Zero();

  //====================
  // First, we need to calculate the
  // coordinates for the midpoints in vector form
  // for the quadrature, based on the shape
  Eigen::MatrixXd points(2,num_nodes);
  // Compute area K based on the shape of the mesh
  double K = 0;
  switch(num_nodes){
    case 3:{
      K = ((vertices(0,1)-vertices(0,0))*(vertices(1,2)-vertices(1,0))-(vertices(1,1)-vertices(1,0))*(vertices(0,2)-vertices(0,0)))/2.; 
      points << vertices(0,0)+vertices(0,1), vertices(0,1)+vertices(0,2),
                vertices(0,2)+vertices(0,0), vertices(1,0)+vertices(1,1),
                vertices(1,1)+vertices(1,2), vertices(1,2)+vertices(1,0);
      break;  
    }
    case 4:{
      K = (vertices(0,1)-vertices(0,0))*(vertices(1,3)-vertices(1,0)); 
      points << vertices(0,0)+vertices(0,1), vertices(0,1)+vertices(0,2),
                vertices(0,2)+vertices(0,3), vertices(0,3)+vertices(0,0),
                vertices(1,0)+vertices(1,1), vertices(1,1)+vertices(1,2),
                vertices(1,2)+vertices(1,3), vertices(1,3)+vertices(1,0);
      break;
    }
  }
  points *= 0.5;

  // According to the midpoint rule, we need to set
  // a suitable vector according to the function f
  Eigen::Vector4d vals = Eigen::Vector4d::Zero();
  for(size_t i = 0; i < num_nodes; i++){
    vals(i) = f(points.col(i));
  }
  // Compute quadrature
  Eigen::Vector4d quadrature = Eigen::Vector4d::Zero();
  for(size_t i = 0; i < num_nodes; i++){
    quadrature(i) += vals(i);

    // Why we need this step?
    quadrature((i+1) % num_nodes) += vals(i);
  }
  quadrature *= K / num_nodes;
  elem_vec = quadrature*0.5;
  //====================

  return elem_vec;
}
/* SAM_LISTING_END_1 */

}  // namespace

Eigen::Vector4d MyLinearLoadVector::Eval(const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  const lf::base::size_type num_nodes{ref_el.NumNodes()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();

  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());

  return computeLoadVector(vertices, f_);
}

}  // namespace ElementMatrixComputation
