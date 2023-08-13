/**
 * @file advectionfv2d.cc
 * @brief NPDE homework AdvectionFV2D code
 * @author Philipp Egg
 * @date 21.06.2020
 * @copyright Developed at ETH Zurich
 */

#include "advectionfv2d.h"

#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <memory>
#include <stdexcept>
#include <vector>

namespace AdvectionFV2D {

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 2, 3> gradbarycoordinates(
    const Eigen::Matrix<double, 2, 3> &triangle) {
  Eigen::Matrix3d X;

  // solve for the coefficients of the barycentric coordinate functions
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = triangle.transpose();
  return X.inverse().block<2, 3>(1, 0);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<Eigen::Matrix<double, 2, Eigen::Dynamic>>>
    computeCellNormals(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
    //====================  
    lf::mesh::utils::CodimMeshDataSet<Eigen::Matrix<double, 2, Eigen::Dynamic>> norm_store(mesh_p,0);
    
    for(const lf::mesh::Entity *cell : mesh_p->Entities(0)){
        const lf::geometry::Geometry *geo_p = cell->Geometry();
        const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);

        if(corners.cols() == 3){
            Eigen::Matrix<double,2,Eigen::Dynamic> normal_vectors(2,3);
            Eigen::Matrix<double,2,3> bary_coord = gradbarycoordinates(corners);
            normal_vectors.col(0) = -bary_coord.col(2).normalized();
            normal_vectors.col(1) = -bary_coord.col(0).normalized();
            normal_vectors.col(2) = -bary_coord.col(1).normalized();
            norm_store(*cell) = normal_vectors;
        } else {
            Eigen::Matrix<double,2,Eigen::Dynamic> normal_vectors(2,4);
            Eigen::Matrix<double,2,3> tri_1;
            Eigen::Matrix<double,2,3> tri_2;
            tri_1 << corners.col(0), corners.col(1), corners.col(2);
            tri_2 << corners.col(0), corners.col(2), corners.col(3);
            Eigen::Matrix<double,2,3> bary_coord_1 = gradbarycoordinates(tri_1);
            Eigen::Matrix<double,2,3> bary_coord_2 = gradbarycoordinates(tri_2);
            normal_vectors.col(0) = -bary_coord_1.col(2).normalized();
            normal_vectors.col(1) = -bary_coord_1.col(0).normalized();
            normal_vectors.col(2) = -bary_coord_2.col(0).normalized();
            normal_vectors.col(3) = -bary_coord_2.col(1).normalized();
            norm_store(*cell) = normal_vectors;
        }
    }
    return std::make_shared<lf::mesh::utils::CodimMeshDataSet<Eigen::Matrix<double,2,Eigen::Dynamic>>>(norm_store);
    //====================
    return nullptr;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity *, 4>>>
    getAdjacentCellPointers(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
    //====================
    lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity *, 2>> prt_reducted(mesh_p,1,{nullptr, nullptr});
    
    for(const lf::mesh::Entity *cell : mesh_p->Entities(0)){
        auto cell_edges = cell->SubEntities(1);
        for(const lf::mesh::Entity *edge : cell_edges){
            if(prt_reducted(*edge)[0] == nullptr) {
                prt_reducted(*edge)[0] = cell;
            } else if(prt_reducted(*edge)[1] == nullptr) {
                prt_reducted(*edge)[1] = cell;
            } else {
               throw std::runtime_error("FUCK");
            }
        }
    }
    lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity *, 4>> ptr_store(mesh_p,0,{nullptr, nullptr, nullptr, nullptr});
    for(const lf::mesh::Entity *cell : mesh_p->Entities(0)){
        auto cell_edges = cell->SubEntities(1);
        unsigned counter = 0;
        for(const lf::mesh::Entity *edge : cell_edges){
            if(prt_reducted(*edge)[0] != cell){
                  ptr_store(*cell)[counter] = prt_reducted(*edge)[0];
            } else {
                  ptr_store(*cell)[counter] = prt_reducted(*edge)[1];
            }
            counter++;
        }
    }
    return std::make_shared<lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity *, 4>>>(ptr_store);
    //====================
    return nullptr;
}
/* SAM_LISTING_END_3 */

// Function returning the barycenter of TRIA or QUAD
/* SAM_LISTING_BEGIN_4 */
Eigen::Vector2d barycenter(const Eigen::MatrixXd corners) {
  Eigen::Vector2d midpoint;
  if (corners.cols() == 3) {
    midpoint = (corners.col(0) + corners.col(1) + corners.col(2)) / 3.0;
  } else if (corners.cols() == 4) {
    midpoint =
        (corners.col(0) + corners.col(1) + corners.col(2) + corners.col(3)) /
        4.0;
  } else {
    throw std::runtime_error("Wrong geometrie in barycenter()");
  }
  return midpoint;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
double computeHmin(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  //====================
  // Your code goes here
  // Replace the dummy return value below:
  return 0.0;
  //====================
}
/* SAM_LISTING_END_5 */

}  // namespace AdvectionFV2D
