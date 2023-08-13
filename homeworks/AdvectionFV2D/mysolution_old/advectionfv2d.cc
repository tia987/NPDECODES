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
std::shared_ptr<lf::mesh::utils::CodimMeshDataSet
               <Eigen::Matrix<double, 2, Eigen::Dynamic>>>
               computeCellNormals(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
    
    //====================
    lf::mesh::utils::CodimMeshDataSet<Eigen::Matrix<double, 2, Eigen::Dynamic>> result(mesh_p,0);
    for(auto *entity : mesh_p->Entities(0)){
        const auto geo = entity->Geometry();
        auto corner = lf::geometry::Corners(*geo);
        if(corner.cols() == 3){
            Eigen::Matrix<double, 2, Eigen::Dynamic> normals(2,3);
            auto bary_coord = gradbarycoordinates(corner);
            normals.col(0) = -(bary_coord.col(2)).normalized();
            normals.col(1) = -(bary_coord.col(0)).normalized();
            normals.col(2) = -(bary_coord.col(1)).normalized();
            result(*entity) = normals;
        } else {
            Eigen::Matrix<double, 2, Eigen::Dynamic> normals(2,4);
            Eigen::Matrix<double,2,3> tri_1;
            Eigen::Matrix<double,2,3> tri_2;
            tri_1 << corner.col(0),corner.col(1),corner.col(2);
            tri_2 << corner.col(0),corner.col(2),corner.col(3);
            auto bary_coord_1 = gradbarycoordinates(tri_1);
            auto bary_coord_2 = gradbarycoordinates(tri_2);
            normals.col(0) = -(bary_coord_1.col(2)).normalized();
            normals.col(1) = -(bary_coord_1.col(0)).normalized();
            normals.col(2) = -(bary_coord_2.col(0)).normalized();
            normals.col(3) = -(bary_coord_2.col(1)).normalized();
            result(*entity) = normals;            
        }        
    }
    return std::make_shared<lf::mesh::utils::CodimMeshDataSet<
  Eigen::Matrix<double, 2, Eigen::Dynamic>>>(result);
    //====================
    
    return nullptr;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::shared_ptr<lf::mesh::utils::CodimMeshDataSet
               <std::array<const lf::mesh::Entity *, 4>>>
               getAdjacentCellPointers(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {

  //====================
  lf::mesh::utils::CodimMeshDataSet<std::array<const lf ::mesh:: Entity *, 2>> aux(mesh_p, 1, {nullptr, nullptr});
  for(auto *over_entity : mesh_p->Entities(0)){
      auto temp = over_entity->SubEntities(1);
      for(auto *entity : temp){
          if(aux(*entity)[0] == nullptr){
              aux(*entity)[0] = over_entity;
          } else if(aux(*entity)[1] == nullptr){
              aux(*entity)[1] = over_entity;
          }
      }
  }

  lf::mesh::utils::CodimMeshDataSet<std::array<const lf ::mesh:: Entity *, 4>> result(mesh_p, 0, {nullptr, nullptr, nullptr, nullptr});
  for(auto *over_entity : mesh_p->Entities(0)){
      auto temp = over_entity->SubEntities(1);
      unsigned counter = 0;
      for(auto *entity : temp){
          if(aux(*entity)[0] != over_entity){
              result(*over_entity)[counter] = aux(*entity)[0];
          } else {
              result(*over_entity)[counter] = aux(*entity)[1];
          }
          counter++;
      }
  }
  return std::make_shared<lf::mesh::utils::CodimMeshDataSet<std::array<const lf::mesh::Entity *, 4>>>(result);
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
    std::vector<double> hmin;
    // std::array<const lf::mesh::Entity *, 4>
    auto AdjacentCellPointers = getAdjacentCellPointers(mesh_p);
    for(auto *entity : mesh_p->Entities(0)){
        const auto geo = entity->Geometry();
        auto corner = lf::geometry::Corners(*geo);
        auto curr_midpoint = barycenter(corner);
        for(auto *cell : (*AdjacentCellPointers)(*entity)){
            if(cell != nullptr){
                auto *geo_neigh = cell->Geometry();
                auto corners_neigh = lf::geometry::Corners(*geo_neigh);
                auto midpoint_neigh = barycenter(corners_neigh);
                double dist = (curr_midpoint-midpoint_neigh).norm();
                hmin.push_back(dist);
            }
        }
    }
    return *std::min_element(std::begin(hmin),std::end(hmin));
    //====================

}
/* SAM_LISTING_END_5 */

}  // namespace AdvectionFV2D
