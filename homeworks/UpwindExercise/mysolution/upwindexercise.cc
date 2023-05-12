/**
* @file upwindquadrature.cc
* @brief NPDE homework template
* @author Philippe Peter
* @date June 2020
* @copyright Developed at SAM, ETH Zurich
*/

#include "upwindexercise.h"

namespace UpwindQuadrature {


lf::mesh::utils::CodimMeshDataSet<double> initializeMasses(
   std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
   lf::mesh::utils::CodimMeshDataSet<double> masses(mesh_p, 2, 0.0);
   // compute masses using a cell-based approach.
   for (const lf::mesh::Entity *entity : mesh_p->Entities(0)) {
       const lf::geometry::Geometry *geo_ptr = entity->Geometry();
       double area = lf::geometry::Volume(*geo_ptr);
       for (const lf::mesh::Entity *corner : entity->SubEntities(2)) {
           // TODO: 2.2.1 compute the masses accordingly
           masses(*corner) += 0.;
       }
   }
   return masses;
}

}  // namespace UpwindQuadrature
