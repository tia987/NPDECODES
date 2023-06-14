/**
 * @file
 * @brief NPDE homework "Handling degrees of freedom (DOFs) in LehrFEM++"
 * @author Julien Gacon
 * @date March 1st, 2019
 * @copyright Developed at ETH Zurich
 */

#include "lfppdofhandling.h"

#include <Eigen/Dense>
#include <array>
#include <memory>

#include "lf/assemble/assemble.h"
#include "lf/base/base.h"
#include "lf/geometry/geometry.h"
#include "lf/mesh/mesh.h"
#include "lf/mesh/utils/utils.h"

namespace LFPPDofHandling {

/* SAM_LISTING_BEGIN_1 */
std::array<std::size_t, 3> countEntityDofs(
    const lf::assemble::DofHandler &dofhandler) {
    std::array<std::size_t, 3> entityDofs;
    //====================
    auto mesh = dofhandler.Mesh();
    for(unsigned codim = 0; codim < 3; codim++){
        entityDofs[codim] = 0;
        for(const auto *entity : mesh->Entities(codim)){
            entityDofs[codim] += dofhandler.NumInteriorDofs(*entity);            
        }
    }
    //====================
    return entityDofs;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::size_t countBoundaryDofs(const lf::assemble::DofHandler &dofhandler) {
    std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
    // given an entity, bd\_flags(entity) == true, if the entity is on the
    // boundary
    lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags(lf::mesh::utils::flagEntitiesOnBoundary(mesh));
    std::size_t no_dofs_on_bd = 0;
    //====================
    for(unsigned codim = 0; codim < 3; codim++){
        for(const auto *entity : mesh->Entities(codim)){
            if(bd_flags(*entity)){
               no_dofs_on_bd += dofhandler.NumInteriorDofs(*entity);
            }
        }
    }
    //====================
    return no_dofs_on_bd;
}
/* SAM_LISTING_END_2 */

// clang-format off
/* SAM_LISTING_BEGIN_3 */
double integrateLinearFEFunction(
    const lf::assemble::DofHandler& dofhandler,
    const Eigen::VectorXd& mu) {
    double I = 0;
    //====================
    std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
    for(unsigned codim = 0; codim < 1; codim++){
        for(const auto *entity : mesh->Entities(codim)){
            double area = lf::geometry::Volume(*(entity->Geometry())); // In this case only the area of the triangle matters
            for(auto idx : dofhandler.GlobalDofIndices(*entity)){
                I += area/3.*mu(idx);
            }
        }
    }
    //====================
    return I;
}
/* SAM_LISTING_END_3 */
// clang-format on

/* SAM_LISTING_BEGIN_4 */
double integrateQuadraticFEFunction(const lf::assemble::DofHandler &dofhandler,
                                    const Eigen::VectorXd &mu) {
    double I = 0;
    //====================
    std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
    for(unsigned codim = 0; codim < 1; codim++){
        for(const auto *entity : mesh->Entities(codim)){
            double area = lf::geometry::Volume(*(entity->Geometry())); // In this case only the area of the triangle matters
            for(unsigned idx = 3; idx < 6; idx++){
                I += area/3.*mu(idx);
            }
        }
    }
    //====================
    return I;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd convertDOFsLinearQuadratic(
    const lf::assemble::DofHandler &dofh_Linear_FE,
    const lf::assemble::DofHandler &dofh_Quadratic_FE,
    const Eigen::VectorXd &mu) {
    if (dofh_Linear_FE.Mesh() != dofh_Quadratic_FE.Mesh()) {
      throw "Underlying meshes must be the same for both DOF handlers!";
    }
    std::shared_ptr<const lf::mesh::Mesh> mesh =
        dofh_Linear_FE.Mesh();                          // get the mesh
    Eigen::VectorXd zeta(dofh_Quadratic_FE.NumDofs());  // initialise empty zeta
    // safety guard: always set zero if you're not sure to set every entry later
    // on for us this shouldn't be a problem, but just to be sure
    zeta.setZero();

    for (const auto *cell : mesh->Entities(0)) {
        // check if the spaces are actually linear and quadratic
        //====================
        
        //====================
        // get the global dof indices of the linear and quadratic FE spaces, note
        // that the vectors obey the LehrFEM++ numbering, which we will make use of
        // lin\_dofs will have size 3 for the 3 dofs on the nodes and
        // quad\_dofs will have size 6, the first 3 entries being the nodes and
        // the last 3 the edges
        //====================
        // Your code goes here
        // assign the coefficients of mu to the correct entries of zeta, use
        // the previous subproblem 2-9.a
        //====================
    }
    return zeta;
}
/* SAM_LISTING_END_5 */

}  // namespace LFPPDofHandling
