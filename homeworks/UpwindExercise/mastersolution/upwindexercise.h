// #ifndef UPWIND_QUADRATURE_H
// #define UPWIND_QUADRATURE_H


// #include <lf/base/base.h>
// #include <lf/geometry/geometry.h>
// #include <lf/mesh/mesh.h>
// #include <lf/mesh/utils/utils.h>

// #include <Eigen/Core>
// #include <Eigen/LU>
// #include <memory>
// #include <vector>

// namespace UpwindQuadrature {

// lf::mesh::utils::CodimMeshDataSet<double> initializeMasses(
//     std::shared_ptr<const lf::mesh::Mesh> mesh_p);


// template <typename FUNCTOR>
// class UpwindConvectionElementMatrixProvider {
// public:
//     explicit UpwindConvectionElementMatrixProvider(
//         FUNCTOR v, lf::mesh::utils::CodimMeshDataSet<double> masses)
//         : v_(v), masses_(masses) {}

//     Eigen::Matrix3d Eval(const lf::mesh::Entity &entity);

//     bool isActive(const lf::mesh::Entity & /*entity*/) const { return true; }

// private:
//     FUNCTOR v_;  // velocity field
//     lf::mesh::utils::CodimMeshDataSet<double>
//         masses_;  // masses of all vertices of the mesh.
// };

// /* SAM_LISTING_BEGIN_1 */
// template <typename FUNCTOR>
// Eigen::Matrix3d UpwindConvectionElementMatrixProvider<FUNCTOR>::Eval(
//     const lf::mesh::Entity &entity) {
//     LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
//                   "Function only defined for triangular cells");

//     const lf::geometry::Geometry *geo_ptr = entity.Geometry();
//     const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
//     Eigen::Matrix3d loc_mat;

//     Eigen::Matrix3d grad_helper;
//     grad_helper.col(0) = Eigen::Vector3d::Ones();
//     grad_helper.rightCols(2) = corners.transpose();
//     // Matrix with gradients of the local shape functions in its columns
//     const Eigen::MatrixXd grad_basis = grad_helper.inverse().bottomRows(2);

//     // masses of the corners.
//     std::vector<double> local_masses;
//     // TODO: 2.2.2 collect the local masses
//     for (...) {
//         local_masses.push_back(masses_(...));
//     }

//     // Compute velocities at corners
//     Eigen::MatrixXd velocities(2, 3);
//     // TODO: 2.2.3 calculate the velocity field at the corners
//     velocities << v_(...), v_(...), v_(...);

//     // Matrix of outer normals(not normalized)
//     // TODO: 2.2.4 why does this work?
//     Eigen::MatrixXd n = -grad_basis;

//     // TODO: 2.2.5
//     // compute rows of the local matrix according to the upwind quadrature scheme.
//     // $-v(a^j)$ points into the triangle K
//     // iff the inner product of $v(a^j)$ with the two adjecant outer normals is
//     // positive.
//     for (int i = 0; i < 3; ++i) {
//         Eigen::Vector2d v = velocities.col(i);
//         int left = (i + 2) % 3, right = (i + 1) % 3;
//         if (...) {
//             // if parallel to one edge
//             if (... == 0) {
//                 loc_mat.row(i) = 0.5 * ...;
//             } else {
//                 loc_mat.row(i) = ...;
//             }
//         } else {
//             loc_mat.row(i) = ...;
//         }
//     }

//     return loc_mat;
// }
// /* SAM_LISTING_END_1 */

// }  // namespace UpwindQuadrature

// #endif  // UPWIND_QUADRATURE_H
