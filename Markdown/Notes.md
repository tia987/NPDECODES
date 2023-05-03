# Questions
## (3-13.e) 
```C++
 lf::uscalfe::ReactionDiffusionElementMatrixProvider<double,
                                                      decltype(mf_sigma), 
                                                      decltype(mf_gamma)> 
                                                      basis_expansion(fe_space,mf_sigma,mf_gamma);
  lf::assemble::AssembleMatrixLocally(0,dofh,dofh,basis_expansion,A);
  double contactPoints = voltvals.size();
  auto selector = [&](auto id)->std::pair<bool,double>{
    const lf::mesh::Entity &node{dofh.Entity(id)}; // Get the entity of the node
    int lim = nodeflags(node); // return 1,0 or -1
    if ((lim >= 0) && (lim < contactPoints)) return {true, voltvals[lim]}; // If nodeflags values are valid return true
    else return {false, -1.}; 
  };
  phi.setZero();
  lf::assemble::FixFlaggedSolutionComponents<double>(selector,A,phi);
```
- why do we need to set the selector in such way?
- Also for the AssembleMatrixLocally, when do we write only once the argument dofh? and what does AssembleMatrixLocally / AssembleVectorLocally do?

## (3-13.j)
```C++
Eigen::Vector2d locGrad = coordinates*(Eigen::Vector3d() << sol_vec[globDof[0]],sol_vec[globDof[1]],sol_vec[globDof[2]]).finished();
``` 
- What does the .finished() actually do?

## (5-3.g)
```C++
for (int i = 0; i < nvals; ++i) {
const Eigen : : Vector2d g{ gradvals [ i ] } ; const double norms_g = g . squaredNorm ( ) ; Avals[i] =
                    }
return
1.0 / std::sqrt(1.0 + norms_g) *
(Eigen::Matrix2d::Identity() − g * g.transpose() / (1.0 + norms_g));
```
- How did they come up with this Solution??

## (9-1.f)
```C++
// Deduction guide for TrapRuleLinFEElemVecProvider
template <typename FUNCTOR>
TrapRuleLinFEElemVecProvider(FUNCTOR) -> TrapRuleLinFEElemVecProvider<FUNCTOR>;

// TrapRuleLinFEElemVecProvider
/* Implementing member function Eval of class TrapRuleLinFEElemVecProvider*/
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR>
Eigen::Vector3d TrapRuleLinFEElemVecProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity &tria) {
  Eigen::Vector3d ElemVec;
  //====================
  // Throw error in case no triangular cell
  LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
		  "Unsupported cell type " << tria.RefEl());
  // Obtain vertex coordinates of the triangle in a 2x3 matrix
  const auto corners{lf::geometry::Corners(*(tria.Geometry()))};
  const double area_third = lf::geometry::Volume(*(tria.Geometry())) / 3.0;
  LF_ASSERT_MSG((corners.cols() == 3) && (corners.rows() == 2),
		  "Invalid vertex coordinate " << corners.rows() << "x"
		  << corners.cols() << " matrix");
  ElemVec << area_third*f_(corners.col(0)),
             area_third*f_(corners.col(1)),
             area_third*f_(corners.col(2));
  //====================
  return ElemVec;
}
```

- Shouldn't we use the trapezoidal rule in this case? Why it doesn't hange from the normal LinFEElemVecProvider?

```C++
auto mesh_p = dofh.Mesh();
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  for(auto vertex : mesh_p->Entities(1)){
    if(bd_flags(*vertex)){
      auto index = dofh.GlobalDofIndices(*vertex);
      phi(index[0]) = 0.;
    }
  }
```

- `mesh_p->Entities(1)` why this works even if we set 1? What should it represent?
- Also for bd_flags, how does this work?
- I also still struggle to understand a little bit the function `GlobalDofIndices`, so basically this gives us the index related to the RHS vector. But what is exactly its type? And also why at 0 we should find the exact position for phi?

```C++
/* Implementing member function Eval of class LinFEMassMatrixProvider*/
Eigen::Matrix<double, 3, 3> LinFEMassMatrixProvider::Eval(
    const lf::mesh::Entity &tria) {
  Eigen::Matrix<double, 3, 3> elMat;
  //====================
  // Throw error in case no triangular cell
  LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
		  "Unsupported cell type " << tria.RefEl());
  // Obtain vertex coordinates of the triangle in a 2x3 matrix
  const double area = lf::geometry::Volume(*(tria.Geometry()));
  LF_ASSERT_MSG((corners.cols() == 3) && (corners.rows() == 2),
		  "Invalid vertex coordinate " << corners.rows() << "x"
		  << corners.cols() << " matrix");
   elMat << 2.0, 1.0, 1.0,
           1.0, 2.0, 1.0,
           1.0, 1.0, 2.0;
  // clang-format on
  elMat *= area / 12.0;;
  //====================
  return elMat;  // return the local mass element matrix
```

- Why do we set the `elMat` in this way? In the exercise questions there is no indication about it.

