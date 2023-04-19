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
- Also for the AssembleMatrixLocally, when do we write only once the argument dofh?

## (3-13.j)
```C++
Eigen::Vector2d locGrad = coordinates*(Eigen::Vector3d() << sol_vec[globDof[0]],sol_vec[globDof[1]],sol_vec[globDof[2]]).finished();
``` 
- What does the .finished() actually do?