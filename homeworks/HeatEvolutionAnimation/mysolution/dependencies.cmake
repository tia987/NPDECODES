#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# C++ sources for this problem 
# DIR will be provided by the calling file.
set(SOURCES
  ${DIR}/heatevolutionanimation_main.cc
  ${DIR}/heatevolutionanimation.h
  ${DIR}/heatevolutionanimation.cc
)

# Libraries to be used. If the code does not rely on LehrFEM++
# all the libraries LF:* can be removed 
set(LIBRARIES
  Eigen3::Eigen
  LF::lf.assemble
  LF::lf.base
  LF::lf.geometry
  LF::lf.io
  LF::lf.quad
  LF::lf.mesh.hybrid2d
  LF::lf.mesh.test_utils
  LF::lf.mesh.utils
  LF::lf.refinement
  LF::lf.fe
  LF::lf.uscalfe
)
