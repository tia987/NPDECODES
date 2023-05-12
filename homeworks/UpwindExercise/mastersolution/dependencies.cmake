#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# C++ sources for this problem 
# DIR will be provided by the calling file.
set(SOURCES
  ${DIR}/upwindexercise_main.cc
  ${DIR}/upwindexercise.h
  ${DIR}/upwindexercise.cc
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
