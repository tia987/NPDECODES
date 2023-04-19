set(SOURCES
${DIR}/minimalgraphsurface_main.cc
${DIR}/minimalgraphsurface.cc
${DIR}/minimalgraphsurface.h
${DIR}/preptrimesh.h
${DIR}/preptrimesh.cc)
set(LIBRARIES Eigen3::Eigen LF::lf.base LF::lf.mesh  LF::lf.geometry  LF::lf.mesh.hybrid2d  LF::lf.mesh.utils  LF::lf.mesh.test_utils  LF::lf.refinement  LF::lf.assemble  LF::lf.quad  LF::lf.io  LF::lf.fe  LF::lf.uscalfe)
