SET(MISCMATHS_SOURCES
  miscmaths.cc
  optimise.cc
  miscprob.cc
  kernel.cc
  histogram.cc
  base2z.cc
  t2z.cc
  f2z.cc
  minimize.cc
  cspline.cc
  sparse_matrix.cc
  sparsefn.cc
  rungekutta.cc)

INCLUDE_DIRECTORIES( ${FSL_SOURCE_DIR}/extras/src/libprob ${FSL_SOURCE_DIR}/extras/src/newmat)

ADD_LIBRARY(miscmaths STATIC ${MISCMATHS_SOURCES})

install(DIRECTORY ./ DESTINATION ${FSL_SOURCE_DIR}/include/  FILES_MATCHING PATTERN "*.h")
install(TARGETS miscmaths DESTINATION ${FSL_SOURCE_DIR}/lib)
