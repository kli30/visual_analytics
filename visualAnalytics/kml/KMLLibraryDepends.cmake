SET(triSurface_LIB_DEPENDS "vtkRendering;vtkGraphics;vtkImaging;vtkIO;vtkFiltering;vtkCommon;")
SET(SurfacePatchSampler_LIB_DEPENDS "triSurface;newmat;")

SET(fibers_LIB_DEPENDS "kmc;indexer;")
SET(kmc_LIB_DEPENDS "newmat;vtkRendering;vtkGraphics;vtkImaging;vtkIO;vtkFiltering;vtkCommon;gsl")
SET(kmPCA_LIB_DEPENDS "alglib;")
SET(gsl_LIB_DEPENDS "gslcblas;")
SET(traceMap_LIB_DEPENDS "kmPCA;")
SET(jointModel_LIB_DEPENDS "triSurface;fibers;newimage;colorScheme;")
SET(colorScheme_LIB_DEPENDS "newran;")

 
