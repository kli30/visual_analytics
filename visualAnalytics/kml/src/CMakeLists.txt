PROJECT(KML)
cmake_minimum_required(VERSION 2.6)



 
FIND_PACKAGE(VTK REQUIRED)

IF(VTK_FOUND)

	INCLUDE( ${USE_VTK_FILE} )

ENDIF(VTK_FOUND)

FIND_PACKAGE(FSL REQUIRED)

IF(FSL_FOUND)

	INCLUDE( ${USE_FSL_FILE} )

ENDIF(FSL_FOUND)

FIND_PACKAGE(KML REQUIRED)

IF(KML_FOUND)

	INCLUDE( ${USE_KML_FILE} )

ENDIF(KML_FOUND)

FIND_PACKAGE(GSL REQUIRED)

IF(GSL_FOUND)

	INCLUDE( ${USE_GSL_FILE} )

ENDIF(GSL_FOUND)


set (srcs_alglib alglibinternal.cpp
alglibmisc.cpp
ap.cpp
dataanalysis.cpp
diffequations.cpp
fasttransforms.cpp
integration.cpp
interpolation.cpp
linalg.cpp
optimization.cpp
solvers.cpp
specialfunctions.cpp
statistics.cpp)

ADD_LIBRARY(alglib ${srcs_alglib} )
ADD_LIBRARY(kmc kaimingCommon.cpp)
ADD_LIBRARY(gcm grangerCausality.cpp)
ADD_LIBRARY(xtPCACCA pcacca.cpp)
ADD_LIBRARY(kmPCA kmPCA.cpp)
ADD_LIBRARY(traceMap traceMap.cpp)
ADD_LIBRARY(indexer indexer.cpp)
ADD_LIBRARY(fibers fibers.cpp )
ADD_LIBRARY(triSurface triSurface.cpp )
ADD_LIBRARY(jointModel jointModel.cpp )
ADD_LIBRARY(colorScheme colorScheme.cpp)


install(TARGETS alglib kmc gcm xtPCACCA kmPCA indexer traceMap fibers triSurface jointModel colorScheme DESTINATION ${KML_DIR}/lib)
install(DIRECTORY ./ DESTINATION ${KML_DIR}/include/ FILES_MATCHING PATTERN "*.h")


 


