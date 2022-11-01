PROJECT(VisualAnalyticTkt)

cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

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

# Use the include path and library for Qt that is used by VTK.
INCLUDE_DIRECTORIES(
        ${QT_INCLUDE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
	${VTK_INCLUDE_DIRS}
)
SET(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")
SET(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
FIND_PACKAGE(Qt4)
IF(QT_USE_FILE)
  INCLUDE(${QT_USE_FILE})
ELSE(QT_USE_FILE)
  SET(QT_LIBRARIES   ${QT_QT_LIBRARY})
ENDIF(QT_USE_FILE)

INCLUDE_DIRECTORIES(
        ${QT_INCLUDE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
	${VTK_INCLUDE_DIRS}
)
 
QT4_WRAP_UI(UI_SRCS GUIV2.ui unit.ui roiAdj.ui defineROI.ui valueChanger.ui) 
QT4_WRAP_CPP(MOC_SRCS GUIV2.h unit.h valueChanger.h roiAdj.h roiDefDial.h SlidingStackedWidget.h)
SET(SRCS  ${MOC_SRCS}  GUIV2.cpp unit.cpp valueChanger.cxx main.cpp roiAdj.cxx roiDefDial.cxx  SlidingStackedWidget.cpp)
 

ADD_EXECUTABLE( VisualAnalyticTkt MACOSX_BUNDLE ${SRCS} ${UI_SRCS})
TARGET_LINK_LIBRARIES( VisualAnalyticTkt
  QVTK
  ${QT_LIBRARIES}
  vtkCharts
  jointModel
  traceMap
  gcm
)

 


