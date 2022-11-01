#
# This file sets up include directories, link directories, and
# compiler settings for a project to use KML.  It should not be
# included directly, but rather through the KML_USE_FILE setting
# obtained from KMLConfig.cmake.
#

  # Load the compiler settings used for KML.

IF(KML_BUILD_SETTINGS_FILE)
  INCLUDE(${CMAKE_ROOT}/Modules/CMakeImportBuildSettings.cmake)
  CMAKE_IMPORT_BUILD_SETTINGS(${KML_BUILD_SETTINGS_FILE})
ENDIF(KML_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use KML.

  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${KML_REQUIRED_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${KML_REQUIRED_CXX_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${KML_REQUIRED_EXE_LINKER_FLAGS}")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${KML_REQUIRED_SHARED_LINKER_FLAGS}")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${KML_REQUIRED_MODULE_LINKER_FLAGS}")

# Add include directories needed to use KML.
INCLUDE_DIRECTORIES( ${kml_INCLUDE_DIRS})

# Add link directories needed to use KML.
LINK_DIRECTORIES(${kml_LIBRARY_DIRS})


