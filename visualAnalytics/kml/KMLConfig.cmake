#-----------------------------------------------------------------------------
#
# kmlConfig.cmake - kml CMake configuration file for external projects.
#
# This file is configured by kml and used by the UseKML.cmake module
# to load kml's settings for an external project.

# The kml include file directories.
SET(kml_INCLUDE_DIRS "${KML_DIR}/include;" "${KML_DIR}/src;")

# The kml library directories.
SET(kml_LIBRARY_DIRS "${KML_DIR}/lib;")

# The C and C++ flags added by kml to the cmake-configured flags.
SET(kml_REQUIRED_C_FLAGS "-DNOMINMAX -DHAVE_ZLIB")
SET(kml_REQUIRED_CXX_FLAGS "-DNOMINMAX -DHAVE_ZLIB ")
SET(kml_REQUIRED_LINK_FLAGS "")
  
SET(kml_VERSION_PATCH "0")

# The location of the UseKML.cmake file.
SET(KML_USE_FILE "${KML_DIR}/UseKML.cmake")
SET(USE_KML_FILE "${KML_DIR}/UseKML.cmake")

# The build settings file.
SET(kml_BUILD_SETTINGS_FILE "${KML_DIR}/KMLBuildSettings.cmake")

# Whether kml was built with shared libraries.
SET(kml_BUILD_SHARED "OFF")


# The kml library dependencies.
IF(NOT kml_NO_LIBRARY_DEPENDS AND
    EXISTS "${KML_DIR}/KMLLibraryDepends.cmake")
  INCLUDE("${KML_DIR}/KMLLibraryDepends.cmake")
ENDIF(NOT kml_NO_LIBRARY_DEPENDS AND
  EXISTS "${KML_DIR}/KMLLibraryDepends.cmake")

 

