#-----------------------------------------------------------------------------
#
# FSLConfig.cmake - FSL CMake configuration file for external projects.
#
# This file is configured by FSL and used by the UseFSL.cmake module
# to load FSL's settings for an external project.

# The FSL include file directories.
SET(FSL_INCLUDE_DIRS ${FSL_DIR}/include;${FSL_DIR}/src/;${FSL_DIR}/extras/include/boost;${FSL_DIR}/extras/src/;)

# The FSL library directories.
SET(FSL_LIBRARY_DIRS ${FSL_DIR}/lib)

# The C and C++ flags added by FSL to the cmake-configured flags.
SET(FSL_REQUIRED_C_FLAGS "-DNOMINMAX -DHAVE_ZLIB")
SET(FSL_REQUIRED_CXX_FLAGS "-DNOMINMAX -DHAVE_ZLIB ")
SET(FSL_REQUIRED_LINK_FLAGS "")
 


# The FSL version number
SET(FSL_VERSION_MAJOR "4")
SET(FSL_VERSION_MINOR "0")
SET(FSL_VERSION_PATCH "0")

# The location of the UseFSL.cmake file.
SET(FSL_USE_FILE ${FSL_DIR}/UseFSL.cmake)
SET(USE_FSL_FILE ${FSL_DIR}/UseFSL.cmake)


# Whether FSL was built with shared libraries.
SET(FSL_BUILD_SHARED "OFF")


# The FSL library dependencies.
IF(NOT FSL_NO_LIBRARY_DEPENDS AND
    EXISTS ${FSL_DIR}/FSLLibraryDepends.cmake)
  INCLUDE(${FSL_DIR}/FSLLibraryDepends.cmake)
ENDIF(NOT FSL_NO_LIBRARY_DEPENDS AND
  EXISTS ${FSL_DIR}/FSLLibraryDepends.cmake)

 

