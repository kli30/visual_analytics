#-----------------------------------------------------------------------------
#
# GSLConfig.cmake - GSL CMake configuration file for external projects.
#
# This file is configured by GSL and used by the UseGSL.cmake module
# to load GSL's settings for an external project.

# The GSL include file directories.
SET(GSL_INCLUDE_DIRS "${GSL_DIR}/include;")

# The GSL library directories.
SET(GSL_LIBRARY_DIRS "${GSL_DIR}/lib64;${GSL_DIR}/lib;")

# The C and C++ flags added by GSL to the cmake-configured flags.
SET(GSL_REQUIRED_C_FLAGS "-DNOMINMAX -DHAVE_ZLIB")
SET(GSL_REQUIRED_CXX_FLAGS "-DNOMINMAX -DHAVE_ZLIB ")
SET(GSL_REQUIRED_LINK_FLAGS "")
  
# The location of the UseGSL.cmake file.
SET(GSL_USE_FILE "${GSL_DIR}/UseGSL.cmake")
SET(USE_GSL_FILE "${GSL_DIR}/UseGSL.cmake")

# Whether GSL was built with shared libraries.
SET(GSL_BUILD_SHARED "OFF")


# The GSL library dependencies.
IF(NOT GSL_NO_LIBRARY_DEPENDS AND
    EXISTS "${GSL_DIR}/GSLLibraryDepends.cmake")
  INCLUDE("${GSL_DIR}/GSLLibraryDepends.cmake")
ENDIF(NOT GSL_NO_LIBRARY_DEPENDS AND
  EXISTS "${GSL_DIR}/GSLLibraryDepends.cmake")

 

