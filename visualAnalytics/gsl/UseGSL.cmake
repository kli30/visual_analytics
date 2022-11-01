#
# This file sets up include directories, link directories, and
# compiler settings for a project to use GSL.  It should not be
# included directly, but rather through the GSL_USE_FILE setting
# obtained from GSLConfig.cmake.
#

  # Load the compiler settings used for GSL.

IF(GSL_BUILD_SETTINGS_FILE)
  INCLUDE(${CMAKE_ROOT}/Modules/CMakeImportBuildSettings.cmake)
  CMAKE_IMPORT_BUILD_SETTINGS(${GSL_BUILD_SETTINGS_FILE})
ENDIF(GSL_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use GSL.

  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GSL_REQUIRED_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GSL_REQUIRED_CXX_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GSL_REQUIRED_EXE_LINKER_FLAGS}")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${GSL_REQUIRED_SHARED_LINKER_FLAGS}")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${GSL_REQUIRED_MODULE_LINKER_FLAGS}")

# Add include directories needed to use GSL.
INCLUDE_DIRECTORIES( ${GSL_INCLUDE_DIRS})

# Add link directories needed to use GSL.
LINK_DIRECTORIES(${GSL_LIBRARY_DIRS})


