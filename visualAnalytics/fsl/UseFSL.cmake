#
# This file sets up include directories, link directories, and
# compiler settings for a project to use FSL.  It should not be
# included directly, but rather through the FSL_USE_FILE setting
# obtained from FSLConfig.cmake.
#

  # Load the compiler settings used for FSL.

IF(FSL_BUILD_SETTINGS_FILE)
  INCLUDE(${CMAKE_ROOT}/Modules/CMakeImportBuildSettings.cmake)
  CMAKE_IMPORT_BUILD_SETTINGS(${FSL_BUILD_SETTINGS_FILE})
ENDIF(FSL_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use FSL.

  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${FSL_REQUIRED_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FSL_REQUIRED_CXX_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${FSL_REQUIRED_EXE_LINKER_FLAGS}")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${FSL_REQUIRED_SHARED_LINKER_FLAGS}")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${FSL_REQUIRED_MODULE_LINKER_FLAGS}")

# Add include directories needed to use FSL.
INCLUDE_DIRECTORIES( ${FSL_INCLUDE_DIRS})

# Add link directories needed to use FSL.
LINK_DIRECTORIES(${FSL_LIBRARY_DIRS})


