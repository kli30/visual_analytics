# Install script for directory: /home/kaiming/programs/visualAnalytics/kml/src

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/kaiming/programs/visualAnalytics/kml/lib/libalglib.a")
FILE(INSTALL DESTINATION "/home/kaiming/programs/visualAnalytics/kml/lib" TYPE STATIC_LIBRARY FILES "/home/kaiming/programs/visualAnalytics/kml/src/build/libalglib.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/kaiming/programs/visualAnalytics/kml/lib/libkmc.a")
FILE(INSTALL DESTINATION "/home/kaiming/programs/visualAnalytics/kml/lib" TYPE STATIC_LIBRARY FILES "/home/kaiming/programs/visualAnalytics/kml/src/build/libkmc.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/kaiming/programs/visualAnalytics/kml/lib/libgcm.a")
FILE(INSTALL DESTINATION "/home/kaiming/programs/visualAnalytics/kml/lib" TYPE STATIC_LIBRARY FILES "/home/kaiming/programs/visualAnalytics/kml/src/build/libgcm.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/kaiming/programs/visualAnalytics/kml/lib/libxtPCACCA.a")
FILE(INSTALL DESTINATION "/home/kaiming/programs/visualAnalytics/kml/lib" TYPE STATIC_LIBRARY FILES "/home/kaiming/programs/visualAnalytics/kml/src/build/libxtPCACCA.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/kaiming/programs/visualAnalytics/kml/lib/libkmPCA.a")
FILE(INSTALL DESTINATION "/home/kaiming/programs/visualAnalytics/kml/lib" TYPE STATIC_LIBRARY FILES "/home/kaiming/programs/visualAnalytics/kml/src/build/libkmPCA.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/kaiming/programs/visualAnalytics/kml/lib/libindexer.a")
FILE(INSTALL DESTINATION "/home/kaiming/programs/visualAnalytics/kml/lib" TYPE STATIC_LIBRARY FILES "/home/kaiming/programs/visualAnalytics/kml/src/build/libindexer.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/kaiming/programs/visualAnalytics/kml/lib/libtraceMap.a")
FILE(INSTALL DESTINATION "/home/kaiming/programs/visualAnalytics/kml/lib" TYPE STATIC_LIBRARY FILES "/home/kaiming/programs/visualAnalytics/kml/src/build/libtraceMap.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/kaiming/programs/visualAnalytics/kml/lib/libfibers.a")
FILE(INSTALL DESTINATION "/home/kaiming/programs/visualAnalytics/kml/lib" TYPE STATIC_LIBRARY FILES "/home/kaiming/programs/visualAnalytics/kml/src/build/libfibers.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/kaiming/programs/visualAnalytics/kml/lib/libtriSurface.a")
FILE(INSTALL DESTINATION "/home/kaiming/programs/visualAnalytics/kml/lib" TYPE STATIC_LIBRARY FILES "/home/kaiming/programs/visualAnalytics/kml/src/build/libtriSurface.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/kaiming/programs/visualAnalytics/kml/lib/libjointModel.a")
FILE(INSTALL DESTINATION "/home/kaiming/programs/visualAnalytics/kml/lib" TYPE STATIC_LIBRARY FILES "/home/kaiming/programs/visualAnalytics/kml/src/build/libjointModel.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/kaiming/programs/visualAnalytics/kml/lib/libcolorScheme.a")
FILE(INSTALL DESTINATION "/home/kaiming/programs/visualAnalytics/kml/lib" TYPE STATIC_LIBRARY FILES "/home/kaiming/programs/visualAnalytics/kml/src/build/libcolorScheme.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/kaiming/programs/visualAnalytics/kml/include//")
FILE(INSTALL DESTINATION "/home/kaiming/programs/visualAnalytics/kml/include/" TYPE DIRECTORY FILES "/home/kaiming/programs/visualAnalytics/kml/src/./" FILES_MATCHING REGEX "/[^/]*\\.h$")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/home/kaiming/programs/visualAnalytics/kml/src/build/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/home/kaiming/programs/visualAnalytics/kml/src/build/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
