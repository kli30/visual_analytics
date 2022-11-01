PROJECT(fsl)
cmake_minimum_required(VERSION 2.8)


IF(NOT EXISTS   ${FSL_SOURCE_DIR})
    message("FSL_SOURCE_DIR: ${FSL_SOURCE_DIR}")
	set(FSL_SOURCE_DIR "" CACHE PATH "Please choose the root directory of fsl package, which is under root direcotry of visual analytics" )
	MESSAGE(FATAL_ERROR "Please choose the root directory of fsl package, which is under root direcotry of visual analytics")
ENDIF(NOT EXISTS  ${FSL_SOURCE_DIR})

INCLUDE_DIRECTORIES(${FSL_SOURCE_DIR}/src ${FSL_SOURCE_DIR}/extras/src ${FSL_SOURCE_DIR}/extras/src/zlib ${FSL_SOURCE_DIR}/src/znzlib )

SUBDIRS(extras src)
