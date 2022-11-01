SET(NEWIMAGE_SOURCES
	complexvolume.cc
costfns.cc
generalio.cc
imfft.cc
lazy.cc
newimage.cc
newimagefns.cc
skeletonapp.cc)

INCLUDE_DIRECTORIES( ${FSL_SOURCE_DIR}/extras/src/libprob )
INCLUDE_DIRECTORIES( ${FSL_SOURCE_DIR}/extras/src/newmat )
ADD_LIBRARY(newimage STATIC ${NEWIMAGE_SOURCES})

 
install(DIRECTORY ./ DESTINATION ${FSL_SOURCE_DIR}/include/  FILES_MATCHING PATTERN "*.h")
install(TARGETS newimage DESTINATION ${FSL_SOURCE_DIR}/lib)
