SET(NEWRAN_SOURCES
extreal.cpp
hist.cpp
newran.cpp)
INCLUDE_DIRECTORIES("${FSL_SOURCE_DIR}/extras/src/newmat")
ADD_LIBRARY(newran STATIC ${NEWRAN_SOURCES})



install(DIRECTORY ./ DESTINATION ${FSL_SOURCE_DIR}/include/  FILES_MATCHING PATTERN "*.h")
install(TARGETS newran DESTINATION ${FSL_SOURCE_DIR}/lib)
  
