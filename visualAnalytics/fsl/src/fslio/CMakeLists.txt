ADD_LIBRARY(fslio STATIC fslio.c)
include_directories(${FSL_SOURCE_DIR}/extras/src/zlib ${FSL_SOURCE_DIR}/src/znzlib)

install(DIRECTORY ./ DESTINATION ${FSL_SOURCE_DIR}/include/  FILES_MATCHING PATTERN "*.h")
install(TARGETS fslio DESTINATION ${FSL_SOURCE_DIR}/lib)
