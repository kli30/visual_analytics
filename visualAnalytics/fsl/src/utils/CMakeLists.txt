PROJECT(utils)

SET(utils_SOURCES
  matches.cc functions.cc usage.cc check.cc parse.cc)

ADD_LIBRARY(utils ${utils_SOURCES})

install(DIRECTORY ./ DESTINATION ${FSL_SOURCE_DIR}/include/  FILES_MATCHING PATTERN "*.h")
install(TARGETS utils DESTINATION ${FSL_SOURCE_DIR}/lib)
 
