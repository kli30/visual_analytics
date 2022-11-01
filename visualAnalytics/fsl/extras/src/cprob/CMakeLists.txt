SET(LIBPROB_SOURCES
  bdtr.c btdtr.c chdtr.c drand.c fdtr.c gamma.c gdtr.c igam.c
  igami.c incbet.c incbi.c mtherr.c nbdtr.c ndtr.c ndtri.c pdtr.c
  stdtr.c unity.c polevl.c const.c)

ADD_LIBRARY(prob STATIC ${LIBPROB_SOURCES})

install(DIRECTORY ./ DESTINATION ${FSL_SOURCE_DIR}/include/  FILES_MATCHING PATTERN "*.h")
install(TARGETS prob DESTINATION ${FSL_SOURCE_DIR}/lib)
