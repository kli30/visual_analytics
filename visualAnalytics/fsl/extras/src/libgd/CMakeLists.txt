SET(LIBPROB_SOURCES
  bdtr.c btdtr.c chdtr.c drand.c fdtr.c gamma.c gdtr.c igam.c
  igami.c incbet.c incbi.c mtherr.c nbdtr.c ndtr.c ndtri.c pdtr.c
  stdtr.c unity.c polevl.c const.c)

ADD_LIBRARY(libgd STATIC ${LIBPROB_SOURCES})

install(DIRECTORY ./ DESTINATION C:/kmp/fsl/include/  FILES_MATCHING PATTERN "*.h")
install(TARGETS libgd DESTINATION C:/kmp/fsl/lib)