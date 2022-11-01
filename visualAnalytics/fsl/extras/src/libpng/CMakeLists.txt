SET(PNG_SOURCES
png.c
pngerror.c
pnggccrd.c
pngget.c
pngmem.c
pngpread.c
pngread.c
pngrio.c
pngrtran.c
pngrutil.c
pngset.c
pngtrans.c
pngvcrd.c
pngwio.c
pngwrite.c
pngwtran.c
pngwutil.c)

INCLUDE_DIRECTORIES("C:/kmp/fsl/extras/src/zlib")
ADD_LIBRARY(libpng STATIC ${PNG_SOURCES})

IF(FSL_BUILD)
  INSTALL_FILES(/extras/include/libpng .*//.h$)
  INSTALL_TARGETS(/extras/lib libpng)
ENDIF(FSL_BUILD)
