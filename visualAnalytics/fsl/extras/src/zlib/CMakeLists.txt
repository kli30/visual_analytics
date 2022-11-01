SET(ZLIB_SOURCES
adler32.c
compress.c
crc32.c
deflate.c
example.c
gzio.c
infback.c
inffast.c
inflate.c
inftrees.c
minigzip.c
trees.c
uncompr.c
zutil.c)

ADD_LIBRARY(zlib STATIC ${ZLIB_SOURCES})


install(DIRECTORY ./ DESTINATION ${FSL_SOURCE_DIR}/include/  FILES_MATCHING PATTERN "*.h")
install(TARGETS zlib DESTINATION ${FSL_SOURCE_DIR}/lib)
