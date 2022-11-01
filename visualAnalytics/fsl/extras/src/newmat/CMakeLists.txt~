SET(NEWMAT_SOURCES
	bandmat.cpp cholesky.cpp evalue.cpp fft.cpp hholder.cpp
	jacobi.cpp myexcept.cpp newmat1.cpp newmat2.cpp newmat3.cpp
	newmat4.cpp newmat5.cpp newmat6.cpp newmat7.cpp newmat8.cpp
	newmat9.cpp newmatex.cpp newmatnl.cpp newmatrm.cpp solution.cpp
	sort.cpp submat.cpp svd.cpp newfft.cpp)

ADD_LIBRARY(newmat STATIC ${NEWMAT_SOURCES})

install(DIRECTORY ./ DESTINATION ${FSL_SOURCE_DIR}//include/  FILES_MATCHING PATTERN "*.h")
install(TARGETS newmat DESTINATION ${FSL_SOURCE_DIR}//lib)