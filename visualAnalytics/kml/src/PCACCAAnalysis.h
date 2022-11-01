/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/


#ifndef __PCACCA__H
#define __PCACCA__H

#include "kaimingCommon.h"

//edit by kaiming based on xintao's code; 
namespace KML{

void principleComponentAnalysis(MatrixShen *var, float *enginValues, MatrixShen *enginVectors);

void canonicalCorrelationAnalysis(MatrixShen *var1, MatrixShen *var2, MatrixShen *latentVar);

float pearsonCorrelation(int i, MatrixShen *U, int j, MatrixShen *V);

}

#endif
