#ifndef PCACCA__H
#define PCACCA__H

#include "kaimingCommon.h"
#include "newimageall.h"

using namespace NEWIMAGE; 

namespace KML{

void principleComponentAnalysis(MatrixShen *var, float *enginValues, MatrixShen *enginVectors);

void kmPCA(volume4D<float>& bolds, volume<float>& roi, int label, vector<float>& eval, vector<vector<float> >& evec); 

void kmPCA(vector<vector<float> >& data,vector<float>& eval, vector<vector<float> >& evec ); 

void canonicalCorrelationAnalysis(MatrixShen *var1, MatrixShen *var2, MatrixShen *latentVar);

float pearsonCorrelation(vector <float> U, vector <float> V);

void normalizeVector(vector <float> &V);

}
#endif