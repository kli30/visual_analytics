/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include "grangerCausality.h"

#include <vector>
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_fit.h"
#include "gsl/gsl_multifit.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_sf_log.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"

using namespace std;

namespace KML{

const int const_int_LAGMAX = 10;
const int const_int_LAGMIN = 1;
const float const_float_PVAL = 0.01;

#pragma region GCM Declare
void GCM_Regress(std::vector< float >& vector_X, std::vector< float >& vector_Y, int int_Lag, std::vector< float >& vector_Output);
void GCM_FindSignificance(std::vector< float >& vector_Prb, float float_Pval, float& float_Q, std::vector< float >& vector_PR);
int GCM_FindModelOrder(std::vector< float >& vector_X, std::vector< float >& vector_Y, int int_MinLag, int int_MaxLag);
float GCM_FindModelOrder_Regress(std::vector< float >& vector_X, std::vector< float >& vector_Y, int int_Lag);
vector<float> GCM_Diff(vector<float>& vector_Input, int int_DiffFlag);
#pragma endregion

#pragma region Stat Declare
float Stat_Covariance(std::vector< float >& vector_X, std::vector< float >& vector_Y);
float Stat_Mean(std::vector< float >& vector_Input);
float Stat_Variance(std::vector< float >& vector_Input);
float Stat_Stdev(std::vector< float >& vector_Input);
void Stat_Normalize(vector<float>& vector_Input);
#pragma endregion

#pragma region Stat
void Stat_Normalize(std::vector< float >& vector_Input)
{
	vector<float> vector_Result;

	int int_Row = vector_Input.size();
	float float_Mean = Stat_Mean(vector_Input);
	float float_Stdev = Stat_Stdev(vector_Input);
	for (int i = 0; i < int_Row; i++)
	{
		vector_Input[i]=((vector_Input[i] - float_Mean) / float_Stdev);
	}
}

float Stat_Covariance(vector<float>& vector_X, vector<float>& vector_Y)
{
	float float_Result = 0;

	float float_MeanX = Stat_Mean(vector_X);
	float float_MeanY = Stat_Mean(vector_Y);
	int int_Row = (int) vector_X.size();
	for (int i = 0; i < int_Row; i++)
	{
		float_Result += vector_X[i] * vector_Y[i];
	}
	float_Result = float_Result / int_Row;
	float_Result = float_Result - (float_MeanX * float_MeanY);

	return float_Result;
}

float Stat_Stdev(vector<float>& vector_Input)
{
	float float_Result = 0;

	float float_Mean = Stat_Mean(vector_Input);
	int int_Row = (int) vector_Input.size();
	for (int i = 0; i < int_Row; i++)
	{
		float_Result += (vector_Input[i] - float_Mean) * (vector_Input[i] - float_Mean);
	}
	float_Result = float_Result / (int_Row - 1);
	float_Result = sqrt(float_Result);

	return float_Result;
}

float Stat_Variance(vector<float>& vector_Input)
{
	float float_Result = 0;

	float float_Mean = Stat_Mean(vector_Input);
	int int_Row = (int) vector_Input.size();
	for (int i = 0; i < int_Row; i++)
	{
		float_Result += (vector_Input[i] - float_Mean) * (vector_Input[i] - float_Mean);
	}
	float_Result = float_Result / (int_Row - 1);

	return float_Result;
}

float Stat_Mean(vector<float>& vector_Input)
{
	float float_Result = 0;

	int int_Row = (int) vector_Input.size();
	for (int i = 0; i < int_Row; i++)
	{
		float_Result += vector_Input[i];
	}
	float_Result = float_Result / int_Row;

	return float_Result;
}
#pragma endregion

#pragma region GCM
vector<float> GCM_Diff(vector<float>& vector_Input, int int_DiffFlag)
{
	vector<float> vector_Result(vector_Input); 
	switch (int_DiffFlag)
	{
	case 0:
// 		for (int i = 0; i < (int)vector_Input.size(); i++)
// 		{
// 			vector_Result.push_back(vector_Input[i]);
// 		}
		break;
	case 1:
		for (int i = 1; i < (int)vector_Input.size(); i++)
		{
			vector_Result[i]=(vector_Input[i] - vector_Input[i - 1]);
		}
		break;
	case 2:
		for (int i = 2; i < (int)vector_Input.size(); i++)
		{
			vector_Result[i]=((vector_Input[i] - vector_Input[i - 1]) - (vector_Input[i - 1] - vector_Input[i - 2]));
		}
		break;
	default:
		cout<<"using 0 diffflag! "<<endl;
// 		cout << "Wrong arument provided in \"int_DiffFlag\"!";
// 		for (int i = 0; i < (int)vector_Input.size(); i++)
// 		{
// 			vector_Result.push_back(vector_Input[i]);
// 		}
	}

	return vector_Result;
}

int GCM_FindModelOrder(vector<float>& vector_X, vector<float>& vector_Y, int int_MinLag, int int_MaxLag)
{
	int int_Result = 0;
	float float_Z = 0;
	float float_Error = 0;
	int int_Nest = 0;
	float float_BIC = 0;
	int int_Row = (int) vector_X.size();
	gsl_vector * vec_BIC = gsl_vector_alloc(int_MaxLag - int_MinLag + 1);

	for (int i = int_MinLag; i <= int_MaxLag; i++)
	{
		float_Z = GCM_FindModelOrder_Regress(vector_X, vector_Y, i);
		if (float_Z > 0)
		{
			float_Error = gsl_sf_log(float_Z);
		}
		else
		{
			float_Error = 0;
		}
		int_Nest = 4 * i;
		float_BIC = float_Error + gsl_sf_log(int_Row) * int_Nest / int_Row;
		gsl_vector_set(vec_BIC, i - int_MinLag, float_BIC);
	}

	int_Result = (int)gsl_vector_min_index(vec_BIC) + int_MinLag;

	//cout<<"optimal model lag: "<<int_Result<<endl;
	return int_Result;
}

float GCM_FindModelOrder_Regress(vector<float>& vector_X, vector<float>& vector_Y, int int_Lag)
{
	float float_Result = 0;
	int int_Row = (int) vector_X.size();

	//REMOVE SAMPLE MEAN START
	float float_MeanX = Stat_Mean(vector_X);
	float float_MeanY = Stat_Mean(vector_Y);
	for (int i = 0; i < int_Row; i++)
	{
		vector_X[i] = vector_X[i] - float_MeanX;
		vector_Y[i] = vector_Y[i] - float_MeanY;
	}
	//REMOVE SAMPLE MEAN END

	//MLR START
	gsl_matrix * mtx_Regressor = gsl_matrix_alloc(int_Row - int_Lag, int_Lag * 2);
	for (int ii = 0; ii < int_Lag; ii++)
	{
		for (int z = 0; z < (int_Row  - int_Lag); z++)
		{
			gsl_matrix_set(mtx_Regressor, z, ii, vector_X[z - ii + int_Lag - 1]);
		}
	}
	for (int ii = 0; ii < int_Lag; ii++)
	{
		for (int z = 0; z < (int_Row  - int_Lag); z++)
		{
			gsl_matrix_set(mtx_Regressor, z, ii + int_Lag, vector_Y[z - ii + int_Lag - 1]);
		}
	}
	gsl_multifit_linear_workspace * linear_Work_Unrestricted = gsl_multifit_linear_alloc (int_Row - int_Lag, int_Lag * 2);
	float float_Residual = 0;
	vector<float> vector_C;

	//MLR -> X START
	gsl_vector * vec_Target_X = gsl_vector_alloc(int_Row - int_Lag);
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		gsl_vector_set(vec_Target_X, z, vector_X[z + int_Lag]);
	}
	gsl_vector * vec_Beta_X = gsl_vector_alloc(int_Lag * 2);
	gsl_matrix * mtx_Cov_X = gsl_matrix_alloc(int_Lag * 2, int_Lag * 2);
	double float_chi_X;
	gsl_multifit_linear (mtx_Regressor, vec_Target_X, vec_Beta_X, mtx_Cov_X, &float_chi_X, linear_Work_Unrestricted);
	vector<float> vector_Residual_X;
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		float_Residual = 0;
		for (int i = 0; i < (int_Lag * 2); i++)
		{
			float_Residual += gsl_matrix_get(mtx_Regressor, z, i) * gsl_vector_get(vec_Beta_X, i);
		}
		float_Residual = gsl_vector_get(vec_Target_X, z) - float_Residual;
		vector_Residual_X.push_back(float_Residual);
	}
	//MLR -> X END

	//MLR -> Y START
	gsl_vector * vec_Target_Y = gsl_vector_alloc(int_Row - int_Lag);
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		gsl_vector_set(vec_Target_Y, z, vector_Y[z + int_Lag]);
	}
	gsl_vector * vec_Beta_Y = gsl_vector_alloc(int_Lag * 2);
	gsl_matrix * mtx_Cov_Y = gsl_matrix_alloc(int_Lag * 2, int_Lag * 2);
	double float_chi_Y;
	gsl_multifit_linear (mtx_Regressor, vec_Target_Y, vec_Beta_Y, mtx_Cov_Y, &float_chi_Y, linear_Work_Unrestricted);
	vector<float> vector_Residual_Y;
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		float_Residual = 0;
		for (int i = 0; i < (int_Lag * 2); i++)
		{
			float_Residual += gsl_matrix_get(mtx_Regressor, z, i) * gsl_vector_get(vec_Beta_Y, i);
		}
		float_Residual = gsl_vector_get(vec_Target_Y, z) - float_Residual;
		vector_Residual_Y.push_back(float_Residual);
	}
	//MLR -> Y END
	//MLR END		

	gsl_multifit_linear_free (linear_Work_Unrestricted);
	gsl_vector_free(vec_Target_X);
	gsl_vector_free(vec_Beta_X);
	gsl_matrix_free(mtx_Cov_X);
	gsl_vector_free(vec_Target_Y);
	gsl_vector_free(vec_Beta_Y);
	gsl_matrix_free(mtx_Cov_Y);
	gsl_matrix_free(mtx_Regressor);
	float_Result = Stat_Covariance(vector_Residual_X, vector_Residual_X) * Stat_Covariance(vector_Residual_Y, vector_Residual_Y);
	float_Result = float_Result - (Stat_Covariance(vector_Residual_X, vector_Residual_Y) * Stat_Covariance(vector_Residual_Y, vector_Residual_X));

	return float_Result;
}

void GCM_FindSignificance(vector<float>& vector_Prb, float float_Pval, float &float_Q, vector<float> &vector_PR)
{
	float_Q = float_Pval / 2;
	if (vector_Prb[0] <= float_Q)
	{
		vector_PR.push_back(1);
	}
	else
	{
		vector_PR.push_back(0);
	}
	if (vector_Prb[1] <= float_Q)
	{
		vector_PR.push_back(1);
	}
	else
	{
		vector_PR.push_back(0);
	}
}

void GCM_Regress(vector<float>& vector_X, vector<float>& vector_Y, int int_Lag, vector<float> &vector_Output)
{
	int int_Row = (int) vector_X.size();

	//REMOVE SAMPLE MEAN START
	float float_MeanX = Stat_Mean(vector_X);
	float float_MeanY = Stat_Mean(vector_Y);
	for (int i = 0; i < int_Row; i++)
	{
		vector_X[i] = vector_X[i] - float_MeanX;
		vector_Y[i] = vector_Y[i] - float_MeanY;
	}
	//REMOVE SAMPLE MEAN END

	//UNRESTRICTED MLR START
	gsl_matrix * mtx_Regressor = gsl_matrix_alloc(int_Row - int_Lag, int_Lag * 2);
	for (int ii = 0; ii < int_Lag; ii++)
	{
		for (int z = 0; z < (int_Row  - int_Lag); z++)
		{
			gsl_matrix_set(mtx_Regressor, z, ii, vector_X[z - ii + int_Lag - 1]);
		}
	}
	for (int ii = 0; ii < int_Lag; ii++)
	{
		for (int z = 0; z < (int_Row  - int_Lag); z++)
		{
			gsl_matrix_set(mtx_Regressor, z, ii + int_Lag, vector_Y[z - ii + int_Lag - 1]);
		}
	}
	gsl_multifit_linear_workspace * linear_Work_Unrestricted = gsl_multifit_linear_alloc (int_Row - int_Lag, int_Lag * 2);
	float float_Residual = 0;
	float float_RSS = 0;
	vector<float> vector_RSSU;
	vector<float> vector_C;

	//MLR -> X START
	gsl_vector * vec_Target_X = gsl_vector_alloc(int_Row - int_Lag);
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		gsl_vector_set(vec_Target_X, z, vector_X[z + int_Lag]);
	}
	gsl_vector * vec_Beta_X = gsl_vector_alloc(int_Lag * 2);
	gsl_matrix * mtx_Cov_X = gsl_matrix_alloc(int_Lag * 2, int_Lag * 2);
	double float_chi_X;
	gsl_multifit_linear (mtx_Regressor, vec_Target_X, vec_Beta_X, mtx_Cov_X, &float_chi_X, linear_Work_Unrestricted);
	vector<float> vector_Residual_X;
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		float_Residual = 0;
		for (int i = 0; i < (int_Lag * 2); i++)
		{
			float_Residual += gsl_matrix_get(mtx_Regressor, z, i) * gsl_vector_get(vec_Beta_X, i);
		}
		float_Residual = gsl_vector_get(vec_Target_X, z) - float_Residual;
		vector_Residual_X.push_back(float_Residual);
	}
	float_RSS = 0;
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		float_RSS += vector_Residual_X[z] * vector_Residual_X[z];
	}
	vector_RSSU.push_back(float_RSS);
	vector_C.push_back(Stat_Covariance(vector_Residual_X, vector_Residual_X));
	//MLR -> X END

	//MLR -> Y START
	gsl_vector * vec_Target_Y = gsl_vector_alloc(int_Row - int_Lag);
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		gsl_vector_set(vec_Target_Y, z, vector_Y[z + int_Lag]);
	}
	gsl_vector * vec_Beta_Y = gsl_vector_alloc(int_Lag * 2);
	gsl_matrix * mtx_Cov_Y = gsl_matrix_alloc(int_Lag * 2, int_Lag * 2);
	double float_chi_Y;
	gsl_multifit_linear (mtx_Regressor, vec_Target_Y, vec_Beta_Y, mtx_Cov_Y, &float_chi_Y, linear_Work_Unrestricted);
	vector<float> vector_Residual_Y;
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		float_Residual = 0;
		for (int i = 0; i < (int_Lag * 2); i++)
		{
			float_Residual += gsl_matrix_get(mtx_Regressor, z, i) * gsl_vector_get(vec_Beta_Y, i);
		}
		float_Residual = gsl_vector_get(vec_Target_Y, z) - float_Residual;
		vector_Residual_Y.push_back(float_Residual);
	}
	float_RSS = 0;
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		float_RSS += vector_Residual_Y[z] * vector_Residual_Y[z];
	}
	vector_RSSU.push_back(float_RSS);	
	vector_C.push_back(Stat_Covariance(vector_Residual_Y, vector_Residual_Y));
	//MLR -> Y END

	gsl_multifit_linear_free (linear_Work_Unrestricted);
	gsl_vector_free(vec_Target_X);
	gsl_vector_free(vec_Target_Y);
	gsl_vector_free(vec_Beta_X);
	gsl_matrix_free(mtx_Cov_X);
	gsl_vector_free(vec_Beta_Y);
	gsl_matrix_free(mtx_Cov_Y);
	gsl_matrix_free(mtx_Regressor);
	vector_Residual_X.clear();
	float_chi_X = 0;
	vector_Residual_Y.clear();
	float_chi_Y = 0;
	//UNRESTRICTED MLR END


	//RESTRICTED MLR START
	gsl_multifit_linear_workspace * linear_Work_Restricted = gsl_multifit_linear_alloc (int_Row - int_Lag, int_Lag);
	float_Residual = 0;
	float_RSS = 0;
	vector<float> vector_RSSS;
	vector<float> vector_S;
	gsl_matrix * mtx_Regressor_X = gsl_matrix_alloc (int_Row - int_Lag, int_Lag);
	gsl_matrix * mtx_Regressor_Y = gsl_matrix_alloc (int_Row - int_Lag, int_Lag);

	//MLR -> X START
	vec_Target_X = gsl_vector_alloc(int_Row - int_Lag);
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		gsl_vector_set(vec_Target_X, z, vector_X[z + int_Lag]);
	}
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		for (int ii = 0; ii < int_Lag; ii++)
		{
			gsl_matrix_set(mtx_Regressor_X, z, ii, vector_X[z + int_Lag - ii - 1]);
		}
	}
	vec_Beta_X = gsl_vector_alloc(int_Lag);
	mtx_Cov_X = gsl_matrix_alloc(int_Lag, int_Lag);
	gsl_multifit_linear (mtx_Regressor_X, vec_Target_X, vec_Beta_X, mtx_Cov_X, &float_chi_X, linear_Work_Restricted);
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		float_Residual = 0;
		for (int i = 0; i < int_Lag; i++)
		{
			float_Residual += gsl_matrix_get(mtx_Regressor_X, z, i) * gsl_vector_get(vec_Beta_X, i);
		}
		float_Residual = gsl_vector_get(vec_Target_X, z) - float_Residual;
		vector_Residual_X.push_back(float_Residual);
	}
	vector_S.push_back(Stat_Covariance(vector_Residual_X, vector_Residual_X));
	float_RSS = 0;
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		float_RSS += vector_Residual_X[z] * vector_Residual_X[z];
	}
	vector_RSSS.push_back(float_RSS);
	//MLR -> X END

	//MLR -> Y START
	vec_Target_Y = gsl_vector_alloc(int_Row - int_Lag);
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		gsl_vector_set(vec_Target_Y, z, vector_Y[z + int_Lag]);
	}
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		for (int ii = 0; ii < int_Lag; ii++)
		{
			gsl_matrix_set(mtx_Regressor_Y, z, ii, vector_Y[z + int_Lag - ii - 1]);
		}
	}
	vec_Beta_Y = gsl_vector_alloc(int_Lag);
	mtx_Cov_Y = gsl_matrix_alloc(int_Lag, int_Lag);
	gsl_multifit_linear (mtx_Regressor_Y, vec_Target_Y, vec_Beta_Y, mtx_Cov_Y, &float_chi_Y, linear_Work_Restricted);
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		float_Residual = 0;
		for (int i = 0; i < int_Lag; i++)
		{
			float_Residual += gsl_matrix_get(mtx_Regressor_Y, z, i) * gsl_vector_get(vec_Beta_Y, i);
		}
		float_Residual = gsl_vector_get(vec_Target_Y, z) - float_Residual;
		vector_Residual_Y.push_back(float_Residual);
	}
	vector_S.push_back(Stat_Covariance(vector_Residual_Y, vector_Residual_Y));
	float_RSS = 0;
	for (int z = 0; z < (int_Row - int_Lag); z++)
	{
		float_RSS += vector_Residual_Y[z] * vector_Residual_Y[z];
	}
	vector_RSSS.push_back(float_RSS);
	//MLR -> Y END

	gsl_multifit_linear_free (linear_Work_Restricted);
	gsl_vector_free(vec_Target_X);
	gsl_vector_free(vec_Target_Y);
	gsl_vector_free(vec_Beta_X);
	gsl_matrix_free(mtx_Cov_X);
	gsl_vector_free(vec_Beta_Y);
	gsl_matrix_free(mtx_Cov_Y);
	gsl_matrix_free(mtx_Regressor_X);
	gsl_matrix_free(mtx_Regressor_Y);
	vector_Residual_X.clear();
	vector_Residual_Y.clear();
	//RESTRICTED MLR END


	//GRANGER VALUE CALCULATION START
	int int_N2 = int_Row - int_Lag - (2 * int_Lag);

	vector<float> vector_Ftest;
	float float_Ftest = 0;
	float_Ftest = ((vector_RSSS[0] - vector_RSSU[0]) / int_Lag) / (vector_RSSU[0] / int_N2);
	vector_Ftest.push_back(float_Ftest);
	float_Ftest = ((vector_RSSS[1] - vector_RSSU[1]) / int_Lag) / (vector_RSSU[1] / int_N2);
	vector_Ftest.push_back(float_Ftest);

	vector<float> vector_Prb;
	float float_Prb = 0;
	float_Prb = 1 - gsl_cdf_fdist_P(vector_Ftest[0], int_Lag, int_N2);

//	cout<<"lower cumulative probability x-> y :"<<1-float_Prb<<endl;
//	cout<<"upper cumulative probability x-> y :"<<float_Prb <<endl;
	vector_Prb.push_back(float_Prb);
	float_Prb = 1 - gsl_cdf_fdist_P(vector_Ftest[1], int_Lag, int_N2);
	vector_Prb.push_back(float_Prb);
//	cout<<"lower cumulative probability y-> x :"<<1-float_Prb<<endl;
//	cout<<"upper cumulative probability y-> x :"<< float_Prb<<endl;


	vector<float> vector_GC;
	float float_GC = 0;
	float_GC = gsl_sf_log(vector_S[0] / vector_C[0]);
	vector_GC.push_back(float_GC);
	float_GC = gsl_sf_log(vector_S[1] / vector_C[1]);
	vector_GC.push_back(float_GC);

	vector<float> vector_DOI;
	float float_DOI = 0;
	float_DOI = vector_GC[0] - vector_GC[1];
	vector_DOI.push_back(float_DOI);
	float_DOI = vector_GC[1] - vector_GC[0];
	vector_DOI.push_back(float_DOI);
	//GRANGER VALUE CALCULATION END


	//OUTPUT START
	vector_Output.clear();
	vector_Output.push_back(vector_Prb[0]);
	vector_Output.push_back(vector_Prb[1]);
	vector_Output.push_back(vector_GC[0]);
	vector_Output.push_back(vector_GC[1]);
	vector_Output.push_back(vector_Ftest[0]);
	vector_Output.push_back(vector_Ftest[1]);
	vector_Output.push_back(vector_DOI[0]);
	vector_Output.push_back(vector_DOI[1]);
	//OUTPUT END
}
#pragma endregion

vector<float> GCA(vector<float>& vector_X,  vector<float>& vector_Y, int int_DiffFlag, float float_PVAL)
{
	vector<float> vector_X_Diff;
	vector<float> vector_Y_Diff;
	vector<float> vector_Output;
	int int_BIC = 0;

	vector_Output.clear();
	vector_X_Diff = GCM_Diff(vector_X, int_DiffFlag);
	vector_Y_Diff = GCM_Diff(vector_Y, int_DiffFlag);
	int_BIC = GCM_FindModelOrder(vector_X_Diff, vector_Y_Diff, const_int_LAGMIN, const_int_LAGMAX);
	GCM_Regress(vector_X_Diff, vector_Y_Diff, int_BIC, vector_Output);
	vector<float> vector_Prb;
	vector_Prb.push_back(vector_Output[0]);
	vector_Prb.push_back(vector_Output[1]);
	
	vector_Prb[0]=vector_Prb[0]>1-vector_Prb[0]? 1-vector_Prb[0] :  vector_Prb[0];
	vector_Prb[1]=vector_Prb[1]> 1-vector_Prb[1] ? 1-vector_Prb[1] : vector_Prb[1]; 

	float pVal=0.5*float_PVAL;
	vector_Prb[0] =vector_Prb[0] < pVal ? 1 : 0;
	vector_Prb[1]= vector_Prb[1] <pVal ? 1 : 0; 
	return vector_Prb;
}

vector<float> GCA(vector<float>& vector_X, vector<float>& vector_Y, int int_DiffFlag)
{
	vector<float> vector_X_Diff;
	vector<float> vector_Y_Diff;
	vector<float> vector_Output;
	int int_BIC = 0;

	vector_Output.clear();
	vector_X_Diff = GCM_Diff(vector_X, int_DiffFlag);
	vector_Y_Diff = GCM_Diff(vector_Y, int_DiffFlag);
	int_BIC = GCM_FindModelOrder(vector_X_Diff, vector_Y_Diff, const_int_LAGMIN, const_int_LAGMAX);
	GCM_Regress(vector_X_Diff, vector_Y_Diff, int_BIC, vector_Output);
	return vector_Output;
}

} //end of KML;
