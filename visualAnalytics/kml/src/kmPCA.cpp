/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#include "kmPCA.h"
#include "ap.h"
#include "dataanalysis.h"
#include <ctime>
#include <iostream>

using namespace std; 
using namespace alglib_impl;


namespace KML{
  
  class CPrivacyPCA{
  public:
    CPrivacyPCA()
    {
      ae_state_init(&_alglib_env_state);
      ae_frame_make(&_alglib_env_state, &_frame_block);
      ae_vector_init(&variances, 0, DT_REAL, &_alglib_env_state, ae_true);
      ae_matrix_init(&basisVectorsInOrder, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
      ae_matrix_init(&rawData, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
      
	
    };
    ~CPrivacyPCA()
    {
      ae_frame_leave(&_alglib_env_state);
      ae_state_clear(&_alglib_env_state);
    };
    
    ae_frame _frame_block;
    ae_int_t nDims;
    ae_int_t nSamplePoints;
    ae_int_t info;
    ae_vector means;
    ae_vector variances;
    ae_matrix basisVectorsInOrder;
    ae_matrix rawData;
    ae_state _alglib_env_state;

  }; 
  
CPCA::CPCA()
{
  data = new CPrivacyPCA();
}
CPCA::~CPCA()
{
  delete data; 
}
int CPCA::DoPCA(vector< std::vector< double > >& basis, std::vector< double >& variance)
{
    pcabuildbasis(&(data->rawData), data->nSamplePoints, data->nDims, &(data->info), &(data->variances), &(data->basisVectorsInOrder), &(data->_alglib_env_state));
    basis.resize(data->basisVectorsInOrder.rows);
   
    for(int row=0; row< data->basisVectorsInOrder.rows; ++row)
    {
      basis[row].resize(data->basisVectorsInOrder.cols);
      for(int col=0; col< data->basisVectorsInOrder.cols; ++col)
      {
	basis[row][col]=data->basisVectorsInOrder.ptr.pp_double[row][col]; 
      }
    }
    
    variance.resize(data->variances.cnt);
    for(int row=0; row< data->variances.cnt; ++row)
    {
      variance[row]= data->variances.ptr.p_double[row]; 
    }
    
    return data->info; 
}

int CPCA::DoPCA(vector< std::vector< float > >& basis, std::vector< float >& variance)
{
    pcabuildbasis(&(data->rawData), data->nSamplePoints, data->nDims, &(data->info), &(data->variances), &(data->basisVectorsInOrder), &(data->_alglib_env_state));
    basis.resize(data->basisVectorsInOrder.rows);
   
    for(int row=0; row< data->basisVectorsInOrder.rows; ++row)
    {
      basis[row].resize(data->basisVectorsInOrder.cols);
      for(int col=0; col< data->basisVectorsInOrder.cols; ++col)
      {
	basis[row][col]=data->basisVectorsInOrder.ptr.pp_double[row][col]; 
      }
    }
    
    variance.resize(data->variances.cnt);
    for(int row=0; row< data->variances.cnt; ++row)
    {
      variance[row]= data->variances.ptr.p_double[row]; 
    }
    
    return data->info; 
}

int CPCA::DoPCA(vector< std::vector< double > >& basis)
{
  vector<double> variance; 
  return this->DoPCA(basis,variance);
}

int CPCA::DoPCA(vector< std::vector< float > >& basis)
{
  vector<float> variance; 
  return this->DoPCA(basis,variance);
}



void CPCA::SetData(vector< std::vector< double > >& inputData)
{
   
    data->nSamplePoints= inputData.size();
    data->nDims= inputData[0].size(); 
    
    ae_matrix_set_length(&(data->rawData), inputData.size(), inputData[0].size(), &(data->_alglib_env_state));
  
    for(int i=0; i<inputData.size(); i++)
    {
	for(int j=0; j<inputData[0].size(); j++)
	{
	    data->rawData.ptr.pp_double[i][j] = inputData[i][j];
	}
    }
}
void CPCA::SetData(vector< std::vector< float > >& inputData)
{
   
    data->nSamplePoints= inputData.size();
    data->nDims= inputData[0].size(); 
    
    ae_matrix_set_length(&(data->rawData), inputData.size(), inputData[0].size(), &(data->_alglib_env_state));
  
    for(int i=0; i<inputData.size(); i++)
    {
	for(int j=0; j<inputData[0].size(); j++)
	{
	    data->rawData.ptr.pp_double[i][j] = inputData[i][j];
	}
    }
}
   
}; //endl of namespace KML; 
