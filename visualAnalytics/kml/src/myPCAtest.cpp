/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#include "ap.h"
#include "dataanalysis.h"
#include <ctime>
#include <iostream>
#include "kmPCA.h"
#include "Vector3D.h"
using namespace std; 
using namespace KML; 
using namespace alglib_impl;

ostream& operator<< (ostream& out, ae_matrix& mat)
{
  for(int row=0; row< mat.rows; ++row)
  {
   for(int col=0; col< mat.cols; ++col)
   {
     out<<mat.ptr.pp_double[row][col]<<" "; 
   }
   cout<<endl;
  }
  return out; 
}
ostream& operator<< (ostream& out, ae_vector& vec)
{
  for(int row=0; row< vec.cnt; ++row)
  {
     out<<vec.ptr.p_double[row]<<" "; 
  }
  out<<endl;
  return out; 
}

//oriData: each line is a sample point; 
//basisVectorsInOrder;

ae_bool testpca( )
{
     ae_state _alglib_env_state;
   ae_state_init(&_alglib_env_state);
   ae_state *_state= &_alglib_env_state;
  
  
    ae_frame _frame_block;
    ae_int_t nDims;
    ae_int_t nSamplePoints;
    ae_int_t info;
    ae_vector means;
    ae_vector variances;
    ae_matrix basisVectorsInOrder;
    ae_matrix data;
 
    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&means, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&variances, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&basisVectorsInOrder, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&data, 0, 0, DT_REAL, _state, ae_true);

 
  
    nDims=3; 
    nSamplePoints=30;
 
  
  vector<vector<double> > rawData(nSamplePoints,vector<double>(nDims,0)); 
  vector<double> myVariance; 
  vector<vector<double> > basis; 
  vector<double> myMeans(nDims); 
      ae_matrix_set_length(&data, nSamplePoints, nDims, _state);
    ae_vector_set_length(&means, nDims, _state);

  for(int j=0; j<=nDims-1; j++)
    {
	myMeans[j] = 1.5*ae_randomreal(&_alglib_env_state)-0.75;
    }
    
    for(int i=0; i<=nSamplePoints-1; i++)
  {
	for(int j=0; j<=nDims-1; j++)
	{
	    rawData[i][j] = myMeans[j]+(2*ae_randomreal(&_alglib_env_state)-1);
	}
  }
  
  CPCA myTst; 
  myTst.SetData(rawData);
  myTst.DoPCA(basis,myVariance); 
  for(int index=0; index < basis.size(); ++index)
  cout<< basis[index]<<endl;
  
  cout<<myVariance<<endl; 
  

    for(int j=0; j<=nDims-1; j++)
    {
	means.ptr.p_double[j] = myMeans[j];
    }
    for(int i=0; i<=nSamplePoints-1; i++)
    {
	for(int j=0; j<=nDims-1; j++)
	{
	    data.ptr.pp_double[i][j] = rawData[i][j];
	}
    }

    cout<<data<<endl;

    pcabuildbasis(&data, nSamplePoints, nDims, &info, &variances, &basisVectorsInOrder, _state);
    if( info!=1 )
    {
	cout<<"Error in PCA Solution :"<<info<<endl;
    } 
    
    cout<<basisVectorsInOrder<<endl;
    cout<<variances<<endl;
    
    ae_frame_leave(_state);
}

int main()
{
  unsigned seed;
  time_t t;
  seed = (unsigned)time(&t);
  
  srand(seed);

  int nDims=3; 
  int nSamplePoints=30; 
  
  vector<vector<double> > rawData(nSamplePoints,vector<double>(nDims,0)); 
  vector<double> variances; 
  vector<vector<double> > basis; 
  vector<double> means(nDims); 
  
   ae_state _alglib_env_state;
    ae_state_init(&_alglib_env_state);
/*  
  for(int j=0; j<=nDims-1; j++)
    {
	means[j] = 1.5*ae_randomreal(&_alglib_env_state)-0.75;
    }
    
    for(int i=0; i<=nSamplePoints-1; i++)
  {
	for(int j=0; j<=nDims-1; j++)
	{
	    rawData[i][j] = means[j]+(2*ae_randomreal(&_alglib_env_state)-1);
	}
  }
  
  CPCA myTst; 
  myTst.SetData(rawData);
  myTst.DoPCA(basis,variances); 
  for(int index=0; index < basis.size(); ++index)
   cout<< basis[index]<<endl;
  
  cout<<variances<<endl; */
  
  
  testpca();
//   testpca(data,basisVectorsInOrder,variances,&_alglib_env_state);
 
	
  return 0; 
  
}
