/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#ifndef ___KMPCA__H
#define __KMPCA__H
#include <vector>
using namespace std;

namespace KML{

class CPrivacyPCA; 
class CPCA
{
private: 
  CPrivacyPCA* data;
public:
  CPCA();
  ~CPCA();
  void SetData(vector<vector<double> >& data); 
  void SetData(vector<vector<float> >& data); 
  /*****************************
  * -4, if SVD subroutine haven't converged
  *  1, if task is solved
  * ***************************/
  int DoPCA(vector<vector<double> >& basis,vector<double>& variance);
  int DoPCA(vector<vector<double> >& basis); 
  int DoPCA(vector<vector<float> >& basis,vector<float>& variance);
  int DoPCA(vector<vector<float> >& basis); 
};




};//end of namespace; 


#endif
