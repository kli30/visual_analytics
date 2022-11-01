/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#pragma  once
#include <vector>
//modified by kaiming based on xiang li's code; 
using namespace std;
namespace KML{
vector<float> GCA(std::vector< float >& vector_X, std::vector< float >& vector_Y, int int_DiffFlag, float float_PVAL);
vector<float> GCA(vector<float>& vector_X, vector<float>& vector_Y, int int_DiffFlag);
}; //end of namespace kml;
