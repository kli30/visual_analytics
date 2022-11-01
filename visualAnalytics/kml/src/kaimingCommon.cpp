/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/


#include "kaimingCommon.h"
#include "newmat.h"
#include <sstream>
#include <vtkArrowSource.h>
#include <vtkConeSource.h>
#include <vtkCylinderSource.h>
#include <vtkSphereSource.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkActor.h>
#include <vtkTransform.h>
#include "vtkTransformPolyDataFilter.h"
#include <vtkCell.h>
#include <list>
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include "ColorSchemeRGB.h"
#include "cfloat"

using namespace NEWMAT;
using namespace std;

namespace KML{

void GetOffsetFromMHD( const char* fileName, Vector3D<float>& offset)
{
	fstream headerStrm;
	OpenReadStrmAscii( headerStrm, fileName);
	string line;
	strstream lineStrm;
	do
	{
		getline(headerStrm,line);
	} while (line.find("Offset")==string::npos);
	lineStrm.clear();
	lineStrm<<line;
	lineStrm>>line>>line; //remove "Offset = "
	lineStrm>>offset.x>>offset.y>>offset.z;

	headerStrm.close();

}

void GetDimsFromMHD( const char* fileName, Vector3D<float>& dims){

	fstream headerStrm;
	OpenReadStrmAscii( headerStrm, fileName);
	string line;
	strstream lineStrm;
	do
	{
		getline(headerStrm,line);
	} while (line.find("ElementSpacing")==string::npos);
	lineStrm.clear();
	lineStrm<<line;
	lineStrm>>line>>line; //remove "Offset = "
	lineStrm>>dims.x>>dims.y>>dims.z;

	headerStrm.close();
}

void GetSizesFromMHD(const char* fileName, Vector3D<size_t>& sizes)
{
	fstream headerStrm;
	OpenReadStrmAscii( headerStrm, fileName);
	string line;
	strstream lineStrm;
	do
	{
		getline(headerStrm,line);
	} while (line.find("DimSize")==string::npos);
	lineStrm.clear();
	lineStrm<<line;
	lineStrm>>line>>line; //remove "Offset = "
	lineStrm>>sizes.x>>sizes.y>>sizes.z;

	headerStrm.close();
}

void  GetSizesFromMHD(const char* fileName,  Vector3D< float >& sizes)
{
  Vector3D<size_t> sizeInt;
  GetSizesFromMHD(fileName, sizeInt);
  sizes.x=sizeInt.x;
  sizes.y=sizeInt.y; 
  sizes.z=sizeInt.z; 
}

void GetMetaInfoFromMHD(string fileName, Vector3D< float >& offset, Vector3D< float >& dims, Vector3D< float >& sizes)
{

}

void GetImgGridFromPhysicalCoord( const Vector3D<float>& phisicalCoord, Vector3D<size_t>&  imageCoord, const Vector3D<float>& offSet, const Vector3D<float>& dims){
    imageCoord.x = size_t((phisicalCoord.x - offSet.x) / dims.x + 0.5);
    imageCoord.y = size_t((phisicalCoord.y - offSet.y) / dims.y + 0.5);
    imageCoord.z = size_t((phisicalCoord.z - offSet.z) / dims.z + 0.5);
}

void TransposeIntMatrix(string infile,string outFile)
{
  vector<vector<int> > mat,matT;  
  KML::ReadIntMatrix(infile,mat);
  KML::TransposeMatrix<>(mat,matT);
  
  KML::SaveIntMatrix(outFile,matT);  
}

void TransposeFloatMatrix(string infile, string outFile)
{
  vector<vector<float> > mat,matT; 
  KML::ReadFloatMatrix(infile,mat);
  KML::TransposeMatrix<>(mat,matT);
  KML::SaveFloatMatrix(outFile,matT); 

}


////////////////////////////////////
//file operations;
void OpenReadStrmAscii(fstream& strm, const char* fileName){
	strm.open(fileName,ios::in);
	if(NULL==strm){
		cerr<<"Error while opening file :"<<fileName<< endl;
		exit(1);
	}
}

void OpenWriteStrmAscii(fstream& strm, const char* fileName){
	strm.open(fileName,ios::out);
	if(NULL==strm){
		cerr<<"Error while opening file :"<<fileName<< endl;
		exit(1);
	}
}

void OpenReadStrmBinary(fstream& strm, const char* fileName){
	strm.open(fileName,ios::in|ios::binary);
	if(NULL==strm){
		cerr<<"Error while opening file :"<<fileName<< endl;
		exit(1);
	}
}

void OpenWriteStrmBinary(fstream& strm, const char*  fileName){
	strm.open(fileName,ios::binary|ios::out);
	if(NULL==strm){
		cerr<<"Error while opening file :"<<fileName<< endl;
		exit(1);
	}
}

void OpenReadStrmAscii(fstream& strm, const string& fileName){
	strm.open(fileName.c_str(),ios::in);
	if(NULL==strm){
		cerr<<"Error while opening file :"<<fileName<< endl;
		exit(1);
	}
}

void OpenWriteStrmAscii(fstream& strm, const string& fileName){
	strm.open(fileName.c_str(),ios::out);
	if(NULL==strm){
		cerr<<"Error while opening file :"<<fileName<< endl;
		exit(1);
	}
}

void OpenReadStrmBinary(fstream& strm, const string& fileName){
	strm.open(fileName.c_str(),ios::in|ios::binary);
	if(NULL==strm){
		cerr<<"Error while opening file :"<<fileName<< endl;
		exit(1);
	}
}

void OpenWriteStrmBinary(fstream& strm, const string& fileName){
	strm.open(fileName.c_str(),ios::binary|ios::out);
	if(NULL==strm){
		cerr<<"Error while opening file :"<<fileName<< endl;
		exit(1);
	}
}

///file operation ends;

void StatsNormalizeHist(std::vector< float >& vector_Input)
{
  float sum_= FLT_MIN; 
  for(int idx = 0; idx < vector_Input.size(); ++idx )
  {  
    sum_+=vector_Input[idx];
  } // end for::idx;
  
  for(int idx = 0; idx < vector_Input.size(); ++idx )
  {  
    vector_Input[idx]/=sum_;
  } // end for::idx;
}
float StatsSum(std::vector< float >& vector_Input)
{
  float sum_=0; 
  for(int idx = 0; idx < vector_Input.size(); ++idx )
  {  
    sum_+=vector_Input[idx];
  } // end for::idx;
  return sum_; 
}

///stats starts
void StatsNormalizeFisher(std::vector< float >& vector_Input)
{
        vector<float> vector_Result;

        int int_Row = vector_Input.size();
        float float_Mean = StatsMean(vector_Input);
        float float_Stdev = StatsStdev(vector_Input);
        float_Stdev+=FLT_MIN ;  
        for (int i = 0; i < int_Row; i++)
        {
                vector_Input[i]=((vector_Input[i] - float_Mean) / float_Stdev);
        }
}

float StatsCovariance(vector<float>& vector_X, vector<float>& vector_Y)
{
        float float_Result = 0;

        float float_MeanX = StatsMean(vector_X);
        float float_MeanY = StatsMean(vector_Y);
        int int_Row = (int) vector_X.size();
        for (int i = 0; i < int_Row; i++)
        {
                float_Result += vector_X[i] * vector_Y[i];
        }
        float_Result = float_Result / int_Row;
        float_Result = float_Result - (float_MeanX * float_MeanY);

        return float_Result;
}

float StatsStdev(vector<float>& vector_Input)
{
        float float_Result = 0;

        float float_Mean = StatsMean(vector_Input);
        int int_Row = (int) vector_Input.size();
        for (int i = 0; i < int_Row; i++)
        {
                float_Result += (vector_Input[i] - float_Mean) * (vector_Input[i] - float_Mean);
        }
        float_Result = float_Result / (int_Row - 1);
        float_Result = sqrt(float_Result);

        return float_Result;
}

float StatsVariance(vector<float>& vector_Input)
{
        float float_Result = 0;

        float float_Mean = StatsMean(vector_Input);
        int int_Row = (int) vector_Input.size();
        for (int i = 0; i < int_Row; i++)
        {
                float_Result += (vector_Input[i] - float_Mean) * (vector_Input[i] - float_Mean);
        }
        float_Result = float_Result / (int_Row - 1);

        return float_Result;
}

float StatsMean(vector<float>& vector_Input)
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


void StatsTTestEqualVariance(vector<float>& group1, vector<float>& group2,double& pvalueG1Bigger, double& pvalueG2Bigger)
{
    float meanG1 = KML::StatsMean(group1);
    float meanG2 = KML::StatsMean(group2);
    float stdevG1 = KML::StatsStdev(group1);
    float stdevG2 = KML::StatsStdev(group2);
    
    float meanDiff = meanG1-meanG2; 
    double degree = group1.size()+group2.size()-2; 
    float S12 = sqrt( ((group1.size()-1)*stdevG1*stdevG1+(group2.size()-1)*stdevG2*stdevG2)/degree);
    double tValue = meanDiff*sqrt(group1.size()*group2.size())/S12/sqrt(group2.size()+group1.size());
    
    pvalueG1Bigger = gsl_cdf_tdist_P(tValue,degree);
    pvalueG2Bigger = gsl_cdf_tdist_Q(tValue,degree);      
}

///stats ends
//////////////////////////////////////////////////////////
///swap endians
bool IsBigEndian()
{
   short word = 0x4321;
   if((*(char *)& word) != 0x21 )
     return true;
   else
     return false;
}

void ByteSwap(unsigned char * b, int n)
{
   register int i = 0;
   register int j = n-1;
   while (i<j)
   {
      std::swap(b[i], b[j]);
      i++, j--;
   }
}
///swap endians end;

float FindMaxValueMatrix(Matrix & matrix, int& line, int& col){
	float max=-1e30;
	for (int lineIndex = 0; lineIndex < matrix.Nrows(); ++lineIndex) {
		for (int colIndex = 0; colIndex < matrix.Ncols(); ++colIndex) {
			if( matrix(lineIndex+1,colIndex+1)>max){
				line=lineIndex+1;
				col=colIndex+1;
				max=matrix(lineIndex+1,colIndex+1);
			}
		} //end of for loop:: colIndex
	} //end of for loop:: lineIndex
	return max;
}

float FindMinValueMatrix(Matrix & matrix, int& line, int& col){
	float min=1e30;
	for (int lineIndex = 0; lineIndex < matrix.Nrows(); ++lineIndex) {
		for (int colIndex = 0; colIndex < matrix.Ncols(); ++colIndex) {
			if( matrix(lineIndex+1,colIndex+1)< min){
				line=lineIndex+1;
				col=colIndex+1;
				min=matrix(lineIndex+1,colIndex+1);
			}
		} //end of for loop:: colIndex
	} //end of for loop:: lineIndex
	return min;
}

string operator + ( string base, int index){
	string newString;
	stringstream	 tmpStrm;
	tmpStrm<<index;
	tmpStrm>>newString;
	return base+newString;
}

string&  operator += ( string&  base, int index){
  string tmpString =base; 
  base = tmpString+ index; 
  return base;
}


void ReadIntMatrix(string ipmiGroup, vector<vector<int> >& mat)
{
  fstream inStrm; 
  KML::OpenReadStrmAscii(inStrm,ipmiGroup);
  string tmpString; 
  while(getline(inStrm,tmpString))
  {
    if(tmpString.size())
    {
      vector<int> tmpVec; 
      KML::ReadVectorFromString<>(tmpVec,tmpString);
      mat.push_back(tmpVec);      
    }
  }
  inStrm.close();
}
void ReadFloatMatrix(string fileName, vector< std::vector< float > >& mat)
{
  fstream inStrm; 
  KML::OpenReadStrmAscii(inStrm,fileName);
  string tmpString; 
  while(getline(inStrm,tmpString))
  {
    if(tmpString.size())
    {
      vector<float> tmpVec; 
      KML::ReadVectorFromString<>(tmpVec,tmpString);
      mat.push_back(tmpVec);      
    }
  }
  inStrm.close();
}

void SaveIntMatrix(string outFileName, vector< std::vector< int > >& mat)
{
  fstream outStrm; 
  KML::OpenWriteStrmAscii(outStrm,outFileName);
  for(int idx1 = 0; idx1 < mat.size(); ++idx1 )
  {  
    for(int idx2 = 0; idx2 < mat[idx1].size(); ++idx2 )
    {  
      outStrm<<mat[idx1][idx2]<<" ";
    } // end for::idx2;
    outStrm<<endl;
  } // end for::idx1;
 
  outStrm.close();
}

void SaveFloatMatrix(string fileName, vector< std::vector< float > >& mat)
{
  fstream outStrm; 
  KML::OpenWriteStrmAscii(outStrm,fileName);
  for(int idx1 = 0; idx1 < mat.size(); ++idx1 )
  {  
    for(int idx2 = 0; idx2 < mat[idx1].size(); ++idx2 )
    {  
      outStrm<<mat[idx1][idx2]<<" ";
    } // end for::idx2;
    outStrm<<endl;
  } // end for::idx1;
 
  outStrm.close();
}

void ReadNameList(string fileName, vector<string > & allNames)
{
  fstream inStrm; 
  KML::OpenReadStrmAscii(inStrm,fileName);
  string tmpString; 
  while(inStrm>>tmpString)
  {
    if(tmpString.size())
    {
      allNames.push_back(tmpString);      
    }
  }
  inStrm.close();
}

void SaveNameList(string fileName, vector< string >& allNames)
{
  fstream outStrm; 
  KML::OpenWriteStrmAscii(outStrm,fileName);
  for(int idx = 0; idx < allNames.size(); ++idx )
  {  
    outStrm<<allNames[idx]<<endl;
  } // end for::idx;

  outStrm.close();
}
void VisualizePointsUsingSphereByVTKWithColors(vector< KML::Vector3D< float > >& roiCenters, vector< KML::Vector3D< float > >& allPointColors, string fileName, float radius, float resolution)
{
        vector<vtkPolyData*> allSpheresData(roiCenters.size());
        vector<vtkSphereSource*> allSphereSource(roiCenters.size());

        double sphereRadius= radius;
        int sphereResolution=resolution;
        int numCenters=roiCenters.size();

        for (int roiIndex = 0 ; roiIndex < numCenters; ++roiIndex) {
                //for each roi center , generate a sphere centering at this center;
                allSphereSource[roiIndex]=vtkSphereSource::New();
                allSpheresData[roiIndex]=vtkPolyData::New();
                allSphereSource[roiIndex]->SetThetaResolution(sphereResolution);
                allSphereSource[roiIndex]->SetPhiResolution(sphereResolution);
                allSphereSource[roiIndex]->SetCenter(roiCenters[roiIndex].x,roiCenters[roiIndex].y,roiCenters[roiIndex].z);
                allSphereSource[roiIndex]->SetRadius(sphereRadius);
                allSpheresData[roiIndex]= allSphereSource[roiIndex]->GetOutput();
                allSpheresData[roiIndex]->Update();
        } //end of for loop:: roiIndex

        vector<Vector3D<float> > allInOnePoints;
        vector<list<int> > allInOneCells;
        int numCells=0;
        int allCellSize=0;
        
        
        vector<Vector3D<float> > allInOnePointsColors;


        {
                //put all spheres  together;
                for (int centerIndex = 0 ; centerIndex < numCenters; ++centerIndex) {
                        //points;
                        for (int pointId = 0 ; pointId < allSpheresData[centerIndex]->GetNumberOfPoints(); ++pointId) {
                                Vector3D<float> pointCoord;
                                pointCoord.x=allSpheresData[centerIndex]->GetPoint(pointId)[0];
                                pointCoord.y=allSpheresData[centerIndex]->GetPoint(pointId)[1];
                                pointCoord.z=allSpheresData[centerIndex]->GetPoint(pointId)[2];
                                allInOnePoints.push_back(pointCoord);
                                allInOnePointsColors.push_back(allPointColors[centerIndex]);
                        } //end of for loop:: pointId

                        //cells;
                        for (int cellId = 0; cellId < allSpheresData[centerIndex]->GetNumberOfCells(); ++cellId) {
                                list<int> newCell;
                                vtkCell* currentCell= allSpheresData[centerIndex]->GetCell(cellId);
                                for (int cellElem = 0 ; cellElem < currentCell->GetNumberOfPoints(); ++cellElem) {
                                        int pointidIncurrentSphere= currentCell->GetPointId(cellElem);
                                        int pointIdInAll=pointidIncurrentSphere+allSpheresData[centerIndex]->GetNumberOfPoints()*centerIndex;
                                        newCell.push_back(pointIdInAll);
                                } //end of for loop:: cellElem
                                allInOneCells.push_back(newCell);
                                allCellSize+=newCell.size();
                                allCellSize+=1;
                        } //end of for loop:: cellId

                } //end of for loop:: centerIndex

        }

        // WriteOutAllInVTK(allInOnePoints,allInOneCells,allCellSize,subjectROIname);
        cout<<"write all objects in one vtk file: "<<fileName<<endl;
        fstream vtkStrm;
        OpenWriteStrmAscii(vtkStrm,fileName.c_str());
        vtkStrm<<"# vtk DataFile Version 3.0\n"<<"vtk output\n"<<"ASCII\n"<<"DATASET POLYDATA"<<endl;
        vtkStrm<<"POINTS "<<allInOnePoints.size()<<" FLOAT"<<endl;
        for (int pointId = 0 ; pointId < allInOnePoints.size(); ++pointId) {
                vtkStrm<<allInOnePoints[pointId]<<endl;
        } //end of for loop:: pointId

        vtkStrm<<"POLYGONS "<<allInOneCells.size()<<" "<<allCellSize<<endl;

        for (int cellId = 0; cellId < allInOneCells.size(); ++cellId) {
                vtkStrm<<allInOneCells[cellId].size()<<" ";
                for(list<int>::iterator itCellElems=allInOneCells[cellId].begin(); itCellElems!=allInOneCells[cellId].end(); itCellElems++)
                        vtkStrm<<*itCellElems<<" ";
                vtkStrm<<endl;
        } //end of for loop:: cellId
        
        
        vtkStrm << "POINT_DATA " << allInOnePoints.size() << endl;
        vtkStrm << "COLOR_SCALARS " << "PointColor  3 " << endl;
        for (int idx = 0; idx < allInOnePoints.size() ; ++idx) {
            vtkStrm<<allInOnePointsColors[idx]<<endl;
        }//end of for loop::idx
        
}

void VisualizePointsUsingSphereByVTK(vector< KML::Vector3D< float > >& roiCenters, string fileName, float radius, float resolution ){

	vector<vtkPolyData*> allSpheresData(roiCenters.size());
	vector<vtkSphereSource*> allSphereSource(roiCenters.size());

	double sphereRadius= radius;
	int sphereResolution=resolution;
	int numCenters=roiCenters.size();

	for (int roiIndex = 0 ; roiIndex < numCenters; ++roiIndex) {
		//for each roi center , generate a sphere centering at this center;
		allSphereSource[roiIndex]=vtkSphereSource::New();
		allSpheresData[roiIndex]=vtkPolyData::New();
		allSphereSource[roiIndex]->SetThetaResolution(sphereResolution);
		allSphereSource[roiIndex]->SetPhiResolution(sphereResolution);
		allSphereSource[roiIndex]->SetCenter(roiCenters[roiIndex].x,roiCenters[roiIndex].y,roiCenters[roiIndex].z);
		allSphereSource[roiIndex]->SetRadius(sphereRadius);
		allSpheresData[roiIndex]= allSphereSource[roiIndex]->GetOutput();
		allSpheresData[roiIndex]->Update();
	} //end of for loop:: roiIndex

	vector<Vector3D<float> > allInOnePoints;
	vector<list<int> > allInOneCells;
	int numCells=0;
	int allCellSize=0;


	{
		//put all spheres  together;
		for (int centerIndex = 0 ; centerIndex < numCenters; ++centerIndex) {
			//points;
			for (int pointId = 0 ; pointId < allSpheresData[centerIndex]->GetNumberOfPoints(); ++pointId) {
				Vector3D<float> pointCoord;
				pointCoord.x=allSpheresData[centerIndex]->GetPoint(pointId)[0];
				pointCoord.y=allSpheresData[centerIndex]->GetPoint(pointId)[1];
				pointCoord.z=allSpheresData[centerIndex]->GetPoint(pointId)[2];
				allInOnePoints.push_back(pointCoord);
			} //end of for loop:: pointId

			//cells;
			for (int cellId = 0; cellId < allSpheresData[centerIndex]->GetNumberOfCells(); ++cellId) {
				list<int> newCell;
				vtkCell* currentCell= allSpheresData[centerIndex]->GetCell(cellId);
				for (int cellElem = 0 ; cellElem < currentCell->GetNumberOfPoints(); ++cellElem) {
					int pointidIncurrentSphere= currentCell->GetPointId(cellElem);
					int pointIdInAll=pointidIncurrentSphere+allSpheresData[centerIndex]->GetNumberOfPoints()*centerIndex;
					newCell.push_back(pointIdInAll);
				} //end of for loop:: cellElem
				allInOneCells.push_back(newCell);
				allCellSize+=newCell.size();
				allCellSize+=1;
			} //end of for loop:: cellId

		} //end of for loop:: centerIndex

	}

	// WriteOutAllInVTK(allInOnePoints,allInOneCells,allCellSize,subjectROIname);
	cout<<"write all objects in one vtk file: "<<fileName<<endl;
	fstream vtkStrm;
	OpenWriteStrmAscii(vtkStrm,fileName.c_str());
	vtkStrm<<"# vtk DataFile Version 3.0\n"<<"vtk output\n"<<"ASCII\n"<<"DATASET POLYDATA"<<endl;
	vtkStrm<<"POINTS "<<allInOnePoints.size()<<" FLOAT"<<endl;
	for (int pointId = 0 ; pointId < allInOnePoints.size(); ++pointId) {
		vtkStrm<<allInOnePoints[pointId]<<endl;
	} //end of for loop:: pointId

	vtkStrm<<"POLYGONS "<<allInOneCells.size()<<" "<<allCellSize<<endl;

	for (int cellId = 0; cellId < allInOneCells.size(); ++cellId) {
		vtkStrm<<allInOneCells[cellId].size()<<" ";
		for(list<int>::iterator itCellElems=allInOneCells[cellId].begin(); itCellElems!=allInOneCells[cellId].end(); itCellElems++)
			vtkStrm<<*itCellElems<<" ";
		vtkStrm<<endl;
	} //end of for loop:: cellId

/*
	if(needColor)
	{
	  CColorSchemeRGB myColors; 
	  vtkStrm << "POINT_DATA " << allInOnePoints.size() << endl;
	  vtkStrm << "COLOR_SCALARS " << "PointColor  3 " << endl;
	  for (int idx = 0; idx < allInOnePoints.size() ; ++idx) {
	    int sphereIndex= idx/nSpheres; 
	    sphereIndex=sphereIndex%255; 
	    const RGBTYPE& theColor= myColors.GetColorByIndex(sphereIndex);
	    vtkStrm<<theColor<<endl;
	  }//end of for loop::idx
        }
        
        */
	vtkStrm.close();
}


void VisualizePointsUsingVerticeByVTK(vector< KML::Vector3D< float > >& allInOnePoints, string fileName)
{
  	// WriteOutAllInVTK(allInOnePoints,allInOneCells,allCellSize,subjectROIname);
	cout<<"write all objects in one vtk file: "<<fileName<<endl;
	fstream vtkStrm;
	OpenWriteStrmAscii(vtkStrm,fileName.c_str());
	vtkStrm<<"# vtk DataFile Version 3.0\n"<<"vtk output\n"<<"ASCII\n"<<"DATASET POLYDATA"<<endl;
	vtkStrm<<"POINTS "<<allInOnePoints.size()<<" FLOAT"<<endl;
	for (int pointId = 0 ; pointId < allInOnePoints.size(); ++pointId) {
		vtkStrm<<allInOnePoints[pointId]<<endl;
	} //end of for loop:: pointId

	vtkStrm<<"VERTICES "<<allInOnePoints.size()<<" "<<allInOnePoints.size()*2<<endl;

	for (int cellId = 0; cellId < allInOnePoints.size(); ++cellId) {
		vtkStrm<<1<<" "<<cellId<<endl;
	} //end of for loop:: cellId
	

	vtkStrm.close();
}


void VisualizePointsUsingVerticeByVTKWithColors(vector< KML::Vector3D< float > >& allInOnePoints, vector< KML::Vector3D< float > >& allPointColors, string fileName)
{
        // WriteOutAllInVTK(allInOnePoints,allInOneCells,allCellSize,subjectROIname);
        cout<<"write all objects in one vtk file: "<<fileName<<endl;
        fstream vtkStrm;
        OpenWriteStrmAscii(vtkStrm,fileName.c_str());
        vtkStrm<<"# vtk DataFile Version 3.0\n"<<"vtk output\n"<<"ASCII\n"<<"DATASET POLYDATA"<<endl;
        vtkStrm<<"POINTS "<<allInOnePoints.size()<<" FLOAT"<<endl;
        for (int pointId = 0 ; pointId < allInOnePoints.size(); ++pointId) {
                vtkStrm<<allInOnePoints[pointId]<<endl;
        } //end of for loop:: pointId

        vtkStrm<<"VERTICES "<<allInOnePoints.size()<<" "<<allInOnePoints.size()*2<<endl;

        for (int cellId = 0; cellId < allInOnePoints.size(); ++cellId) {
                vtkStrm<<1<<" "<<cellId<<endl;
        } //end of for loop:: cellId
        
        vtkStrm<<"POINT_DATA "<<allInOnePoints.size()<<endl;
        vtkStrm << "COLOR_SCALARS " << "PointColor  3 " << endl;
        
        for(int idx = 0; idx < allInOnePoints.size(); ++idx )
        {  
          vtkStrm<< allPointColors[idx]<<endl;
        } // end for::idx;
        
        vtkStrm.close();
}


















Fvector3d *** Fvector3dalloc3d(int i_size,int j_size,int k_size)
{
  Fvector3d ***array;
  int i,j,k;

  array=(Fvector3d ***) calloc(k_size,sizeof(Fvector3d **));

  for(k=0;k<k_size;k++)
    array[k]=(Fvector3d **) calloc(i_size,sizeof(Fvector3d *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(Fvector3d *) calloc(j_size,sizeof(Fvector3d ));

  return(array);
}

void Fvector3dfree3d(Fvector3d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

void nrerrorSHEN(const char *error_text)
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *dvectorSHEN(int nl, int nh)
{
	double *v;

	v = (double *) calloc((unsigned) (nh-nl+1), sizeof(double));
	if (!v) nrerrorSHEN("allocation failure in dvectorSHEN()");
	return v-nl;
}

void free_dvectorSHEN(double *v, int nl, int nh)
{
	free((char*) (v+nl));
}


float *vectorSHEN(int nl, int nh)
{
	float *v;

	v = (float *) calloc((unsigned) (nh-nl+1), sizeof(float));
	if (!v) nrerrorSHEN("allocation failure in dvectorSHEN()");
	return v-nl;
}

void free_vectorSHEN(float *v, int nl, int nh)
{
	free((char*) (v+nl));
}

float **matrixSHEN(int nrl,int nrh,int ncl,int nch)
{
	int i;
	float **m;

	m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
	if (!m) nrerrorSHEN("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
		if (!m[i]) nrerrorSHEN("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_matrixSHEN(float **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


double **dmatrixSHEN(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m = (double **) calloc((unsigned) (nrh-nrl+1), sizeof(double*));
	if (!m) nrerrorSHEN("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i] = (double *) calloc((unsigned) (nch-ncl+1), sizeof(double));
		if (!m[i]) nrerrorSHEN("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_dmatrixSHEN(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


/*double log2(double a)
{
	return log10(a)/log10(2.0);
}*/

void sort(double *Y, int *I, double *A, int length)
{
	int i, j;
	double max, *tmp;

	tmp = (double *) calloc(length, sizeof(double));

	for (i=0;i<length;i++)
		tmp[i] = A[i];

	max = tmp[0];
	for (i=1;i<length;i++) {
		if (tmp[i] > max)
			max = tmp[i];
	}

	max = fabs(10*max);

	for (i=0;i<length;i++) {
		Y[i] = tmp[0];
		I[i] = 0;
		for (j=1;j<length;j++) {
			if (tmp[j] < Y[i]) {
				Y[i] = tmp[j];
				I[i] = j;
			}
		}

		tmp[I[i]] = max;
	}

	free(tmp);
}

void minimun(double *Y, int *I, double *A, int length)
{
	int i, index;
	double min;

	min = A[0];
	index = 0;
	for (i=1;i<length;i++)
		if (A[i] < min) {
			min = A[i];
			index = i;
		}

	*Y = min;
	*I = index;
}

void Mat_Abs(MatrixShen *A)
{
	int h, w;

	for (h=0;h<A->height;h++)
		for (w=0;w<A->width;w++) {
			if (A->data[h][w] < 0)
			    A->data[h][w] = -1.0*(A->data[h][w]);
		}
}

void Mat_Mean(double *mean, MatrixShen *A)
{
	int h, w;
	double tmp;

	tmp = 0.0;
	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			tmp += A->data[h][w];
		}
	}

	*mean = tmp/(double) (A->height*A->width);
}

void Mat_Variance(double *variance, MatrixShen *A)
{
	int h, w;
	double mean, tmp;

	Mat_Mean(&mean, A) ;

	tmp = 0.0;
	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			tmp += pow(A->data[h][w]-mean,2.0) ;
		}
	}

	*variance = sqrt(tmp/(double)(A->height*A->width));
}


void Mat_Vector(MatrixShen *A, float *a)
{
	int h, w;

	for (h=0;h<A->height;h++)
		for (w=0;w<A->width;w++)
			a[h*A->width+w] = (float) A->data[h][w];
}

void Mat_Shift(MatrixShen *A, MatrixShen *B, int side)
{
	int h, w;

	for (h=side;h<B->height;h++)
		for (w=side;w<B->width;w++)
			A->data[h-side][w-side] = B->data[h][w];

	for (h=side;h<B->height;h++)
		for (w=0;w<side;w++)
			A->data[h-side][B->width-side+w] = B->data[h][w];

	for (h=0;h<side;h++)
		for (w=side;w<B->width;w++)
			A->data[B->height-side+h][w-side] = B->data[h][w];

	for (h=0;h<side;h++)
		for (w=0;w<side;w++)
			A->data[B->height-side+h][B->width-side+w] = B->data[h][w];
}

void Mat_Zeros(MatrixShen *A)
{
	int h, w;

	for (h=0;h<A->height;h++)
		for (w=0;w<A->width;w++)
			A->data[h][w] = 0;
}

void Mat_Zeros_uc(uc_Matrix *A)
{
	int h, w;

	for (h=0;h<A->height;h++)
		for (w=0;w<A->width;w++)
			A->data[h][w] = 0;
}

void Mat_Zeros_i(i_Matrix *A)
{
	int h, w;

	for (h=0;h<A->height;h++)
		for (w=0;w<A->width;w++)
			A->data[h][w] = 0;
}


void CreateMatrix(MatrixShen **M, int hei, int wid)
{
	int h;

	MatrixShen *tmp;

	tmp = (MatrixShen *) calloc(1, sizeof(MatrixShen));
	tmp->data = (double **) calloc(hei, sizeof(double *));

	if (!(tmp->data)) {
		nrerrorSHEN("allocation failure in CreateMatrix()");
		exit(1);
	}

	for (h=0; h<hei; h++) {
		tmp->data[h] = (double *) calloc(wid, sizeof(double));
		if (!(tmp->data[h])) {
			nrerrorSHEN("allocation failure in CreateMatrix()");
			exit(1);
		}
	}

	tmp->height = hei;
	tmp->width = wid;

	*M = tmp;
}


void FreeMatrix(MatrixShen *M)
{
	int h, hei = M->height;

	for (h=0; h<hei; h++) {
	     free(M->data[h]);
	}
	free(M->data);
	free(M);
}

void Create_i_Matrix(i_Matrix **M, int hei, int wid)
{
	int h;

	i_Matrix *tmp;

	tmp = (i_Matrix *) calloc(1, sizeof(i_Matrix));
	tmp->data = (int **) calloc(hei, sizeof(int *));
	if (!(tmp->data)) {
		nrerrorSHEN("allocation failure in Create_i_Matrix()");
		exit(1);
	}

	for (h=0; h<hei; h++) {
		tmp->data[h] = (int *) calloc(wid, sizeof(int));
		if (!(tmp->data[h])) {
			nrerrorSHEN("allocation failure in Create_i_Matrix()");
			exit(1);
		}
	}

	tmp->height = hei;
	tmp->width = wid;

	*M = tmp;
}

void Free_i_Matrix(i_Matrix *M)
{
	int h;

	for (h=0; h<M->height; h++)
	     free(M->data[h]);
	free(M->data);
	free(M);
}

void Create_uc_Matrix(uc_Matrix **M, int hei, int wid)
{
	int h;

	uc_Matrix *tmp;

	tmp = (uc_Matrix *) calloc(1, sizeof(uc_Matrix));
	tmp->data = (unsigned char **) calloc(hei, sizeof(unsigned char *));
	if (!(tmp->data)) {
		nrerrorSHEN("allocation failure in Create_uc_Matrix()");
		exit(1);
	}

	for (h=0; h<hei; h++) {
		tmp->data[h] = (unsigned char *) calloc(wid, sizeof(unsigned char));
		if (!(tmp->data[h])) {
			nrerrorSHEN("allocation failure in Create_uc_Matrix()");
			exit(1);
		}
	}

	tmp->height = hei;
	tmp->width = wid;

	*M = tmp;
}

void Free_uc_Matrix(uc_Matrix *M)
{
	int h;

	for (h=0; h<M->height; h++)
	     free(M->data[h]);
	free(M->data);
	free(M);
}


void Mat_FFT2(MatrixShen *Output_real, MatrixShen *Output_imag, MatrixShen *Input_real, MatrixShen *Input_imag)
{
	int xs, ys, i, j;
	double **R, **I, **Fr, **Fi;

	xs = Input_real->height;
	ys = Input_real->width;

	R  = dmatrixSHEN(1,xs,1,ys);
	I  = dmatrixSHEN(1,xs,1,ys);
	Fr = dmatrixSHEN(1,xs,1,ys);
	Fi = dmatrixSHEN(1,xs,1,ys);

	for (i=1;i<=Input_real->height;i++)
	    for (j=1;j<=Input_real->width;j++) {
		R[i][j] = Input_real->data[i-1][j-1];
		I[i][j] = Input_imag->data[i-1][j-1];
	    }

	four2(Fr, Fi, R, I, xs, ys, 1);         /* 2-D FFT */

	for (i=1;i<=Input_real->height;i++)
	    for (j=1;j<=Input_real->width;j++) {
		Output_real->data[i-1][j-1] = Fr[i][j];
		Output_imag->data[i-1][j-1] = Fi[i][j];
	    }

	free_dmatrixSHEN(R,1,xs,1,ys);
	free_dmatrixSHEN(I,1,xs,1,ys);
	free_dmatrixSHEN(Fr,1,xs,1,ys);
	free_dmatrixSHEN(Fi,1,xs,1,ys);
}

void Mat_IFFT2(MatrixShen *Output_real, MatrixShen *Output_imag, MatrixShen *Input_real, MatrixShen *Input_imag)
{
	int xs, ys, i, j;
	double **R, **I, **Fr, **Fi, NN;

	xs = Input_real->height;
	ys = Input_real->width;

	R  = dmatrixSHEN(1,xs,1,ys);
	I  = dmatrixSHEN(1,xs,1,ys);
	Fr = dmatrixSHEN(1,xs,1,ys);
	Fi = dmatrixSHEN(1,xs,1,ys);

	for (i=1;i<=Input_real->height;i++)
	    for (j=1;j<=Input_real->width;j++) {
		R[i][j] = Input_real->data[i-1][j-1];
		I[i][j] = Input_imag->data[i-1][j-1];
	    }

	four2(Fr, Fi, R, I, xs, ys, -1);         /* 2-D IFFT */

	NN = (double) (xs*ys);

	for (i=1;i<=Input_real->height;i++)
	    for (j=1;j<=Input_real->width;j++) {
		Output_real->data[i-1][j-1] = Fr[i][j]/NN;
		Output_imag->data[i-1][j-1] = Fi[i][j]/NN;
	    }

	free_dmatrixSHEN(R,1,xs,1,ys);
	free_dmatrixSHEN(I,1,xs,1,ys);
	free_dmatrixSHEN(Fr,1,xs,1,ys);
	free_dmatrixSHEN(Fi,1,xs,1,ys);
}

void four2(double **fftr, double **ffti, double **rdata, double **idata, int rs, int cs, int isign)
/************************************************************

   2-D fourier transform of data with real part stored in
   "rdata" and imaginary part in "idata" with size "rs" x
   "cs". The result is in "fftr" and "ffti". The isign is
   "isign" =  1 forward, and "isign" = -1 inverse

*************************************************************/
{
	double **T, *tmp1, *tmp2;
	int i, j;

	tmp1 = dvectorSHEN(1,2*cs);
	tmp2 = dvectorSHEN(1,2*rs);
	T = dmatrixSHEN(1,2*rs,1,cs);

	for (i=1;i<=rs;i++) {
	    for (j=1;j<=cs;j++) {
		tmp1[j*2-1] = rdata[i][j];
		tmp1[j*2] = idata[i][j];
	    }
	    four1(tmp1, cs, isign);
	    for (j=1;j<=cs;j++) {
		T[i*2-1][j] = tmp1[j*2-1];
		T[i*2][j] = tmp1[j*2];
	    }
	}

	for (i=1;i<=cs;i++) {
	    for (j=1;j<=rs;j++) {
		tmp2[j*2-1] = T[j*2-1][i];
		tmp2[j*2] = T[j*2][i];
	    }
	    four1(tmp2,rs,isign);
	    for (j=1;j<=rs;j++) {
		fftr[j][i] = tmp2[j*2-1];
		ffti[j][i] = tmp2[j*2];
	    }
	}
	free_dvectorSHEN(tmp1, 1, 2*cs);
	free_dvectorSHEN(tmp2, 1, 2*rs);
	free_dmatrixSHEN(T, 1, 2*rs, 1, cs);
}

void four1(double *data, int nn, int isign)
{
	int n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	n = nn << 1;
	j = 1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (n > mmax) {
		istep = 2*mmax;
		theta = 6.28318530717959/(isign*mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j = i+mmax;
				tempr = wr*data[j]-wi*data[j+1];
				tempi = wr*data[j+1]+wi*data[j];
				data[j] = data[i]-tempr;
				data[j+1] = data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr = (wtemp=wr)*wpr-wi*wpi+wr;
			wi = wi*wpr+wtemp*wpi+wi;
		}
		mmax = istep;
	}
}

void Mat_Copy(MatrixShen *A, MatrixShen *B, int h_target, int w_target, int h_begin, int w_begin, int h_end, int w_end)
{
	int i, j, h, w, h_done, w_done;

	if ((h_target >= 0)&&(h_target < A->height)&&(w_target >= 0)&&(w_target < A->width)) {
		if ((h_begin >= 0)&&(h_begin < B->height)&&(w_begin >= 0)&&(w_begin < B->width)) {
			h = h_end-h_begin+1;
			w = w_end-w_begin+1;
			if ((h >= 1)&&(w >= 1)) {
				h_done = h_target+h-1;
				w_done = w_target+w-1;
				if ((h_done < A->height)&&(w_done < A->width)) {
					for (i=0;i<h;i++) {
						for (j=0;j<w;j++) {
							A->data[i+h_target][j+w_target] = B->data[i+h_begin][j+w_begin];
						}
					}
				}
			}
		}
	}
	else {
		printf("matrix dimension error!\n");
		exit(1);
	}
}

void Mat_uc_Copy(uc_Matrix *A, uc_Matrix *B, int h_target, int w_target, int h_begin, int w_begin, int h_end, int w_end)
{
	int i, j, h, w, h_done, w_done;

	if ((h_target >= 0)&&(h_target < A->height)&&(w_target >= 0)&&(w_target < A->width)) {
		if ((h_begin >= 0)&&(h_begin < B->height)&&(w_begin >= 0)&&(w_begin < B->width)) {
			h = h_end-h_begin+1;
			w = w_end-w_begin+1;
			if ((h >= 1)&&(w >= 1)) {
				h_done = h_target+h-1;
				w_done = w_target+w-1;
				if ((h_done < A->height)&&(w_done < A->width)) {
					for (i=0;i<h;i++) {
						for (j=0;j<w;j++) {
							A->data[i+h_target][j+w_target] = B->data[i+h_begin][j+w_begin];
						}
					}
				}
			}
		}
	}
	else {
		printf("matrix dimension error!\n");
		exit(1);
	}
}

void Mat_i_Copy(i_Matrix *A, i_Matrix *B, int h_target, int w_target, int h_begin, int w_begin, int h_end, int w_end)
{
	int i, j, h, w, h_done, w_done;

	if ((h_target >= 0)&&(h_target < A->height)&&(w_target >= 0)&&(w_target < A->width)) {
		if ((h_begin >= 0)&&(h_begin < B->height)&&(w_begin >= 0)&&(w_begin < B->width)) {
			h = h_end-h_begin+1;
			w = w_end-w_begin+1;
			if ((h >= 1)&&(w >= 1)) {
				h_done = h_target+h-1;
				w_done = w_target+w-1;
				if ((h_done < A->height)&&(w_done < A->width)) {
					for (i=0;i<h;i++) {
						for (j=0;j<w;j++) {
							A->data[i+h_target][j+w_target] = B->data[i+h_begin][j+w_begin];
						}
					}
				}
			}
		}
	}
	else {
		printf("matrix dimension error!\n");
		exit(1);
	}
}

void Mat_Product(MatrixShen *A, MatrixShen *B, MatrixShen *C)
{
	int h, w;


	if(A->height!=B->height || A->height!=C->height ||
	   A->width!=B->width || A->width!=C->width )
	  nrerrorSHEN("Mat_Substract fail!");


	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = B->data[h][w]*C->data[h][w];
		}
	}
}

void Mat_Sum(MatrixShen *A, MatrixShen *B, MatrixShen *C)
{
	int h, w;


	if(A->height!=B->height || A->height!=C->height ||
	   A->width!=B->width || A->width!=C->width )
	  nrerrorSHEN("Mat_Substract fail!");


	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = B->data[h][w]+C->data[h][w];
		}
	}
}

void Mat_Substract(MatrixShen *A, MatrixShen *B, MatrixShen *C)
{
	int h, w;

	if(A->height!=B->height || A->height!=C->height ||
	   A->width!=B->width || A->width!=C->width )
	  nrerrorSHEN("Mat_Substract fail!");


	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = B->data[h][w]-C->data[h][w];
		}
	}
}

void Mat_Fliplr(MatrixShen *A)
{
	MatrixShen *tmp;
	int h, w;

	CreateMatrix(&tmp, A->height, A->width);

	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			tmp->data[h][w] = A->data[h][(A->width)-w-1];
		}
	}

	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = tmp->data[h][w];
		}
	}

	FreeMatrix(tmp);
}

void Mat_Flipud(MatrixShen *A)
{
	MatrixShen *tmp;
	int h, w;

	CreateMatrix(&tmp, A->height, A->width);

	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			tmp->data[h][w] = A->data[(A->height)-h-1][w];
		}
	}

	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = tmp->data[h][w];
		}
	}

	FreeMatrix(tmp);
}


void Mat_uc_Fliplr(uc_Matrix *A)
{
	uc_Matrix *tmp;
	int h, w;

	Create_uc_Matrix(&tmp, A->height, A->width);

	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			tmp->data[h][w] = A->data[h][(A->width)-w-1];
		}
	}

	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = tmp->data[h][w];
		}
	}

	Free_uc_Matrix(tmp);
}

void Mat_uc_Flipud(uc_Matrix *A)
{
	uc_Matrix *tmp;
	int h, w;

	Create_uc_Matrix(&tmp, A->height, A->width);

	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			tmp->data[h][w] = A->data[(A->height)-h-1][w];
		}
	}

	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = tmp->data[h][w];
		}
	}

	Free_uc_Matrix(tmp);
}


/*add by SHEN when in JHU*/
int *ivectorSHEN(long nl, long nh)
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerrorSHEN("allocation failure in ivector()");
	return v-nl+NR_END;
}
void free_ivectorSHEN(int *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}


void gaussj(float **a, int n, float **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	float big,dum,pivinv,temp;

	indxc=ivectorSHEN(1,n);
	indxr=ivectorSHEN(1,n);
	ipiv=ivectorSHEN(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) nrerrorSHEN("gaussj: Singular Matrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerrorSHEN("gaussj: Singular Matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivectorSHEN(ipiv,1,n);
	free_ivectorSHEN(indxr,1,n);
	free_ivectorSHEN(indxc,1,n);
}

float pythag(float a, float b)
{
	float absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void svdcmpSHEN(float **a, int m, int n, float w[], float **v)
{
	float pythag(float a, float b);
	int flag,i,its,j,jj,k,l,nm;
	float anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=vectorSHEN(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((float)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((float)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((float)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) nrerrorSHEN("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vectorSHEN(rv1,1,n);
}

void Mat_A_equal_BxC(MatrixShen *A, MatrixShen *B, MatrixShen *C)
{
	int h, w, k;
	double sum ;

	if(B->width!=C->height || A->height!=B->height || A->width!=C->width )
	  {
	    printf("Matrix operation error!\n");
	    exit(1);
	  }


	for (h=0; h<A->height; h++)
	  {
	    for (w=0; w<A->width; w++)
	      {
		sum = 0 ;
		for(k=0; k<B->width; k++)
		  sum += (B->data[h][k] * C->data[k][w]) ;

		A->data[h][w] = sum ;
	      }
	  }

}

void Mat_Print(MatrixShen *A)
{
  int k,l ;

  printf("\n");

  for(k=0; k<A->height; k++)
    {
      for (l=0; l<A->width; l++)
	printf("%12.3 f", A->data[k][l]) ;
      printf("\n");
    }
}


void Mat_Calculate_EigenVectors_EigenValues(MatrixShen *C, float *EigenValue, MatrixShen *EigenVector, int PRNorNOT)
{
  int   j, k, l,   np, mp,  position ;
  float tmp, max ;

  /* for getting eigenVectors (u) and eigenValues (w)*/
  float *w,**u,**v;


  /* Calculate eigenvectors and eigenvalues */
   mp = C->height ;  /* size */
   np = C->width ;
   /* apply memory */
   /* Notice here the beginning index for u (v or w) is 1, not 0.  */
   w = vectorSHEN(1,np);
   u = matrixSHEN(1,mp,1,np);
   v = matrixSHEN(1,np,1,np);

   /* set matrices: u <- C */
   for(k=1; k<=mp; k++)
     for(l=1; l<=np; l++)
       u[k][l] = C->data[k-1][l-1] ;


  if(PRNorNOT==TRUE)
  {
    printf("\nCheck product against original matrix:\n");
    printf("Original matrix:\n");
    for (k=1;k<=mp;k++)
      {
	for (l=1;l<=np;l++)
	  printf("%12.3f ",u[k][l]);
	printf("\n");
      }
  }

  /* perform decomposition */
  svdcmpSHEN(u,mp,np,w,v);


  if(PRNorNOT==TRUE)
  {
    printf("Product u*w*(v-transpose):\n");
    for (k=1;k<=mp;k++)
      {
	for (l=1;l<=np;l++)
	  {
	    tmp=0.0;
	    for (j=1;j<=np;j++)
	      tmp += u[k][j]*w[j]*v[l][j];
	    printf("%12.3 f",tmp) ;
	  }
	printf("\n");
      }

    /* write results */
    printf("Decomposition matrices:\n");
    printf("Matrix u\n");
    for (k=1; k<=mp; k++)
      {
	for (l=1; l<=np; l++)
	  printf("%12.3f ", u[k][l]);
	printf("\n");
      }

    printf("Diagonal of matrix w\n");
    for (k=1;k<=np;k++)
      printf("%12.3f ", w[k]);

    printf("\nMatrix v-transpose\n");
    for (k=1;k<=np;k++)
      {
	for (l=1;l<=np;l++)
	  printf("%12.3f ", v[l][k]);
	printf("\n");
      }
  }


  /* record for returning */
   for(k=1; k<=np; k++)
     {
       max = -100000.0 ;
      for(j=1; j<=np; j++)
       {
	 if( w[j]>max )
	   {
	     max = w[j] ;
	     position = j ;
	   }
       }
      /*eigenvalue*/
      EigenValue[k-1] = w[position] ;
      /*eigenvector*/
      for(l=1;l<=np;l++)
	  EigenVector->data[l-1][k-1] = u[l][position] ;
      w[position] = -1000000.0 ;
     }


   free_vectorSHEN(w,1,np);
   free_matrixSHEN(u,1,mp,1,np);
   free_matrixSHEN(v,1,np,1,np);
}


void Mat_Inverse(MatrixShen *A, MatrixShen *B) /*A in, B out*/
{
  int    i, j, m, n1, n2 ;
  float  **a,**b;
  MatrixShen *At, *AtA, *AtAInv ;

  n1 = A->height ;
  n2 = A->width ;
  m = 1 ;

  if(n1==n2)
    {
      if( B->height != n1 || B->width != n2 )
	nrerrorSHEN("marix size not matching") ;

      /* apply memory */
      /* Notice here the beginning index for u (v or w) is 1, not 0.  */
      a = matrixSHEN(1,n1,1,n2);
      b = matrixSHEN(1,n1,1,m);

      /* set matrices: a <- A */
      for(i=1; i<=n1; i++)
	for(j=1; j<=n2; j++)
	  a[i][j] = A->data[i-1][j-1] ;

      /*cheat b*/
     for(i=1; i<=n1; i++)
       for(j=1; j<=m; j++)
	  b[i][j] = A->data[i-1][j-1] ;

     gaussj(a, n1, b, m) ; /* get inverse of a */

      /* set matrices: A <- a */
      for(i=1; i<=n1; i++)
	for(j=1; j<=n2; j++)
	  B->data[i-1][j-1] = a[i][j] ;

     free_matrixSHEN(a,1,n1,1,n2);
     free_matrixSHEN(b,1,n1,1,m);
    }
  else /* n1 != n2*/
    {
      if( B->height != n2 || B->width != n1 )
	nrerrorSHEN("marix size not matching") ;

      /*printf("n1=%d n2=%d\n", n1, n2) ;*/

      CreateMatrix(&At,  n2, n1);
      CreateMatrix(&AtA, n2, n2);
      CreateMatrix(&AtAInv, n2, n2); /* used as the inverse matrix of C*/

      /* B=A(t)*/
      for(i=0; i<n2; i++)
	for(j=0; j<n1; j++)
	  At->data[i][j] = A->data[j][i] ;

printf("At\n") ; Mat_Print(At) ;
printf("A\n") ; Mat_Print(A) ;

      Mat_A_equal_BxC(AtA, At,  A) ; /*A=BxC*/
printf("AtA\n") ; Mat_Print(AtA) ;

      Mat_Inverse(AtA, AtAInv) ;     /*A in, B out*/
printf("AtAI\n") ; Mat_Print(AtAInv) ;

      Mat_A_equal_BxC(B, AtAInv,  At) ; /*A=BxC*/


      FreeMatrix(At) ;
      FreeMatrix(AtA) ;
      FreeMatrix(AtAInv) ;
    }
}


void Mat_times_Vector(float *Vout, MatrixShen *A, float *Vin)
{
  int    i, j, n1, n2 ;
  double sum ;

  n1 = A->height ;
  n2 = A->width ;

  /* set matrices: a <- A */
  for(i=0; i<n1; i++)
    {
      sum = 0 ;
      for(j=0; j<n2; j++)
	sum += A->data[i][j]*Vin[j] ;
      Vout[i] = sum ;
    }
}


void vector_Print(float *v, int size)
{
  int k ;

  printf("\n");

  for(k=0; k<size; k++)
    printf("%12.3f ", v[k]) ;
}


void Mat_EqualCopy(MatrixShen *A, MatrixShen *B)
{
	int i, j, h, w, h_done, w_done;

	if ( A->height==B->height && A->width==B->width )
	  {
	    h = A->height ;
	    w = A->width  ;

	    for (i=0; i<h; i++)
	      for (j=0; j<w; j++)
		A->data[i][j] = B->data[i][j];
	  }
	else
	  {
	    printf("Matrix copy: matrix dimension error!\n");
	    exit(1);
	  }
}



#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


float gasdev(long *idum)
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


/* June 2001  */
ImgAttribute ****ImgAttributealloc4d(int i_size,int j_size,int k_size, int t_size)
{
  ImgAttribute ****array;
  int i,j,k, t;

  array=(ImgAttribute ****) calloc(t_size,sizeof(ImgAttribute ***));

  for(t=0;t<t_size;t++)
    array[t]=(ImgAttribute ***) calloc(k_size,sizeof(ImgAttribute **));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      array[t][k]=(ImgAttribute **) calloc(i_size,sizeof(ImgAttribute *));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	array[t][k][i]=(ImgAttribute *) calloc(j_size,sizeof(ImgAttribute ));

  return(array);
}

ImgAttribute ***ImgAttributealloc3d(int i_size,int j_size,int k_size)
{
  ImgAttribute ***array;
  int i,j,k;

  array=(ImgAttribute ***) calloc(k_size,sizeof(ImgAttribute **));

  for(k=0;k<k_size;k++)
    array[k]=(ImgAttribute **) calloc(i_size,sizeof(ImgAttribute *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(ImgAttribute *) calloc(j_size,sizeof(ImgAttribute ));

  return(array);
}

ImgAttribute *ImgAttributealloc1d(int k_size)
{
  ImgAttribute *array;
  int i,j,k;

  array=(ImgAttribute *) calloc(k_size,sizeof(ImgAttribute));

  return(array);
}

/*free*/
void ImgAttributefree4d(ImgAttribute ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}

void ImgAttributefree3d(ImgAttribute ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}


/* Nov 2001, for mapping volumetric head images */
HeadImgAttribute ***HeadImgAttributealloc3d(int i_size,int j_size,int k_size)
{
  HeadImgAttribute ***array;
  int i,j,k;

  array=(HeadImgAttribute ***) calloc(k_size,sizeof(HeadImgAttribute **));

  for(k=0;k<k_size;k++)
    array[k]=(HeadImgAttribute **) calloc(i_size,sizeof(HeadImgAttribute *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(HeadImgAttribute *) calloc(j_size,sizeof(HeadImgAttribute ));

  return(array);
}

HeadImgAttribute *HeadImgAttributealloc1d(int k_size)
{
  HeadImgAttribute *array;
  int i,j,k;

  array=(HeadImgAttribute *) calloc(k_size,sizeof(HeadImgAttribute));

  return(array);
}

/*free*/
void HeadImgAttributefree3d(HeadImgAttribute ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}



/* Feb 2002, for warping DTI */
DTIattribute ***DTIattributealloc3d(int i_size,int j_size,int k_size)
{
  DTIattribute ***array;
  int i,j,k;

  array=(DTIattribute ***) calloc(k_size,sizeof(DTIattribute **));

  for(k=0;k<k_size;k++)
    array[k]=(DTIattribute **) calloc(i_size,sizeof(DTIattribute *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(DTIattribute *) calloc(j_size,sizeof(DTIattribute ));

  return(array);
}

DTIattribute *DTIattributealloc1d(int k_size)
{
  DTIattribute *array;
  int i,j,k;

  array=(DTIattribute *) calloc(k_size,sizeof(DTIattribute));

  return(array);
}


/*free*/
void DTIattributefree3d(DTIattribute ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

void ConvertMatrix2Shen(Matrix& newmat, MatrixShen& shenmat){
	if(newmat.Nrows()!= shenmat.height && newmat.Ncols()!=shenmat.width){
		cout<<"dimension error"<<endl;
		exit(-1);
	}

	for (int row = 0; row < newmat.Nrows(); ++row) {
		for (int col = 0; col < newmat.Ncols(); ++col) {
			shenmat.data[row][col]=newmat(row+1,col+1);
		} //end of for loop:: col
	} //end of for loop:: row
}

void OutputMatShen  (  MatrixShen& shenmat){
	for (int row = 0; row < shenmat.height; ++row) {
		for (int col = 0; col < shenmat.width; ++col) {
			cout<<shenmat.data[row][col]<<" ";
		} //end of for loop:: col
		cout<<endl;
	} //end of for loop:: row
}


void ConvertMatrix2Shen(Matrix& newmat, MatrixShen* shenmat){
	ConvertMatrix2Shen(newmat,*shenmat);
}

void OutputMatShen(  MatrixShen* shenmat){
	OutputMatShen( *shenmat);
}

float *Falloc1d(int i_size)
{
  float *array;
  int i;

  array=(float *) calloc(i_size,sizeof(float ));

  return(array);
}


float **Falloc2d(int i_size,int j_size)
{
  float **array;
  int i,j;

  array=(float **) calloc(i_size,sizeof(float *));

  for(i=0;i<i_size;i++)
    array[i]=(float *) calloc(j_size,sizeof(float ));

  return(array);
}

void Ffree2d(float **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}

unsigned char ***UCalloc3d(int i_size,int j_size,int k_size)
{
  unsigned char ***array;
  int i,j,k;

  array=(unsigned char ***) calloc(k_size,sizeof(unsigned char **));

  for(k=0;k<k_size;k++)
    array[k]=(unsigned char **) calloc(i_size,sizeof(unsigned char *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(unsigned char *) calloc(j_size,sizeof(unsigned char ));
	
  return(array);
}

void UCfree3d(unsigned char ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

int ***Ialloc3d(int i_size,int j_size,int k_size)
{
  int ***array;
  int i,j,k;

  array=(int ***) calloc(k_size,sizeof(int **));

  for(k=0;k<k_size;k++)
    array[k]=(int **) calloc(i_size,sizeof(int *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(int *) calloc(j_size,sizeof(int ));
	
  return(array);
}

void Ifree3d(int ***array,int k_size,int i_size)
{
  int k,i;


  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);
	
  free(array);
}

int ****Ialloc4d(int i_size,int j_size,int k_size, int t_size)
{
  int ****array;
  int i,j,k, t;

  array=(int ****) calloc(t_size,sizeof(int ***));

  for(t=0;t<t_size;t++)
    array[t]=(int ***) calloc(k_size,sizeof(int **));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      array[t][k]=(int **) calloc(i_size,sizeof(int *));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	array[t][k][i]=(int *) calloc(j_size,sizeof(int ));
	
  return(array);
}

void Ifree4d(int ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}

float ***Falloc3d(int i_size,int j_size,int k_size)
{
  float ***array;
  int i,j,k;

  array=(float ***) calloc(k_size,sizeof(float **));

  for(k=0;k<k_size;k++)
    array[k]=(float **) calloc(i_size,sizeof(float *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(float *) calloc(j_size,sizeof(float));
	
  return(array);
}

void Ffree3d(float ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

float ****Falloc4d(int i_size,int j_size,int k_size, int t_size)
{
	float ****array;
	int i,k,t;

	array=(float ****) calloc(t_size,sizeof(float ***));
	
	for(t=0;t<t_size;t++)
		array[t]=(float ***) calloc(k_size,sizeof(float **));
	
	for(t=0;t<t_size;t++)
		for(k=0;k<k_size;k++)
			array[t][k]=(float **) calloc(i_size,sizeof(float *));
	
	for(t=0;t<t_size;t++)
		for(k=0;k<k_size;k++)
			for(i=0;i<i_size;i++)
				array[t][k][i]=(float *) calloc(j_size,sizeof(float));
	
	return(array);
}

void Ffree4d(float ****array, int i_size, int j_size, int k_size, int t_size)
{
	int t,k,i;
	
	for(t=0;t<t_size;t++)
		for(k=0;k<k_size;k++)
			for(i=0;i<i_size;i++)
				free(array[t][k][i]);
	
	for(t=0;t<t_size;t++)
		for(k=0;k<k_size;k++)
			free(array[t][k]);
	
	for(t=0;t<t_size;t++)
		free(array[t]);

  free(array);
}

unsigned short ***Salloc3d(int i_size,int j_size,int k_size)
{
  unsigned short ***array;
  int i,j,k;

  array=(unsigned short ***) calloc(k_size,sizeof(unsigned short **));

  for(k=0;k<k_size;k++)
    array[k]=(unsigned short **) calloc(i_size,sizeof(unsigned short *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(unsigned short *) calloc(j_size,sizeof(unsigned short ));
	
  return(array);
}

void Sfree3d(unsigned short ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}
int *Ialloc1d(int i_size)
{
  int *array;
  int i;

  array=(int *) calloc(i_size,sizeof(int ));

  return(array);
}


void sort(float *Y, int *I, float *A, int length)
{
	int i, j;
	float max, *tmp;

	tmp = (float *) calloc(length, sizeof(float));

	for (i=0;i<length;i++) 
		tmp[i] = A[i];

	max = tmp[0];
	for (i=1;i<length;i++) {
		if (tmp[i] > max) 
			max = tmp[i];
	}

	max = fabs(10*max);

	for (i=0;i<length;i++) {
		Y[i] = tmp[0];
		I[i] = 0;
		for (j=1;j<length;j++) {
			if (tmp[j] < Y[i]) {
				Y[i] = tmp[j];
				I[i] = j;
			}
		}

		tmp[I[i]] = max;
	}

	free(tmp);
}


void svdcmp(float **a, int m, int n, float w[], float **v)
{
	float pythag(float a, float b);
	int flag,i,its,j,jj,k,l,nm;
	float anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=vectorSHEN(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((float)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((float)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((float)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 100) nrerrorSHEN("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vectorSHEN(rv1,1,n);
}




}
