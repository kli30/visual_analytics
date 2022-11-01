/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#include "fibers.h"
#include "kaimingCommon.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkCell.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkIdList.h"
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataNormals.h>
#include "vtkPolyDataWriter.h"
#include "vtkCurvatures.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <cassert>
#include <sstream>
#include <map>
#include "jointModel.h"


using namespace std;
using namespace KML;

namespace KML
{
	class CFibersPrivacy
	{
		friend class CFibers;
	private:
		CFibersPrivacy(void):numPoints(0),numFibers(0), fileName(""),isASCII(true),resolution(0.f) {};
		~CFibersPrivacy(){};
		vector<PointCoordType> allPoints;
		vector<vector<int> > allFibers; 
		vector<PointCoordType> allPointsColorsAscii;
		vector<Vector3D<unsigned char> > allPointsColorsBinary;
		vector<PointCoordType> allLineColor;
		vector<int> fiberLabel;
		int numPoints;
		int numFibers;
		string fileName;
		bool needSwap;
		bool isASCII;
		bool hasPointColor;
		bool hasLineColor;
		bool hasLabel;
		float resolution; 
		
	};


	CFibers::CFibers(string fileName )
	{
		m_data= new CFibersPrivacy;
		m_data->fileName=string(fileName);
		m_data->needSwap=true;
		m_data->hasPointColor=true;
		m_data->hasLabel=false;

		fstream headerStrm;
		OpenReadStrmAscii(headerStrm,fileName);
		string modeString;
		getline(headerStrm,modeString);
		getline(headerStrm,modeString);
		getline(headerStrm,modeString);

		if(modeString.find("ASCII")!=string::npos  || modeString.find("ascii")!=string::npos){
			m_data->isASCII=true;
		}else{
			m_data->isASCII=false;
		}

		headerStrm.close();

 


		fstream fiberFileStrm;

		if(m_data->isASCII){
			OpenReadStrmAscii(fiberFileStrm,fileName);

			//read points;
			string line;
			do
			{
				getline(fiberFileStrm,line);
			} while (line.find("POINTS")==string::npos);

			stringstream lineStrm;
			lineStrm.str(line);
			lineStrm>>line>>m_data->numPoints;
//			cout<<m_data->numPoints<<" fiber points "<<"found. "<<endl;

			m_data->allPoints.resize(m_data->numPoints);

			for (int pointIndex=0; pointIndex<m_data->numPoints; ++pointIndex)
			{
				fiberFileStrm>>m_data->allPoints[pointIndex].x>>m_data->allPoints[pointIndex].y>>m_data->allPoints[pointIndex].z;
			}

			// read fibers;
			do
			{
				getline(fiberFileStrm,line);
			} while (line.find("LINES")==string::npos);

			lineStrm.str(line);
			lineStrm>>line>>m_data->numFibers;
//			cout<<m_data->numFibers<<" fibers "<<"found. "<<endl;
			m_data->allFibers.resize(m_data->numFibers);
			for(int fiberIndex=0; fiberIndex<m_data->numFibers; ++fiberIndex)
			{
				int currentFiberSize=0;
				fiberFileStrm>>currentFiberSize;
				m_data->allFibers[fiberIndex].resize(currentFiberSize);
				for(int fiberPointIndex=0; fiberPointIndex<currentFiberSize; ++fiberPointIndex)
				{
					fiberFileStrm>>m_data->allFibers[fiberIndex][fiberPointIndex];
				}
			}


			//read fiber colors;
		 	do
			{
				getline(fiberFileStrm,line);
			} while (line.find("POINT_DATA")==string::npos && !fiberFileStrm.eof());

		 	if(fiberFileStrm.eof()){
		 		m_data->hasPointColor=false;
		 	}
		 	if(m_data->hasPointColor){
			getline(fiberFileStrm,line); // remove the following line; color_scalars

			m_data->allPointsColorsAscii.resize(m_data->numPoints);
			for(int pointIndex=0; pointIndex<m_data->numPoints; ++pointIndex)
			{
				fiberFileStrm>>m_data->allPointsColorsAscii[pointIndex].x>>m_data->allPointsColorsAscii[pointIndex].y
				>>m_data->allPointsColorsAscii[pointIndex].z;
			}
			
					 	
		 	//copy color data;
			m_data->allPointsColorsBinary.resize(m_data->numPoints);
			for (int index = 0; index < m_data->numPoints; ++index) {
				m_data->allPointsColorsBinary[index].x= (unsigned char)(m_data->allPointsColorsAscii[index].x*255.);
				m_data->allPointsColorsBinary[index].y= (unsigned char)(m_data->allPointsColorsAscii[index].y*255.);
				m_data->allPointsColorsBinary[index].z= (unsigned char)(m_data->allPointsColorsAscii[index].z*255.);
			} //end loop:: index
			
		 	}


		} else { // binary mode;
			OpenReadStrmBinary(fiberFileStrm,fileName);
			char tmpLine[80];
			while(!fiberFileStrm.eof()){
				fiberFileStrm.getline(tmpLine,80);
				string lineString(tmpLine);
				if(lineString.find("POINTS")!=string::npos){
					strstream lineStrm;
					lineStrm<<lineString;
					lineStrm>>lineString>>m_data->numPoints;
					break;
			    } // end of if;
			}

			//read points;
			m_data->allPoints.resize(m_data->numPoints);
			fiberFileStrm.read((char*) (& (m_data->allPoints[0].x)),sizeof(float)*m_data->numPoints*3);

			//read fibers;
			{
				while(!fiberFileStrm.eof()){
					fiberFileStrm.getline(tmpLine,80);
					string lineString(tmpLine);
					if(lineString.find("LINES")!=string::npos){
						strstream lineStrm;
						lineStrm<<lineString;
						lineStrm>>lineString>>m_data->numFibers;
						break;
				    } // end of if;
				}// end of while;

				m_data->allFibers.resize(m_data->numFibers);
				for (int fiberIndex = 0; fiberIndex < m_data->numFibers; ++fiberIndex) {
					int currentFiberSize;
					fiberFileStrm.read((char*) &currentFiberSize, sizeof(int));
					if (m_data->needSwap)
						ByteSwapAllType(currentFiberSize);
					if(0> currentFiberSize || currentFiberSize> 10000000){
						cout<<"Byte order error, you need to change the swap setting."<<endl;
					    exit(EXIT_FAILURE);
					}
					m_data->allFibers[fiberIndex].resize(currentFiberSize);
					fiberFileStrm.read((char*) (&(m_data->allFibers[fiberIndex][0])),sizeof(int)*(currentFiberSize));
				} //end of for loop:: fiberIndex

			} // read fibers end

			//read fiber data

			m_data->allPointsColorsBinary.resize(m_data->numPoints);
			while(!fiberFileStrm.eof()){
					fiberFileStrm.getline(tmpLine,80);
					string lineString(tmpLine);
					if(lineString.find("DATA")!=string::npos){
						fiberFileStrm.getline(tmpLine,80); // remove the next line;
						if(lineString.find("CELL_DATA")!=string::npos){
							int sizeofElements=m_data->numFibers*3;
							char* tmpDst= new char[sizeofElements];
							fiberFileStrm.read(tmpDst,sizeofElements);
							delete tmpDst;
						}
						if(lineString.find("POINT_DATA")!=string::npos){
							break;
						};
					};
			};

			fiberFileStrm.read((char*) (& (m_data->allPointsColorsBinary[0].x)), sizeof(unsigned char)*3*m_data->numPoints);

			if (! m_data->isASCII  && m_data->needSwap)
				SwapEndians();
//			cout<<"First point coordinate: "<<m_data->allPoints[0].x<<" "<<m_data->allPoints[0].y<<" "<<m_data->allPoints[0].z<<endl;
//			cout<<"First fiber : ";
//			for (int pointid = 0; pointid < m_data->allFibers[0].size(); ++pointid) {
//				cout<< m_data->allFibers[0][pointid]<<" ";
//			} //end of for loop:: pointid

//			cout<<endl<<"Second fiber : ";
//			for (int pointid = 0; pointid < m_data->allFibers[1].size(); ++pointid) {
//				cout<< m_data->allFibers[1][pointid]<<" ";
//			} //end of for loop:: pointid

//			cout<<endl<<"If not correct, you need to call SwapEndians()"<<endl;

// 			cout<<"Fiber points: "<<m_data->numPoints<<endl;
// 			cout<<"Fiber lines: "<<m_data->numFibers<<endl;
			//copy color data;
			m_data->allPointsColorsAscii.resize(m_data->numPoints);
			for (int index = 0; index < m_data->numPoints; ++index) {
				m_data->allPointsColorsAscii[index].x=m_data->allPointsColorsBinary[index].x/255.;
				m_data->allPointsColorsAscii[index].y=m_data->allPointsColorsBinary[index].y/255.;
				m_data->allPointsColorsAscii[index].z=m_data->allPointsColorsBinary[index].z/255.;
			} //end loop:: index

		}
		fiberFileStrm.close();
		m_data->resolution=this->GetResolution();
	}  //end of constructor; 
	CFibers::CFibers(CFibers& fibers){
		m_data= new CFibersPrivacy;
		m_data->allFibers= fibers.m_data->allFibers;
		m_data->allPoints=fibers.m_data->allPoints;
		m_data->allPointsColorsAscii=fibers.m_data->allPointsColorsAscii;
		m_data->fiberLabel=fibers.m_data->fiberLabel;
		m_data->fileName=fibers.m_data->fileName;
		m_data->hasPointColor=fibers.m_data->hasPointColor;
		m_data->hasLineColor= fibers.m_data->hasLineColor;
		m_data->hasLabel=fibers.m_data->hasLabel;
		m_data->isASCII=fibers.m_data->isASCII;
		m_data->needSwap=fibers.m_data->needSwap;
		m_data->numFibers=fibers.m_data->numFibers;
		m_data->numPoints=fibers.m_data->numPoints;
		m_data->resolution=fibers.m_data->resolution; 

	}
	CFibers::~CFibers()
	{
		delete m_data; 
	} //end of deconstructor; 

	void CFibers::IsSwapNeeded(bool needSwap){
		m_data->needSwap=needSwap;
	}
	vector<int>& CFibers::GetFiber( int fiberId )
	{
		if (fiberId<0 ||(fiberId>=m_data->numFibers) )
		{
			cerr<<"Wrong fiber id: "<<fiberId<<endl;
			cerr<<"Should be within the range of "<<0<<" and "<<m_data->numFibers; 
			exit(EXIT_FAILURE);
		}
		return m_data->allFibers[fiberId];
	}

	PointCoordType& CFibers::GetPointCoords( int pointId )
	{
		if (pointId<0 || (pointId>=m_data->numPoints) )
		{
			cerr<<"Wrong point id: "<<pointId<<endl;
			cerr<<"Should be within the range of "<<0<<" and "<<m_data->numPoints; 
			exit(EXIT_FAILURE);
		}

		return m_data->allPoints[pointId];
	}

	void CFibers::SetPointCoords( int pointId, PointCoordType& newCoords )
	{

		if (pointId<0 || (pointId>=m_data->numPoints) )
		{
			cerr<<"Wrong point id: "<<pointId<<endl;
			cerr<<"Should be within the range of "<<0<<" and "<<m_data->numPoints; 
			exit(EXIT_FAILURE);
		}

		m_data->allPoints[pointId]=newCoords;

	}

	PointCoordType& CFibers::GetPointColorAscii(int pointId)
	{
	  return m_data->allPointsColorsAscii[pointId]; 
	}

	Vector3D<unsigned char>& CFibers::GetPointColorBinary(int pointId)
	{
	  return m_data->allPointsColorsBinary[pointId]; 
	}

	int CFibers::GetNumFibers( void )
	{
		return m_data->numFibers;
	}

	const string& CFibers::GetFiberName() const{
		return m_data->fileName;
	}
	float CFibers::GetFiberLength( int fiberId ) const
	{
		vector<int>& currentFiber=m_data->allFibers[fiberId];
		float length=0;
		for (int fiberPointIndex=0; fiberPointIndex<currentFiber.size()-1; ++fiberPointIndex)
		{
			int pointId1=currentFiber[fiberPointIndex];
			int pointId2=currentFiber[fiberPointIndex+1];
			length+=DistanceVector3D(m_data->allPoints[pointId1],m_data->allPoints[pointId2]);
		}
		return length;
	}

	int CFibers::GetFiberSize( int fiberId ) const
	{
		return m_data->allFibers[fiberId].size();
	}

	void CFibers::SwapEndians(){
		//swap points coordinate;
		for (int pointIndex = 0; pointIndex < m_data->numPoints; ++pointIndex) {
			ByteSwapAllType(m_data->allPoints[pointIndex].x);
			ByteSwapAllType(m_data->allPoints[pointIndex].y);
			ByteSwapAllType(m_data->allPoints[pointIndex].z);
		} //end of for loop:: pointIndex

		//swap fibers;
		for (int fiberIndex = 0; fiberIndex < m_data->numFibers; ++fiberIndex) {
			for (int pointIndex = 0; pointIndex < m_data->allFibers[fiberIndex].size(); ++pointIndex) {
				ByteSwapAllType(m_data->allFibers[fiberIndex][pointIndex]);
			} //end of for loop:: pointIndex
		} //end of for loop:: fiberIndex

		//swap fibers colors;
		for (int pointId = 0; pointId < m_data->numPoints; ++pointId) {
			ByteSwapAllType(m_data->allPointsColorsBinary[pointId].x);
			ByteSwapAllType(m_data->allPointsColorsBinary[pointId].y);
			ByteSwapAllType(m_data->allPointsColorsBinary[pointId].z);
		} //end of for loop:: pointId

	}
	
	

	int CFibers::GetNumPoints(){
		return m_data->numPoints;
	}

	void CFibers::SaveAllFibers(string fileName, bool toAsciiMode){

		fstream fiberStrm;
		if(toAsciiMode)
			OpenWriteStrmAscii(fiberStrm,fileName.c_str());
		else{
			OpenWriteStrmBinary(fiberStrm,fileName.c_str());
		}


		fiberStrm<<"# vtk DataFile Version 2.0"<<endl;
		fiberStrm<<"all fibers"<<endl;
		if(toAsciiMode){
			fiberStrm<<"ASCII"<<endl;
		}else{
			fiberStrm<<"BINARY"<<endl;
		}
		fiberStrm<<"DATASET POLYDATA"<<endl;
		fiberStrm<<"POINTS "<<m_data->numPoints<<" float"<<endl;

		//get line size;
		int allLineSize=0;
		for (int fiberId = 0; fiberId < m_data->numFibers; ++fiberId) {
			allLineSize+=m_data->allFibers[fiberId].size();
		} //end of for loop:: fiberId
		allLineSize+=m_data->allFibers.size();

		if(m_data->isASCII){
			if(toAsciiMode){
				//write points;
				for (int pointID = 0; pointID < m_data->numPoints; ++pointID) {
					fiberStrm<<m_data->allPoints[pointID].x<<" "<<m_data->allPoints[pointID].y<<" "<<m_data->allPoints[pointID].z<<endl;
				} //end of for loop:: pointID

				//write lines;
				fiberStrm<<"LINES "<<m_data->numFibers<<" "<<allLineSize<<endl;
				for (int fiberId = 0; fiberId < m_data->numFibers; ++fiberId) {
					fiberStrm<<m_data->allFibers[fiberId].size()<<" ";
					for (int currentFiberElement = 0; currentFiberElement < m_data->allFibers[fiberId].size(); ++currentFiberElement) {
						fiberStrm<<m_data->allFibers[fiberId][currentFiberElement]<<" ";
					} //end of for loop:: currentFiberElement
					fiberStrm<<endl;
				} //end of for loop:: fiberId

				//write point data;
				if(m_data->hasPointColor){
				fiberStrm<<"POINT_DATA "<<m_data->numPoints<<endl;
				fiberStrm<<"COLOR_SCALARS  pointColor 3"<<endl;
				for (int pointID = 0; pointID < m_data->numPoints; ++pointID) {
					fiberStrm<<m_data->allPointsColorsAscii[pointID].x<<" "<<m_data->allPointsColorsAscii[pointID].y<<" "<<m_data->allPointsColorsAscii[pointID].z<<endl;
				} //end of for loop:: pointID
				}
		    }else{ //ascii->binary;

				//whether need to convert float to unsigned char first;
				bool needConvert=true;
				for (int pointID = 0; pointID < m_data->numPoints; ++pointID) {
					if ( m_data->allPointsColorsAscii[pointID].x>1) {
						needConvert=false;
						break;
					} // end of if ::  m_data->allPointsColorsAscii[pointID].x>1
				} //end of for loop:: pointID

				m_data->allPointsColorsBinary.resize(m_data->numPoints);
				if(needConvert){
					for (int pointID = 0; pointID < m_data->numPoints; ++pointID) {
						m_data->allPointsColorsBinary[pointID].x=(unsigned char)(m_data->allPointsColorsAscii[pointID].x*255.);
						m_data->allPointsColorsBinary[pointID].y=(unsigned char)(m_data->allPointsColorsAscii[pointID].y*255.);
						m_data->allPointsColorsBinary[pointID].z=(unsigned char)(m_data->allPointsColorsAscii[pointID].z*255.);
					} //end of for loop:: pointID
				}else{
					for (int pointID = 0; pointID < m_data->numPoints; ++pointID) {
						m_data->allPointsColorsBinary[pointID].x= m_data->allPointsColorsAscii[pointID].x;
						m_data->allPointsColorsBinary[pointID].y= m_data->allPointsColorsAscii[pointID].y;
						m_data->allPointsColorsBinary[pointID].z= m_data->allPointsColorsAscii[pointID].z;
					} //end of for loop:: pointID
				}


		    this->SwapEndians();
			//write points;
			fiberStrm.write((const char*) & (m_data->allPoints[0].x),sizeof(float)*3*m_data->numPoints);
			//write LINES;
			fiberStrm<<endl<<"LINES "<<m_data->numFibers<<" "<<allLineSize<<endl;
			for (int fiberId = 0; fiberId < m_data->numFibers; ++fiberId) {
				int currentFiberSize=m_data->allFibers[fiberId].size();
				ByteSwapAllType(currentFiberSize);
				fiberStrm.write( (const char*) (&currentFiberSize), sizeof(int));
				fiberStrm.write((const char*) (& m_data->allFibers[fiberId][0]), sizeof(int)*m_data->allFibers[fiberId].size());
			} //end of for loop:: fiberId

			//write point data;
			fiberStrm<<endl<<"POINT_DATA "<<m_data->numPoints<<endl;
			fiberStrm<<"COLOR_SCALARS  pointColor 3"<<endl;
			fiberStrm.write((const char*) &(m_data->allPointsColorsBinary[0].x),sizeof(unsigned char)*3*m_data->numPoints);
			//swap back;
			 this->SwapEndians();
	    	} //end of if-else:: toAsciiMode

		}else{ // original is binary;
			if(toAsciiMode){
				//write points;
				for (int pointID = 0; pointID < m_data->numPoints; ++pointID) {
					fiberStrm<<m_data->allPoints[pointID].x<<" "<<m_data->allPoints[pointID].y<<" "<<m_data->allPoints[pointID].z<<endl;
				} //end of for loop:: pointID

				//write lines;
				fiberStrm<<"LINES "<<m_data->numFibers<<" "<<allLineSize<<endl;
				for (int fiberId = 0; fiberId < m_data->numFibers; ++fiberId) {
					fiberStrm<<m_data->allFibers[fiberId].size()<<" ";
					for (int currentFiberElement = 0; currentFiberElement < m_data->allFibers[fiberId].size(); ++currentFiberElement) {
						fiberStrm<<m_data->allFibers[fiberId][currentFiberElement]<<" ";
					} //end of for loop:: currentFiberElement
					fiberStrm<<endl;
				} //end of for loop:: fiberId

				//write point data;
				fiberStrm<<"POINT_DATA "<<m_data->numPoints<<endl;
				fiberStrm<<"COLOR_SCALARS  pointColor 3"<<endl;
				for (int pointID = 0; pointID < m_data->numPoints; ++pointID) {
					fiberStrm<<m_data->allPointsColorsBinary[pointID].x/255.<<" "<<m_data->allPointsColorsBinary[pointID].y/255.
							<<" "<<m_data->allPointsColorsBinary[pointID].z/255.<<endl;
				} //end of for loop:: pointID
			}else{
				this->SwapEndians(); //paraview recognize big endian;
				//write points;
				fiberStrm.write((const char*) & (m_data->allPoints[0].x),sizeof(float)*3*m_data->numPoints);
				//write LINES;
				fiberStrm<<endl;
				fiberStrm<<"LINES "<< m_data->numFibers <<" "<< allLineSize<<endl;
				for (int fiberId = 0; fiberId < m_data->numFibers; ++fiberId) {
					int currentFiberSize=m_data->allFibers[fiberId].size();
					ByteSwapAllType(currentFiberSize);
					fiberStrm.write( (const char*) & ( currentFiberSize), sizeof(int));
					fiberStrm.write((const char*) & (m_data->allFibers[fiberId][0]), sizeof(int)*m_data->allFibers[fiberId].size());
				} //end of for loop:: fiberId

                //write point data;
				fiberStrm<<endl<<"POINT_DATA "<<m_data->numPoints<<endl;
				fiberStrm<<"COLOR_SCALARS  pointColor 3"<<endl;
				fiberStrm.write((const char*) &(m_data->allPointsColorsBinary[0].x),sizeof(unsigned char)*3*m_data->numPoints);

				//swap back;
				this->SwapEndians();
			}
		}
		fiberStrm.close();
	}

	void CFibers::SaveFibersWithLabel(string fileName, int labelValue){

		cout<<"Saving fibers labeled with "<<labelValue<<"...";
		fstream fiberStrm;
		OpenWriteStrmAscii(fiberStrm,fileName.c_str());

		fiberStrm<<"# vtk DataFile Version 2.0"<<endl;
		fiberStrm<<"all fibers"<<endl;
		 fiberStrm<<"ASCII"<<endl;
		fiberStrm<<"DATASET POLYDATA"<<endl;
		fiberStrm<<"POINTS "<<m_data->numPoints<<" float"<<endl;

		//get line size;
		int allLineSize=0;
		int linesWithLabelvalue=0;
		for (int fiberId = 0; fiberId < m_data->numFibers; ++fiberId) {
			if(labelValue==m_data->fiberLabel[fiberId]){
				allLineSize+=m_data->allFibers[fiberId].size();
				++linesWithLabelvalue;}
		} //end of for loop:: fiberId
		allLineSize+=linesWithLabelvalue;

				//write points;
				for (int pointID = 0; pointID < m_data->numPoints; ++pointID) {
					fiberStrm<<m_data->allPoints[pointID].x<<" "<<m_data->allPoints[pointID].y<<" "<<m_data->allPoints[pointID].z<<endl;
				} //end of for loop:: pointID

				//write lines;
				fiberStrm<<"LINES "<<linesWithLabelvalue<<" "<<allLineSize<<endl;
				for (int fiberId = 0; fiberId < m_data->numFibers; ++fiberId) {

					if (labelValue==m_data->fiberLabel[fiberId]) {
						fiberStrm<<m_data->allFibers[fiberId].size()<<" ";
						for (int currentFiberElement = 0; currentFiberElement < m_data->allFibers[fiberId].size(); ++currentFiberElement) {
							fiberStrm<<m_data->allFibers[fiberId][currentFiberElement]<<" ";
						} //end of for loop:: currentFiberElement
					} // end of if :: labelValue==m_data->fiberLabel[fiberId]
					fiberStrm<<endl;
				} //end of for loop:: fiberId

				//write point data;
				if(m_data->hasPointColor){
				fiberStrm<<"POINT_DATA "<<m_data->numPoints<<endl;
				fiberStrm<<"COLOR_SCALARS  pointColor 3"<<endl;
				for (int pointID = 0; pointID < m_data->numPoints; ++pointID) {
					fiberStrm<<m_data->allPointsColorsAscii[pointID].x<<" "<<m_data->allPointsColorsAscii[pointID].y<<" "<<m_data->allPointsColorsAscii[pointID].z<<endl;
				} //end of for loop:: pointID
		}
		fiberStrm.close();
		cout<<"done."<<endl;
	}
	void CFibers::SaveFibers(list<int>& fiberIds, string fileName, bool toAsciiMode){
 
	    vector<VectorType> newAllPoints; 
	    vector<vector<int> > newAllLines; 
	    vector<PointCoordType> newAllPointsColorsAscii;
	    vector<Vector3D<unsigned char> > newAllPointsColorsBinary;
	    
	    int newPid=0; 
	    
	    for(list<int>::iterator it= fiberIds.begin(); it != fiberIds.end(); it++)
	    {
	      int idLine= *it; 
	      vector<int>& currentLine= m_data->allFibers[idLine]; 
	      {
                vector<int> newLine; 
		for (int newEleIndex = 0; newEleIndex <  currentLine.size(); ++newEleIndex) {
		  int originalId= currentLine[newEleIndex]; 
		  newAllPoints.push_back(m_data->allPoints[originalId]);
		  newAllPointsColorsAscii.push_back(m_data->allPointsColorsAscii[originalId]);
		  newAllPointsColorsBinary.push_back(m_data->allPointsColorsBinary[originalId]);
		  newLine.push_back(newPid++);
		}//end of for loop::newEleIndex
		newAllLines.push_back(newLine);
	      }
	    }//end of for loop::idLine
	 
 

	    fstream fiberStrm;
	    if(toAsciiMode)
		    OpenWriteStrmAscii(fiberStrm,fileName.c_str());
	    else{
		    OpenWriteStrmBinary(fiberStrm,fileName.c_str());
	    }


	    fiberStrm<<"# vtk DataFile Version 2.0"<<endl;
	    fiberStrm<<"all fibers"<<endl;
	    if(toAsciiMode){
		    fiberStrm<<"ASCII"<<endl;
	    }else{
		    fiberStrm<<"BINARY"<<endl;
	    }
	    fiberStrm<<"DATASET POLYDATA"<<endl;
	    fiberStrm<<"POINTS "<< newAllPoints.size()<<" float"<<endl;

	    //get line size;
	    int allLineSize=0;
	    for (int fiberId = 0; fiberId < newAllLines.size(); ++fiberId) {
		    allLineSize+=newAllLines[fiberId].size();
	    } //end of for loop:: fiberId
	    allLineSize+=newAllLines.size();

		if(m_data->isASCII){
			if(toAsciiMode){
				//write points;
				for (int pointID = 0; pointID < newAllPoints.size(); ++pointID) {
					fiberStrm<<newAllPoints[pointID].x<<" "<<newAllPoints[pointID].y<<" "<<newAllPoints[pointID].z<<endl;
				} //end of for loop:: pointID

				//write lines;
				fiberStrm<<"LINES "<<newAllLines.size()<<" "<<allLineSize<<endl;
				for (int fiberId = 0; fiberId < newAllLines.size(); ++fiberId) {
					fiberStrm<<newAllLines[fiberId].size()<<" ";
					for (int currentFiberElement = 0; currentFiberElement < newAllLines[fiberId].size(); ++currentFiberElement) {
						fiberStrm<<newAllLines[fiberId][currentFiberElement]<<" ";
					} //end of for loop:: currentFiberElement
					fiberStrm<<endl;
				} //end of for loop:: fiberId

				//write point data;
				if(m_data->hasPointColor){
				fiberStrm<<"POINT_DATA "<<newAllPoints.size()<<endl;
				fiberStrm<<"COLOR_SCALARS  pointColor 3"<<endl;
				for (int pointID = 0; pointID < newAllPoints.size(); ++pointID) {
					fiberStrm<<newAllPointsColorsAscii[pointID].x<<" "<<newAllPointsColorsAscii[pointID].y<<" "<<newAllPointsColorsAscii[pointID].z<<endl;
				} //end of for loop:: pointID
				}
		    }else{ //ascii->binary;

// 			this->SwapEndians();
		    				for (int pointIndex = 0; pointIndex < newAllPoints.size(); ++pointIndex) {
					ByteSwapAllType(newAllPoints[pointIndex].x);
					ByteSwapAllType(newAllPoints[pointIndex].y);
					ByteSwapAllType(newAllPoints[pointIndex].z);
				} //end of for loop:: pointIndex

				//swap fibers;
				for (int fiberIndex = 0; fiberIndex < newAllLines.size(); ++fiberIndex) {
					for (int pointIndex = 0; pointIndex < newAllLines[fiberIndex].size(); ++pointIndex) {
						ByteSwapAllType(newAllLines[fiberIndex][pointIndex]);
					} //end of for loop:: pointIndex
				} //end of for loop:: fiberIndex

				//swap fibers colors;
				for (int pointId = 0; pointId < newAllPoints.size(); ++pointId) {
					ByteSwapAllType(newAllPointsColorsBinary[pointId].x);
					ByteSwapAllType(newAllPointsColorsBinary[pointId].y);
					ByteSwapAllType(newAllPointsColorsBinary[pointId].z);
				} //end of for loop:: pointId1
			//write points;
			fiberStrm.write((const char*) & (newAllPoints[0].x),sizeof(float)*3*newAllPoints.size());
			//write LINES;
			fiberStrm<<endl<<"LINES "<<newAllLines.size()<<" "<<allLineSize<<endl;
			for (int fiberId = 0; fiberId < newAllLines.size(); ++fiberId) {
				int currentFiberSize=newAllLines[fiberId].size();
				ByteSwapAllType(currentFiberSize);
				fiberStrm.write( (const char*) (&currentFiberSize), sizeof(int));
				fiberStrm.write((const char*) (& newAllLines[fiberId][0]), sizeof(int)*newAllLines[fiberId].size());
			} //end of for loop:: fiberId

			//write point data;
			fiberStrm<<endl<<"POINT_DATA "<<newAllPoints.size()<<endl;
			fiberStrm<<"COLOR_SCALARS  pointColor 3"<<endl;
			fiberStrm.write((const char*) &(newAllPointsColorsBinary[0].x),sizeof(unsigned char)*3*newAllPoints.size());
			//swap back;
// 			 this->SwapEndians();
	    	} //end of if-else:: toAsciiMode

		}else{ // original is binary;
			if(toAsciiMode){
				//write points;
				for (int pointID = 0; pointID < newAllPoints.size(); ++pointID) {
					fiberStrm<<newAllPoints[pointID].x<<" "<<newAllPoints[pointID].y<<" "<<newAllPoints[pointID].z<<endl;
				} //end of for loop:: pointID

				//write lines;
				fiberStrm<<"LINES "<<newAllLines.size()<<" "<<allLineSize<<endl;
				for (int fiberId = 0; fiberId < newAllLines.size(); ++fiberId) {
					fiberStrm<<newAllLines[fiberId].size()<<" ";
					for (int currentFiberElement = 0; currentFiberElement < newAllLines[fiberId].size(); ++currentFiberElement) {
						fiberStrm<<newAllLines[fiberId][currentFiberElement]<<" ";
					} //end of for loop:: currentFiberElement
					fiberStrm<<endl;
				} //end of for loop:: fiberId

				//write point data;
				fiberStrm<<"POINT_DATA "<<newAllPoints.size()<<endl;
				fiberStrm<<"COLOR_SCALARS  pointColor 3"<<endl;
				for (int pointID = 0; pointID < newAllPoints.size(); ++pointID) {
					fiberStrm<<newAllPointsColorsBinary[pointID].x/255.<<" "<<newAllPointsColorsBinary[pointID].y/255.
							<<" "<<newAllPointsColorsBinary[pointID].z/255.<<endl;
				} //end of for loop:: pointID
			}else{
				//this->SwapEndians(); //paraview recognize big endian;
				
				for (int pointIndex = 0; pointIndex < newAllPoints.size(); ++pointIndex) {
					ByteSwapAllType(newAllPoints[pointIndex].x);
					ByteSwapAllType(newAllPoints[pointIndex].y);
					ByteSwapAllType(newAllPoints[pointIndex].z);
				} //end of for loop:: pointIndex

				//swap fibers;
				for (int fiberIndex = 0; fiberIndex < newAllLines.size(); ++fiberIndex) {
					for (int pointIndex = 0; pointIndex < newAllLines[fiberIndex].size(); ++pointIndex) {
						ByteSwapAllType(newAllLines[fiberIndex][pointIndex]);
					} //end of for loop:: pointIndex
				} //end of for loop:: fiberIndex

				//swap fibers colors;
				for (int pointId = 0; pointId < newAllPoints.size(); ++pointId) {
					ByteSwapAllType(newAllPointsColorsBinary[pointId].x);
					ByteSwapAllType(newAllPointsColorsBinary[pointId].y);
					ByteSwapAllType(newAllPointsColorsBinary[pointId].z);
				} //end of for loop:: pointId1
		
				//write points;
				fiberStrm.write((const char*) & (newAllPoints[0].x),sizeof(float)*3*newAllPoints.size());
				//write LINES;
				fiberStrm<<endl;
				fiberStrm<<"LINES "<< newAllLines.size() <<" "<< allLineSize<<endl;
				for (int fiberId = 0; fiberId < newAllLines.size(); ++fiberId) {
					int currentFiberSize=newAllLines[fiberId].size();
					ByteSwapAllType(currentFiberSize);
					fiberStrm.write( (const char*) & ( currentFiberSize), sizeof(int));
					fiberStrm.write((const char*) & (newAllLines[fiberId][0]), sizeof(int)*newAllLines[fiberId].size());
				} //end of for loop:: fiberId

				//write point data;
				fiberStrm<<endl<<"POINT_DATA "<<newAllPoints.size()<<endl;
				fiberStrm<<"COLOR_SCALARS  pointColor 3"<<endl;
				fiberStrm.write((const char*) &(newAllPointsColorsBinary[0].x),sizeof(unsigned char)*3*newAllPoints.size());

				//swap back;
				//this->SwapEndians();
			}
		}
		fiberStrm.close();


	}
	void CFibers::SaveFibers(set<int>& fiberIds, string fileName){
		list<int> allFibers;
		for(set<int>::iterator itAllFibers= fiberIds.begin(); itAllFibers!=fiberIds.end(); itAllFibers++){
			allFibers.push_back(*itAllFibers);
		}
		this->SaveFibers(allFibers,fileName);
	}
	void CFibers::SaveFibers( const vector<int>& fiberIds, string fileName){
		list<int> allFibers( fiberIds.size());
		copy( fiberIds.begin(), fiberIds.end(), allFibers.begin());
		this->SaveFibers( allFibers, fileName);
	}

	void CFibers::SaveFibers(const set<int>& fiberIds, string fileName){
		list<int> allFibers;
		for(set<int>::const_iterator itAllFibers= fiberIds.begin(); itAllFibers!=fiberIds.end(); itAllFibers++){
			allFibers.push_back(*itAllFibers);
		}
		this->SaveFibers(allFibers,fileName);
	}
	void CFibers::ReadFiberLabel(string fileName){
		// works  only in ascii mode;
		if(true!=m_data->isASCII){
			cout<<"Error! Only supports ascii mode now."<<endl;
			exit(EXIT_FAILURE);
		}

		fstream fiberStrm;
		OpenReadStrmAscii(fiberStrm,m_data->fileName);
		string line;
		while(getline(fiberStrm,line)){
			if(line.find(fileName)!=string::npos){
				getline(fiberStrm,line); // remove the LOOKUP_TABLE default line;
				break;
			}
		} //end of while;

		m_data->fiberLabel.resize(m_data->numFibers);
		for (int fiberId = 0; fiberId < m_data->numFibers; ++fiberId) {
			fiberStrm>>m_data->fiberLabel[fiberId];
		} //end of for loop:: fiberId
	}

	vector<int>& CFibers::GetALLFiberLabel(void){
		return m_data->fiberLabel;
	}

	int& CFibers::GetFiberLabel(int fiberId){
		return m_data->fiberLabel[fiberId];
	}
	void CFibers::SetFiberLabel(int fiberId,int fiberLabel){
		m_data->fiberLabel[fiberId]=fiberLabel;
	}
	void CFibers::SetFiberLabel(vector<int>& allLabels){
		assert(allLabels.size()==m_data->fiberLabel.size());
		m_data->fiberLabel=allLabels;
	}
	void CFibers::Translate(float x, float y, float z){

		for (int pointId = 0; pointId < m_data->numPoints; ++pointId) {
			m_data->allPoints[pointId].x+=x;
			m_data->allPoints[pointId].y+=y;
			m_data->allPoints[pointId].z+=z;
		} //end of for loop:: pointId
	}

	void CFibers::ShortenFibers(int numPoints)
	{
		for (int fiberID = 0; fiberID < m_data->numFibers; ++fiberID) {
			vector<int> & currentFiber= m_data->allFibers[fiberID];
			vector<int> newFiber;
			if(currentFiber.size() > 2*numPoints)
			{
				for (int elemIndex = numPoints; elemIndex < currentFiber.size()-numPoints; ++elemIndex) {
					newFiber.push_back(currentFiber[elemIndex]);
				} //end of for loop:: elemIndex
				currentFiber= newFiber;
			}

		} //end of for loop:: fiberID

	}
// 	void CFibers::ShortenFibers(CTriSurface& mySurface, float maxDis)
// 	{
// 		CJointModel myModel( mySurface, *this);
// 		myModel.FindFiberSurfIntersections();
// 		myModel.FindFibersOnCells();
// 		list<int> allInvalidFiber;
// 		vector<bool> isFiberIntersectedWithSurf(this->GetNumFibers(),false);
// 
// 		for each cell on the cortical surface, find the intersecion fibers;
// 		for each fiber , if it has two intersection points; remove the first and last points ;
// 		if no intersections or only one intersection, then discard the fiber; that is to remove the fiber, or copy the fiber from a valid one;
// 		for (int cellID = 0; cellID < mySurface.GetNumOfCells(); ++cellID) {
// 
// 			vector<int>&  allInterFibersCurrentCell = myModel.GetIntersectingFibersOfCell(cellID);
// 			if(allInterFibersCurrentCell.size()){
// 				for( int fiberIndex=0; fiberIndex< allInterFibersCurrentCell.size(); ++fiberIndex){
// 					int fid= allInterFibersCurrentCell[fiberIndex];
// 					isFiberIntersectedWithSurf[fid]=true;
// 					CIntersectionType& frontIntersection= myModel.GetFrontIntersectingCellOfFiber(fid);
// 					CIntersectionType& rearIntersection= myModel.GetRearIntersectingCellOfFiber(fid);
// 
// 					for each fiber , if it has two intersection points; remove the first and last points ;
// 					if( frontIntersection.IsValid() && rearIntersection.IsValid() ){
// 						remove the first and last fiber elements;
// 						this->ShortenCurrentFiber(fid, frontIntersection,rearIntersection,maxDis);
// 					}
// 					else
// 					{
// 						allInvalidFiber.push_back(fid);
// 					}
// 				}
// 			}
// 
// 		} //end of for loop:: cellID
// 
// 		remove all the invalid fibers;
// 		for (int index = 0; index < isFiberIntersectedWithSurf.size(); ++index) {
// 			if(false==isFiberIntersectedWithSurf[index])
// 				allInvalidFiber.push_back(index);
// 		} //end of for loop:: index
// 
// 		cout<<"remove invalid fibers: "<<allInvalidFiber.size();
// 		for( list<int>::iterator it=allInvalidFiber.begin(); it!= allInvalidFiber.end(); ++it){
// 			int fid= *it;
// 			m_data->allFibers[fid].clear();
// 		}
// 
// 	}
// 
// 	int CFibers::ShortenCurrentFiber(int fid, CIntersectionType& front, CIntersectionType& rear, float tolerance){
// 
// 		 vector<int>& allFiberElems= m_data->allFibers[fid];
// 		 Vector3D<float> frontCoord= front.GetCoord();
// 		 Vector3D<float> rearCoord= rear.GetCoord();
// 
// 		 int frontID=0;
// 		 int rearID=0;
// 
// 		 find the two elements on intersections;
// 		 for (int index = 0; index < allFiberElems.size(); ++index) {
// 			 int fid= allFiberElems[index];
// 			 Vector3D<float>& thePoint= this->m_data->allPoints[fid];
// 			 Vector3D<float> disFront= thePoint-frontCoord;
// 			 Vector3D<float> disRear= thePoint-rearCoord;
//  			 if(disFront.Norm() <tolerance)
// 			 {
// 				 break;
// 			 }
//  			 ++frontID;
// 		} //end of for loop:: index
// 
// 		 find the two elements on intersections;
// 		 rearID= allFiberElems.size();
// 		 for (vector<int>::reverse_iterator it=allFiberElems.rbegin(); it!=allFiberElems.rend(); ++it) {
// 			 int fid= *it;
// 			 Vector3D<float>& thePoint= this->m_data->allPoints[fid];
// 			 Vector3D<float> disRear= thePoint-rearCoord;
// 			 rearID--;
// 			 if(disRear.Norm() <tolerance)
// 			 {
// 				 break;
// 			 }
// 		} //end of for loop:: index
// 
// 		 now we have frontID and rearID in the fiber;
// 
// 		 if(abs(frontID-rearID)<2)
// 			 return 0;
// 		 if(frontID>rearID)
// 			 swap(frontID,rearID);
// 
// 		 vector<int> newFiber;
// 		 for (int index = frontID; index <= rearID; ++index) {
// 			 newFiber.push_back(allFiberElems[index]);
// 		} //end of for loop:: index
// 
// 		 m_data->allFibers[fid]= newFiber;
// 		 return 0;
// 
// 	}
// 
// 	void CFibers::ExtendFibers(CTriSurface& mySurface,float max2Extent){
// 
// 
// 		CJointModel myModel( mySurface, *this);
// 		myModel.FindFiberSurfIntersections();
// 		myModel.FindFibersOnCells();
// 
// 		list<int> fibers2rm;
// 		list<int> fibers2ext;
// 
// 
// 		for (int fiberId = 0; fiberId < m_data->numFibers; ++fiberId) {
// 
// 			CIntersectionType&  frontIntersction= myModel.GetFrontIntersectingCellOfFiber(fiberId);
// 			CIntersectionType&  endIntersction= myModel.GetRearIntersectingCellOfFiber(fiberId);
// 
// 			if ( frontIntersction.IsValid() &&  endIntersction.IsValid())
// 				continue;
// 			if( ! frontIntersction.IsValid()  && ! endIntersction.IsValid() )
// 			{
// 				fibers2rm.push_back(fiberId);
// 				continue;
// 			}
// 			if( frontIntersction.IsValid()  || endIntersction.IsValid() ){
// 				fibers2ext.push_back( fiberId);
// 			}
// 
// 		} //end loop:: fiberId
// 
// 		cout<<"number of fibers to remove: "<< fibers2rm.size()<<endl;
// 		cout<<"number of fibers to extend: "<< fibers2ext.size()<<endl;
// 
// 		BOOST_FOREACH( int& fid2ext, fibers2ext)
// 		{
// 
// 		}
// 
// 	}

// 	void CFibers::RemoveInvalidFibers(CTriSurface& mySurface){
// 
// 
// 		CJointModel myModel( mySurface, *this);
// 		myModel.FindFiberSurfIntersections();
// 		myModel.FindFibersOnCells();
// 
// 		list<int> fibers2remain;
// 		list<int> fibers2ext;
// 
// 
// 		for (int fiberId = 0; fiberId < m_data->numFibers; ++fiberId) {
// 
// 			CIntersectionType&  frontIntersction= myModel.GetFrontIntersectingCellOfFiber(fiberId);
// 			CIntersectionType&  endIntersction= myModel.GetRearIntersectingCellOfFiber(fiberId);
// 
// 			if ( frontIntersction.IsValid() &&  endIntersction.IsValid())
// 				fibers2remain.push_back(fiberId);
// 
// 		} //end loop:: fiberId
// 
// 		cout<<"number of fibers to remain: "<< fibers2remain.size()<<endl;
// 		cout<<"number of fibers to remove: "<< m_data->numFibers-fibers2remain.size()<<endl;
// 
// 		vector<vector<int> > fibersLeft;
// 		BOOST_FOREACH( int& fid, fibers2remain)
// 		{
// 			fibersLeft.push_back( m_data->allFibers[fid]);
// 		}
// 
// 		m_data->allFibers= fibersLeft;
// 		m_data->numFibers=fibersLeft.size();
// 	}

	void CFibers::NeedLineColor(void){
		if(0== m_data->allLineColor.size()){
			m_data->allLineColor.resize(m_data->numFibers);
			m_data->hasLineColor=true;
		}
	}

	void CFibers::SetLineColor(int pid, PointCoordType& theColor){
		if(m_data->hasLineColor== false){
			m_data->allLineColor.resize(m_data->numFibers);
			m_data->hasLineColor=true;
		}
		m_data->allLineColor[pid]= theColor;
	}
	
	float CFibers::GetResolution(void )
	{
		if(0.f== m_data->resolution && m_data->numFibers)
		{
		  vector<int>& sampleFiber = m_data->allFibers[0]; 
		  float allDis=0; 
		  
		  for (int idx = 1; idx < sampleFiber.size() ; ++idx) {
		    VectorType& pointBk= m_data->allPoints[ sampleFiber[idx]];
		    VectorType& pointFw= m_data->allPoints[ sampleFiber[idx-1]]; 
		    VectorType dis= pointBk-pointFw; 
		    allDis+= dis.Norm(); 
		  }//end of for loop::idx
		  
		  m_data->resolution= allDis/(sampleFiber.size()-1);
		}
		
		return m_data->resolution; 
		  
	}

	void CFibers::ResamplePoints(int resampelDegree)
	{
	        cout<<"Resampling points by "<<resampelDegree<<" , now resolution: "<<resampelDegree<<"*"<< this->GetResolution()<<endl;
		assert(resampelDegree>1); 
		vector<VectorType> newAllPoints; 
		vector<vector<int> > newAllLines(m_data->numFibers); 
		vector<PointCoordType> newAllPointsColorsAscii;
		vector<Vector3D<unsigned char> > newAllPointsColorsBinary;
		
		int newPid=0; 
		for (int idLine = 0; idLine < m_data->numFibers ; ++idLine) {
		  int newLineLength=m_data->allFibers[idLine].size()/resampelDegree; 
		  vector<int>& currentLine= m_data->allFibers[idLine]; 
		  for (int newEleIndex = 0; newEleIndex <  newLineLength; ++newEleIndex) {
		    int originalId= currentLine[newEleIndex*resampelDegree]; 
		    newAllPoints.push_back(m_data->allPoints[originalId]);
		    newAllPointsColorsAscii.push_back(m_data->allPointsColorsAscii[originalId]);
		    newAllPointsColorsBinary.push_back(m_data->allPointsColorsBinary[originalId]);
		    newAllLines[idLine].push_back(newPid++);
		  }//end of for loop::newEleIndex
		}//end of for loop::idLine
		
		
		//assignment; 
		m_data->allFibers= newAllLines;
		m_data->allPoints= newAllPoints; 
		m_data->allPointsColorsAscii= newAllPointsColorsAscii; 
		m_data->allPointsColorsBinary= newAllPointsColorsBinary; 
		m_data->numFibers= newAllLines.size(); 
		m_data->numPoints= newAllPoints.size();  
		
		m_data->resolution=0.f; 
		m_data->resolution=this->GetResolution(); 
	}
	
	void CFibers::ResampleLines(int resampelDegree)
	{

		assert(resampelDegree>1); 
		
		vector<VectorType> newAllPoints; 
		vector<PointCoordType> newAllPointsColorsAscii;
		vector<Vector3D<unsigned char> > newAllPointsColorsBinary;
		
		int newNLines= m_data->allFibers.size()/resampelDegree; 
		vector<vector<int> > newAllLines(newNLines); 
	  	cout<<"Resampling lines by "<<resampelDegree<<" ,now has "<<newNLines<<" \\"<<m_data->numFibers<<endl;		
		int newPid=0; 
		for (int idLine = 0; idLine < newNLines ; ++idLine) {
		  int oldLineId= idLine*resampelDegree; 
		  vector<int>& currentLine= m_data->allFibers[oldLineId]; 		  
		  for (int newEleIndex = 0; newEleIndex <  currentLine.size(); ++newEleIndex) {
		    int originalId= currentLine[newEleIndex]; 
		    newAllPoints.push_back(m_data->allPoints[originalId]);
		    newAllPointsColorsAscii.push_back(m_data->allPointsColorsAscii[originalId]);
		    newAllPointsColorsBinary.push_back(m_data->allPointsColorsBinary[originalId]);
		    newAllLines[idLine].push_back(newPid++);
		  }//end of for loop::newEleIndex
		}//end of for loop::idLine
		
		
		//assignment; 
		m_data->allFibers= newAllLines;
		m_data->allPoints= newAllPoints; 
		m_data->allPointsColorsAscii= newAllPointsColorsAscii; 
		m_data->allPointsColorsBinary= newAllPointsColorsBinary; 
		m_data->numFibers= newAllLines.size(); 
		m_data->numPoints= newAllPoints.size(); 
		
	}
	
	void CFibers::RemoveShortFibers(float minFiberLength)
	{
	    assert(minFiberLength>0);
	    int minNPoints= (int) (minFiberLength/m_data->resolution); 
	    vector<VectorType> newAllPoints; 
	    vector<vector<int> > newAllLines; 
	    vector<PointCoordType> newAllPointsColorsAscii;
	    vector<Vector3D<unsigned char> > newAllPointsColorsBinary;
	    
	    int newPid=0; 
	    for (int idLine = 0; idLine < m_data->allFibers.size() ; ++idLine) {
	      vector<int>& currentLine= m_data->allFibers[idLine]; 
	      if(currentLine.size()>minNPoints){
                vector<int> newLine; 
		for (int newEleIndex = 0; newEleIndex <  currentLine.size(); ++newEleIndex) {
		  int originalId= currentLine[newEleIndex]; 
		  newAllPoints.push_back(m_data->allPoints[originalId]);
		  newAllPointsColorsAscii.push_back(m_data->allPointsColorsAscii[originalId]);
		  newAllPointsColorsBinary.push_back(m_data->allPointsColorsBinary[originalId]);
		  newLine.push_back(newPid++);
		}//end of for loop::newEleIndex
		newAllLines.push_back(newLine);
	      }
	    }//end of for loop::idLine
	    cout<<"removing short fibers..."<< newAllLines.size()<<"\\"<<m_data->numFibers<<endl;

	    //assignment; 
	    m_data->allFibers= newAllLines;
	    m_data->allPoints= newAllPoints; 
	    m_data->allPointsColorsAscii= newAllPointsColorsAscii; 
	    m_data->allPointsColorsBinary= newAllPointsColorsBinary; 
	    m_data->numFibers= newAllLines.size(); 
	    m_data->numPoints= newAllPoints.size(); 	    
	}
	
	
	

} //end of KML 
