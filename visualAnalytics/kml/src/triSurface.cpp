/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#include "triSurface.h"
#include <iostream>
#include <fstream>
#include <cmath>
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
#include <algorithm>
#include <string>
#include <cassert>
#include <strstream>
#include <map>
#include <queue>
#include "Vector3D.h"
#include <string.h>
#include <stack>
using namespace std;


namespace KML
{

	class CTriSurfacePrivacy
	{
	public:
		//map<int,PointType> points;

		vector<VectorType> points;
		vector<CellType> cells;

		vector< list<int> > neighborPoints;
		vector< list<int> > neighborCells;

		vector<PointDataType> pointsColor;
		vector<VectorType> pointsNormals;
		vector< float> pointsCurvature;

		vector< int> cellsLabel;
		vector<int> pointsLabel;

		int numOfPoints;
		int numOfCells;
		string inputFileName;

		bool savePointsNormal;
		bool savePointsColor;
		bool savePointsCurvature;
		bool saveCellsLabel;
		bool savePointsLabel;
		bool curvatureThreshSet;

		bool scalarPointDataAdded;
		bool vectorPointDataAdded;
		bool scalarCellDataAdded;
		bool vectorCellDataAdded;
		bool EdgesBuilt; 
		bool NeighborsBuilt; 

		bool writeLines;

		map<string, vector<float> > pointDataScalors; //store all the added scalars;
		map<string, vector< VectorType> > pointDataVectors; // store all the added vectors;

		map<string, vector<float> > cellDataScalors;
		map<string, vector< VectorType> > cellDataVectors;

		vector<CEdgeInfoType > edges;
		vector< set<int> > pointEdges;  //record edges for each point;
		vector< vector<int> > cellEdges; //keep 3 edges for each cell;

		double curvatureThreshold;

                vector<int> pointStatus; 
                vector<int> cellStatus; 
	};

	CTriSurface::~CTriSurface(void)
	{
		delete data;
	}
	CTriSurface::CTriSurface( string filename)
	{

		data=new CTriSurfacePrivacy;

		//set input filename;
		data->inputFileName=string(filename);

		data->savePointsColor=false;
		data->saveCellsLabel=false;
		data->savePointsLabel=false;
		data->savePointsNormal=false;
		data->savePointsCurvature=false;

		data->scalarPointDataAdded=false;
		data->vectorPointDataAdded=false;
		data->scalarCellDataAdded=false;
		data->vectorCellDataAdded=false;
		data->writeLines=false;
		data->EdgesBuilt=false;
		data->NeighborsBuilt=false; 
		
		InitData();
// 		cout<<"Surf points: "<< data->numOfPoints<<endl;
// 		cout<<"Surf cells: "<< data->numOfCells<<endl;
                

	}

	CTriSurface::CTriSurface(CTriSurface& mySurface)
	{
		data=new CTriSurfacePrivacy;
		data->cellDataScalors= mySurface.data-> cellDataScalors;
		data->cellDataVectors= mySurface.data->cellDataVectors;
		data->cellEdges=mySurface.data->cellEdges;
		data->cells=mySurface.data->cells;
		data->cellsLabel=mySurface.data->cellsLabel;
		data->edges=mySurface.data->edges;
		data->inputFileName= mySurface.data->inputFileName;
		data->neighborCells= mySurface.data->neighborCells;
		data->neighborPoints=mySurface.data->neighborPoints;
		data->numOfCells=mySurface.GetNumOfCells();
		data->numOfPoints=mySurface.GetNumOfPoints();
		data->pointDataScalors=mySurface.data->pointDataScalors;
		data->pointDataVectors=mySurface.data->pointDataVectors;
		data->pointEdges=mySurface.data->pointEdges;
		data->points=mySurface.data->points;
		data->pointsColor=mySurface.data->pointsColor;
		data->pointsCurvature=mySurface.data->pointsCurvature;
		data->pointsLabel=mySurface.data->pointsLabel;
		data->pointsNormals=mySurface.data->pointsNormals;
		data->saveCellsLabel=mySurface.data->saveCellsLabel;
		data->savePointsColor=mySurface.data->savePointsColor;
		data->savePointsCurvature=mySurface.data->savePointsCurvature;
		data->savePointsLabel=mySurface.data->savePointsLabel;
		data->savePointsNormal=mySurface.data->savePointsNormal;
		data->scalarCellDataAdded= mySurface.data->scalarCellDataAdded;
		data->scalarPointDataAdded=mySurface.data->scalarPointDataAdded;
		data->vectorCellDataAdded=mySurface.data->vectorCellDataAdded;
		data->vectorPointDataAdded=mySurface.data->vectorPointDataAdded;
		data->writeLines=mySurface.data->writeLines;
		data->EdgesBuilt=mySurface.data->EdgesBuilt;
		data->NeighborsBuilt=mySurface.data->NeighborsBuilt; 
                data->pointStatus=mySurface.data->pointStatus;
                data->cellStatus=mySurface.data->cellStatus; 
 
	}
	void CTriSurface::InitData()
	{
		//////////////////////////////////////////////////////////////////////////
		//read file using vtk;
		vtkPolyDataReader* surfaceDataReader=vtkPolyDataReader::New();
		surfaceDataReader->SetFileName(data->inputFileName.c_str());
		surfaceDataReader->Update();

		vtkPolyData* surfaceData=vtkPolyData::New();
		surfaceData=surfaceDataReader->GetOutput();


		data->numOfPoints=surfaceData->GetNumberOfPoints();
		data->numOfCells=surfaceData->GetNumberOfCells();

		data->points.resize(data->numOfPoints);
		data->cells.resize(data->numOfCells);
		data->neighborPoints.resize(data->numOfPoints);
		data->neighborCells.resize(data->numOfPoints);
		data->pointsColor.resize(data->numOfPoints);
		data->pointsNormals.resize(data->numOfPoints);
		data->pointsCurvature.resize(data->numOfPoints);
		data->pointsLabel.resize(data->numOfPoints);
		data->cellsLabel.resize(data->numOfCells);
		data->cellEdges.resize(data->numOfCells);
                data->pointStatus.resize(data->numOfPoints);
                data->cellStatus.resize(data->numOfCells);
		//////////////////////////////////////////////////////////////////////////
		// build vtkstd::map for points;

		double point[3];
			for(int i=0; i<data->numOfPoints; ++i)
		{
			surfaceData->GetPoint(i,point);
			data->points[i].x=point[0];
			data->points[i].y=point[1];
			data->points[i].z=point[2];
		}




		//////////////////////////////////////////////////////////////////////////
		// build vtkstd::map for cells
		for(int i=0; i<data->numOfCells; ++i)
		{
			vtkIdList* pointList=vtkIdList::New();
			surfaceData->GetCellPoints(i,pointList);

			//each cell has three points, here is the id;
			data->cells[i].x=pointList->GetId(0);
			data->cells[i].y=pointList->GetId(1);
			data->cells[i].z=pointList->GetId(2);
			pointList->Delete();
		}



		//////////////////////////////////////////////////////////////////////////
		// build vtkstd::map for point data;
		// point data should use vectors which have 3 components; since color_scalors seems to meed problems;

		vtkPointData* pointData=vtkPointData::New();
		pointData=surfaceData->GetPointData();

// 		/// for colors
// 		vtkDataArray* vectorArray=NULL;
// 		vectorArray=pointData->GetVectors("Colors");
// 		if (NULL!=vectorArray)
// 		{
// 			data->savePointsColor=true;
// 
// 			for (int i=0; i<data->numOfPoints; ++i)
// 			{
// 				data->pointsColor[i].x=vectorArray->GetComponent(i,0);
// 				data->pointsColor[i].y=vectorArray->GetComponent(i,1);
// 				data->pointsColor[i].z=vectorArray->GetComponent(i,2);
// 			}
// 
// 		}

// 		// for Normals if have
// 		vtkDataArray* normalArray=NULL;
// 		normalArray=pointData->GetVectors("Normals");
// 		if (NULL!=normalArray)
// 		{
// 			data->savePointsNormal=true;
// 
// 			for (int pointId=0; pointId<data->numOfPoints; ++pointId)
// 			{
// 				data->pointsNormals[pointId].x=normalArray->GetComponent(pointId,0);
// 				data->pointsNormals[pointId].y=normalArray->GetComponent(pointId,1);
// 				data->pointsNormals[pointId].z=normalArray->GetComponent(pointId,2);
// 			}
// 		}


/////////////////////////////////////////////////////////////////////////
		//end of init;

		surfaceDataReader->Delete();
	}

	void CTriSurface::SaveAs( string filename)
	{
		fstream outFile;
		outFile.open(filename.c_str(),ios::out);
		if (NULL==outFile)
		{
			cerr<<"Can not create file :"<<filename<<endl;
			exit(EXIT_FAILURE);
		}

		//write header
		outFile<<"# vtk DataFile Version 2.0"<<endl;
		outFile<<"Brain surface "<<endl;
		outFile<<"ASCII"<<endl;
		outFile<<"DATASET POLYDATA"<<endl;


		//write points;
		outFile<<"POINTS "<<data->numOfPoints<<" FLOAT"<<endl;
		for (int i=0; i<data->numOfPoints; ++i)
		{
			outFile<<data->points[i].x<<" "<<data->points[i].y<<" "<<data->points[i].z<<endl;
		}

		//write cells;
		outFile<<"POLYGONS "<<data->numOfCells<<" "<<data->numOfCells*4<<endl;
		for(int i=0; i<data->numOfCells; ++i)
		{
			outFile<<"3 "<<data->cells[i].x<<" "<<data->cells[i].y<<" "<<data->cells[i].z<<endl;
		}

		//////////////////////////////////////////////////////////////////////////
		//write point data if need
		if (data->savePointsLabel|| data->savePointsColor || data->savePointsCurvature || data->savePointsNormal || data->scalarPointDataAdded || data->vectorPointDataAdded){
            outFile << "POINT_DATA " << data->numOfPoints << endl;
        }
        //write color information if need;
        if(data->pointsColor.size() > 0 && data->savePointsColor == true){
            outFile << "COLOR_SCALARS " << "PointColor  3 " << endl;
            for(int i = 0;i < data->numOfPoints;++i){
                outFile << data->pointsColor[i].x << " " << data->pointsColor[i].y << " " << data->pointsColor[i].z << endl;
            }
        }

        //write normal information if need;
        if(data->pointsNormals.size() > 0 && data->savePointsNormal == true){
            outFile << "VECTORS " << "PointNormal " << "FLOAT" << endl;
            for(int i = 0;i < data->numOfPoints;++i){
                outFile << data->pointsNormals[i].x << " " << data->pointsNormals[i].y << " " << data->pointsNormals[i].z << endl;
            }
        }

        //write curvatures information if need
        if(data->pointsCurvature.size() > 0 && data->savePointsCurvature == true){
            outFile << "SCALARS Maximudata->Curvature " << "FLOAT" << endl;
            outFile << "LOOKUP_TABLE default" << endl;
            for(int i = 0;i < data->numOfPoints;++i){
                outFile << data->pointsCurvature[i] << endl;
            }
        }

        //write label information of points if needed
        if(data->pointsLabel.size() > 0 && data->savePointsLabel == true){
            outFile << "SCALARS PointLabel " << "int" << endl;
            outFile << "LOOKUP_TABLE default" << endl;
            for(int i = 0;i < data->numOfPoints;++i){
                outFile << data->pointsLabel[i] << endl;
            }
        }

        //write all added scalars of point if needed;
        if (data->scalarPointDataAdded==true)
		{
			for (vtkstd::map<string, vtkstd::vector<float> >::iterator itAllScalars=data->pointDataScalors.begin(); itAllScalars!=data->pointDataScalors.end(); itAllScalars++)
			{
				outFile<<"SCALARS "<<itAllScalars->first.c_str()<<" float"<<endl;
				outFile<<" LOOKUP_TABLE default"<<endl;
				for (int i=0 ; i<data->numOfPoints; ++i)
				{
					outFile<<itAllScalars->second[i]<<endl;
				}
			}
		}
        //write all added vectors of point if needed;
        if (data->vectorPointDataAdded==true)
		{
			for (vtkstd::map<string, vtkstd::vector< VectorType> >::iterator itAllVectors=data->pointDataVectors.begin(); itAllVectors!=data->pointDataVectors.end(); itAllVectors++)
			{
				outFile<<"VECTORS "<<itAllVectors->first.c_str()<<" FLOAT"<<endl;
				for (int i=0; i<data->numOfPoints; ++i)
				{
					outFile<<itAllVectors->second[i].x<<" "<<itAllVectors->second[i].y<<" "<<itAllVectors->second[i].z<<endl;
				}

			}
		}
        if(data->saveCellsLabel || data->scalarCellDataAdded || data->vectorCellDataAdded){
            outFile << "CELL_DATA " << data->numOfCells << endl;
        }
        //write label information of cells if needed
        if(data->cellsLabel.size() > 0 && data->saveCellsLabel == true){
            outFile << "SCALARS CellLabel " << "int" << " 1" << endl;
            outFile << "LOOKUP_TABLE default" << endl;
            for(int i = 0;i < data->numOfCells;++i){
                outFile << data->cellsLabel[i] << endl;
            }
        }

        //write all added scalars  of cells if needed;
        if (data->scalarCellDataAdded==true)
		{
			for (vtkstd::map<string, vtkstd::vector<float> >::iterator itAllScalars=data->cellDataScalors.begin(); itAllScalars!=data->cellDataScalors.end(); itAllScalars++)
			{
				outFile<<"SCALARS "<<itAllScalars->first.c_str()<<" float"<<endl;
				outFile<<" LOOKUP_TABLE default"<<endl;
				for (int i=0 ; i<data->numOfCells; ++i)
				{
					outFile<<itAllScalars->second[i]<<endl;
				}
			}
		}
        //write all added vectors of cells if needed;
        if (data->vectorCellDataAdded==true)
		{
			for (vtkstd::map<string, vtkstd::vector<VectorType> >::iterator itAllVectors=data->cellDataVectors.begin(); itAllVectors!=data->cellDataVectors.end(); itAllVectors++)
			{
				outFile<<"VECTORS "<<itAllVectors->first.c_str()<<" FLOAT"<<endl;
				for (int i=0; i<data->numOfCells; ++i)
				{
					outFile<<itAllVectors->second[i].x<<" "<<itAllVectors->second[i].y<<" "<<itAllVectors->second[i].z<<endl;
				}

			}
		}
        outFile.close();
    }

    void CTriSurface::Save(void)
    {
        SaveAs(data->inputFileName.c_str());
    }

    const string& CTriSurface::GetFileName(void) const {

    	return data->inputFileName;
    }
    int CTriSurface::GetNumOfCells() const
    {
        return data->numOfCells;
    }

    int CTriSurface::GetNumOfPoints() const
    {
        return data->numOfPoints;
    }

    const list<int> & CTriSurface::GetNeighboringCellsOfPoint(int pointId)
    {
        return data->neighborCells[pointId];
    }

void CTriSurface::GetNeighboringCellsOfPoint( int pointId, set<int>& nbrCells,int range/*=1*/ )
    {
            list<int> nbrPoints;
            this->GetNeighboringPointsOfPoint(pointId,nbrPoints,range);
            this->GetNeighboringCellsOfPoints(nbrPoints,nbrCells);
    }
	
	
void CTriSurface::GetNeighboringCellsOfPoint(int pointId, std::vector< int >& nbrCells, int range)
{
                list<int> nbrPoints;
                this->GetNeighboringPointsOfPoint(pointId,nbrPoints,range);
                this->GetNeighboringCellsOfPoints(nbrPoints,nbrCells);
}


void CTriSurface::GetNeighboringCellsOfPoints(std::list< int >& nbrPoints, std::vector< int >& nbrCells)
{
    vector<unsigned char> cellStatus(data->numOfCells,0);
    for (list<int>::iterator itNbrPoints=nbrPoints.begin(); itNbrPoints!=nbrPoints.end(); itNbrPoints++)
    {
            list<int>& cellList=data->neighborCells[*itNbrPoints];
            for (list<int>::iterator itCells  = cellList.begin();  itCells != cellList.end() ; itCells ++)
            {
              if(false==cellStatus[*itCells])
              {
                cellStatus[*itCells]=true; 
                nbrCells.push_back(*itCells);
              }
            } //end of for loop::                                       
    } //end of for loop::index  
}


    void CTriSurface::GetNeighboringCellsOfPoints(list<int>& nbrPoints, set<int>& nbrCells){

    for (list<int>::iterator itNbrPoints=nbrPoints.begin(); itNbrPoints!=nbrPoints.end(); itNbrPoints++)
    {
            list<int>& cellList=data->neighborCells[*itNbrPoints];
            for (list<int>::iterator itCells  = cellList.begin();  itCells != cellList.end() ; itCells ++)
            {
              nbrCells.insert( *itCells);	
            } //end of for loop::					
    } //end of for loop::index

    }
    
	
    void CTriSurface::GetNeighboringPointsOfPoint(int pointId, list<int> & nbrs, int range)
    {
        if(pointId >= data->numOfPoints)
        {
          cout<<"Invalid pointid :"<<pointId<<endl;
          exit(1);
        }
        if (1==range) //output only the nearest nbr points; and push them back to nbrs;
		{
                  for (list<int>::iterator itNbr=data->neighborPoints[pointId].begin(); itNbr!=data->neighborPoints[pointId].end(); itNbr++)
                  {
                          nbrs.push_back(*itNbr);
                  }

		}
		else if ( range> 1)
		{
                  vector<unsigned char> nbrStatus(data->numOfPoints,0);
                  nbrStatus[pointId]=1; 
                  
                  int currentRange=0;
                  
//                   list<int> newCandidate;
//                   newCandidate.push_back(pointId);
// 
//                   while(currentRange<range)
//                   {
//                     ++currentRange;
//                     list<int> nextNewCandi;
//                     for (list<int>::iterator itCandi = newCandidate.begin(); itCandi != newCandidate.end(); itCandi++) // for each new candidate;
//                     {
//                      // for each nbr of the candi; find nextNewCandi, and push back new find nbrs of pointID;
//                       for ( list<int>::iterator itCandiNbr=data->neighborPoints[*itCandi].begin(); itCandiNbr!=data->neighborPoints[*itCandi].end();itCandiNbr++)
//                       {
//                         if (0==nbrStatus[*itCandiNbr])
//                         {
//                                 nextNewCandi.push_back(*itCandiNbr);
//                                 nbrStatus[*itCandiNbr]=true;
//                                 nbrs.push_back(*itCandiNbr);
//                         }
//                       }
//                     }
// 
//                     newCandidate.clear();
//                     for (list<int>::iterator itNextNewCandi=nextNewCandi.begin(); itNextNewCandi!=nextNewCandi.end(); itNextNewCandi++)
//                     {
//                             newCandidate.push_back(*itNextNewCandi);
//                     }
// 
//                   }
                  
                  queue<int> newCandidate;
                  newCandidate.push(pointId);

                  while(currentRange<range)
                  {
                    ++currentRange;
                    int sizeCurrentRange = newCandidate.size();
                    int cadiProcessed =0; 
                    while ( cadiProcessed++ < sizeCurrentRange)
                    {
                      int theCadiId = newCandidate.front();
                      newCandidate.pop();
                      for ( list<int>::iterator itCandiNbr=data->neighborPoints[theCadiId].begin(); itCandiNbr!=data->neighborPoints[theCadiId].end();itCandiNbr++)
                      {
                        if (0==nbrStatus[*itCandiNbr])
                        {
                                nbrStatus[*itCandiNbr]=true;
                                nbrs.push_back(*itCandiNbr);
                                newCandidate.push(*itCandiNbr);
                        }
                      }
                    }// end while ::current range; 
                 }//end while::currentRange<range; 

		}
		else
		{
                    cerr<<"unaccepted range, it should >=1"<<endl;
                    exit(EXIT_FAILURE);
		}
    }
    
    void CTriSurface::GetNeighboringPoints(list<int>&  srcPoints, list<int>& nbrs, vector<bool>&  indicatorROI){
      
	for(list<int>::iterator it= srcPoints.begin(); it!=srcPoints.end(); it++){
	  indicatorROI[*it]=false; 
	}
 
	for (list<int>::iterator itCandi = srcPoints.begin(); itCandi != srcPoints.end(); itCandi++) // for each new candidate;
	{
		// for each nbr of the candi; find nextNewCandi, and push back new find nbrs of pointID;
		for ( list<int>::iterator itCandiNbr=data->neighborPoints[*itCandi].begin(); itCandiNbr!=data->neighborPoints[*itCandi].end();itCandiNbr++)
		{
			if (true==indicatorROI[*itCandiNbr])
			{
				indicatorROI[*itCandiNbr]=false;
				nbrs.push_back(*itCandiNbr);
			}
		}
	}
      
    }

    void CTriSurface::GetNeighboringPointsOfPoint( int pointId, list<int>& nbrs, vector<float>& roiData, float interestedValue, int range)
    {
        if (1==range) //output only the nearest nbr points; and push them back to nbrs;
		{
			for (list<int>::iterator itNbr=data->neighborPoints[pointId].begin(); itNbr!=data->neighborPoints[pointId].end(); itNbr++)
			{
				if (roiData[*itNbr]==interestedValue)
				{ 
					nbrs.push_back(*itNbr);
				}
			}

		}
		else if ( range> 1)
		{
			vector<bool> nbrStatus(data->numOfPoints,false); // true if it is a nbr with range;
			nbrStatus[pointId]=true; // we want to exclude this point from its nbrs;

			int currentRange=0;
			list<int> newCandidate;
			newCandidate.push_back(pointId);

			while(currentRange<range)
			{
				++currentRange;
				list<int> nextNewCandi;
				for (list<int>::iterator itCandi = newCandidate.begin(); itCandi != newCandidate.end(); itCandi++) // for each new candidate;
				{
					// for each nbr of the candi; find nextNewCandi, and push back new find nbrs of pointID;
					for ( list<int>::iterator itCandiNbr=data->neighborPoints[*itCandi].begin(); itCandiNbr!=data->neighborPoints[*itCandi].end();itCandiNbr++)
					{
						if (false==nbrStatus[*itCandiNbr] && roiData[*itCandiNbr]==interestedValue )
						{
							nextNewCandi.push_back(*itCandiNbr);
							nbrStatus[*itCandiNbr]=true;
							nbrs.push_back(*itCandiNbr);
						}
					}
				}

				newCandidate.clear();
				for (list<int>::iterator itNextNewCandi=nextNewCandi.begin(); itNextNewCandi!=nextNewCandi.end(); itNextNewCandi++)
				{
					newCandidate.push_back(*itNextNewCandi);
				}

			}

		}
		else
		{
			cerr<<"unaccepted range, it should >=1"<<endl;
			exit(EXIT_FAILURE);
		}
    }

    void CTriSurface::GetNeighboringPointsOfPoint(int pointId, list<int> & nbrs, vector<int> & roiData, int interestedValue, int range)
    {
        if (1==range) //output only the nearest nbr points; and push them back to nbrs;
		{
			for (list<int>::iterator itNbr=data->neighborPoints[pointId].begin(); itNbr!=data->neighborPoints[pointId].end(); itNbr++)
			{
				if (roiData[*itNbr]==interestedValue)
				{
					nbrs.push_back(*itNbr);
				}
			}

		}
		else if ( range> 1)
		{
			vector<bool> nbrStatus(data->numOfPoints,false); // true if it is a nbr with range;
			nbrStatus[pointId]=true; // we want to exclude this point from its nbrs;

			int currentRange=0;
			list<int> newCandidate;
			newCandidate.push_back(pointId);

			while(currentRange<range)
			{
				++currentRange;
				list<int> nextNewCandi;
				for (list<int>::iterator itCandi = newCandidate.begin(); itCandi != newCandidate.end(); itCandi++) // for each new candidate;
				{
					// for each nbr of the candi; find nextNewCandi, and push back new find nbrs of pointID;
					for ( list<int>::iterator itCandiNbr=data->neighborPoints[*itCandi].begin(); itCandiNbr!=data->neighborPoints[*itCandi].end();itCandiNbr++)
					{
						if (false==nbrStatus[*itCandiNbr] && roiData[*itCandiNbr]==interestedValue )
						{
							nextNewCandi.push_back(*itCandiNbr);
							nbrStatus[*itCandiNbr]=true;
							nbrs.push_back(*itCandiNbr);
						}
					}
				}

				newCandidate.clear();
				for (list<int>::iterator itNextNewCandi=nextNewCandi.begin(); itNextNewCandi!=nextNewCandi.end(); itNextNewCandi++)
				{
					newCandidate.push_back(*itNextNewCandi);
				}

			}

		}
		else
		{
			cerr<<"unaccepted range, it should >=1"<<endl;
			exit(EXIT_FAILURE);
		}
    }

    void CTriSurface::GetNbrsOfPointOnRange(int pointId, list<int> & outerNbrs, int range)
    {
        if (0>=range)
		{
			cerr<<"unaccepted range, it should >0"<<endl;
			exit(EXIT_FAILURE);

		}
        vector<bool> nbrStatus(data->numOfPoints, false); // true if it is a nbr  at and within range;
        nbrStatus[pointId] = true; // we want to exclude this point from its nbrs;
        int currentRange = 0;
        list<int> newCandidate;
        newCandidate.push_back(pointId);
        while(currentRange < range){
            ++currentRange;
            list<int> nextNewCandi;
            for(list<int>::iterator itCandi = newCandidate.begin();itCandi != newCandidate.end();itCandi++){
                // for each nbr of the candi; find nextNewCandi, and push back new find nbrs of pointID;
                for(list<int>::iterator itCandiNbr = data->neighborPoints[*itCandi].begin();itCandiNbr != data->neighborPoints[*itCandi].end();itCandiNbr++){
                    if(false == nbrStatus[*itCandiNbr]){
                        nextNewCandi.push_back(*itCandiNbr);
                        nbrStatus[*itCandiNbr] = true;
                    }
                }

            }

            newCandidate.clear();
            for(list<int>::iterator itNextNewCandi = nextNewCandi.begin();itNextNewCandi != nextNewCandi.end();itNextNewCandi++){
                newCandidate.push_back(*itNextNewCandi);
            }
            //if =, then the newCandidates are outer nbrs;
            if(currentRange == range){
                outerNbrs = newCandidate;
            }
        }

    }

    void CTriSurface::Propagate2All(int pointId, list<int>& groupMembers, map<int,list<int> >& nbrLabels,vector<int>& roiData, int maxSize)
	{
			if (-1==maxSize)
			{
				maxSize=data->numOfPoints+1;
			}
  			vector<bool> isVisited(roiData.size(),false);
			int theLabel=roiData[pointId];
			groupMembers.push_back(pointId);
			isVisited[pointId]=true;
			queue<int> toDoQueue;
			toDoQueue.push(pointId);

			
			//while todo queue not empty && groupmember.size <maxSize
			while(toDoQueue.size() && groupMembers.size()< maxSize)
			{
				int currentQueueSize=toDoQueue.size();
				for (int toDoIndex=0; toDoIndex< currentQueueSize; ++toDoIndex)
				{
					int roiPoint=toDoQueue.front();
					toDoQueue.pop(); // delete the point; 
					list<int>& nbrRoiPoint=data->neighborPoints[roiPoint];

					for (list<int>::iterator itNbrRoiPoint= nbrRoiPoint.begin(); itNbrRoiPoint!= nbrRoiPoint.end(); itNbrRoiPoint++)
					{
						if (! isVisited[*itNbrRoiPoint])
						{
							isVisited[*itNbrRoiPoint]=true;
							if (theLabel==roiData[*itNbrRoiPoint])
							{
								toDoQueue.push(*itNbrRoiPoint);
								groupMembers.push_back(*itNbrRoiPoint);
							}
							else
							{
								nbrLabels[roiData[*itNbrRoiPoint]].push_back(*itNbrRoiPoint);
							}
						}
					}
					
				}
			}
			//do 
			// for each point in the list,
			// get its nbring points;
			// pop the point out; 
			// for each of the nbring points, if not visited{
			// mark as visited;  if interested, putback into the todo queue and putback into groupmembers; 
			// if not interested, put it into nbrLabesl;
			// if visited, continue; 	
    }

    void CTriSurface::Propagate2All(int pointId, list<int>& groupMembers, map<float, list<int> >& nbrLabels,vector<float>& roiData, int maxSize)
	{
		if (-1==maxSize)
		{
			maxSize=data->numOfPoints+1;
		}
		vector<bool> isVisited(roiData.size(),false);
		float theLabel=roiData[pointId];
		groupMembers.push_back(pointId);
		isVisited[pointId]=true;
		queue<int> toDoQueue;
		toDoQueue.push(pointId);


		//while todo queue not empty && groupmember.size <maxSize
		while(toDoQueue.size() && groupMembers.size()< maxSize)
		{
			int currentQueueSize=toDoQueue.size();
			for (int toDoIndex=0; toDoIndex< currentQueueSize; ++toDoIndex)
			{
				int roiPoint=toDoQueue.front();
				toDoQueue.pop(); // delete the point; 
				list<int>& nbrRoiPoint=data->neighborPoints[roiPoint];

				for (list<int>::iterator itNbrRoiPoint= nbrRoiPoint.begin(); itNbrRoiPoint!= nbrRoiPoint.end(); itNbrRoiPoint++)
				{
					if (! isVisited[*itNbrRoiPoint])
					{
						isVisited[*itNbrRoiPoint]=true;
						if (theLabel==roiData[*itNbrRoiPoint])
						{
							toDoQueue.push(*itNbrRoiPoint);
							groupMembers.push_back(*itNbrRoiPoint);
						}
						else
						{
							nbrLabels[roiData[*itNbrRoiPoint]].push_back(*itNbrRoiPoint);
						}
					}
				}

			}
		}

    }


    const VectorType & CTriSurface::GetPointNormal(int pointId) const
    {
        return data->pointsNormals[pointId];
    }

    float CTriSurface::GetPointCurvature(int pointId)
    {
        return data->pointsCurvature[pointId];
    }

    const VectorType & CTriSurface::GetPointCoords(int pointId)
    {
        return data->points[pointId];
    }

    const VectorType & CTriSurface::GetPointColor(int pointId) const
    {
        return data->pointsColor[pointId];
    }

    void CTriSurface::SetPointColor(int pointId, const VectorType & floatColor)
    {
        data->pointsColor[pointId] = floatColor;
    }

    int CTriSurface::GetPointLabel(int pointId)
    {
        return data->pointsLabel[pointId];
    }

   vector<int>& CTriSurface::GetAllPointLabel(){
	   return data->pointsLabel;
   }
   vector<int>& CTriSurface::GetAllCellLabel(){
	   return data->cellsLabel;
   }
    void CTriSurface::SetPointLabel(int pointId, int label)
    {
        data->pointsLabel[pointId] = label;
    }

    void CTriSurface::GetPointsWithLabel(int label, list<int> & targetPoints)
    {
        if(0==targetPoints.size()){
    		for (int pointIndex = 0; pointIndex < data->numOfPoints; ++pointIndex) {
				if(label==data->pointsLabel[pointIndex]){
					targetPoints.push_back(pointIndex);
				}
			}
    	}else{
    		cerr<<"Error: target not empty!"<<endl;
    		exit(EXIT_FAILURE);
    	}
    }

    int CTriSurface::GetCellLabel(int cellId)
    {
        return data->cellsLabel[cellId];
    }

    void CTriSurface::SetCellLabel(int cellId, int label)
    {
        data->cellsLabel[cellId] = label;
    }

    void CTriSurface::GetCellsWithLabel(int label, list<int> & targetCells)
    {
        if(0==targetCells.size()){
    		for (int cellIndex = 0; cellIndex < data->numOfCells; ++cellIndex) {
				if(label==data->cellsLabel[cellIndex]){
					targetCells.push_back(cellIndex);
				}
			}
    	}else{
    		cerr<<"Error: target not empty"<<endl;
    		exit(EXIT_FAILURE);
    	}
    }

    VectorType CTriSurface::GetCellNormal(int cellId)
    {
        VectorType edge1;
        VectorType edge2;
        edge1 = data->points[data->cells[cellId].x] - data->points[data->cells[cellId].y];
        edge2 = data->points[data->cells[cellId].x] - data->points[data->cells[cellId].z];
        VectorType norm;
        norm = CrossProduct(edge1, edge2);
        norm.Normalize();
        return norm;
    }

    float CTriSurface::GetCellArea(int cellId)
    {
        VectorType edge1;
        VectorType edge2;
        int pointId1 = data->cells[cellId].x;
        int pointId2 = data->cells[cellId].y;
        int pointId3 = data->cells[cellId].z;
        edge1 = data->points[pointId2] - data->points[pointId1];
        edge2 = data->points[pointId3] - data->points[pointId1];
        float cosAngle = CosAngleBetweenVectors(edge1, edge2);
        float sinAngle = sqrt(1 - cosAngle * cosAngle);
        float area = 0.5 * edge1.Norm() * edge2.Norm() * sinAngle;
        return area;
    }

    float CTriSurface::GetCellsTotalArea(const CellIDsType & cellIds)
    {
        float area = 0;
        for(CellIDsType::const_iterator itCellIDs = cellIds.begin();itCellIDs != cellIds.end();itCellIDs++){
            area += this->GetCellArea(*itCellIDs);
        }
        if(area == 0){
            cerr << "Warning, area is 0 !!!!" << endl;
        }
        return area;
    }

    float CTriSurface::GetCellsTotalArea(int startId, int numOfCells)
    {
        float area = 0;
        for(int i = startId;i < numOfCells;++i){
            area += this->GetCellArea(i);
        }
        if(area == 0){
            cerr << "Warning, area is 0 !!!!" << endl;
        }
        return area;
    }

    const CellType & CTriSurface::GetCell(int cellId)
    {
        return data->cells[cellId];
    }

    void CTriSurface::NeedCellsLabel(bool bLabel)
    {
        data->saveCellsLabel = bLabel;
    }

    void CTriSurface::NeedPointsNormal(bool bNormal)
    {
        data->savePointsNormal = bNormal;
        this->GeneratePointNormal();
    }

    void CTriSurface::NeedPointsCurvature(bool bCurvature)
    {
        data->savePointsCurvature = bCurvature;
        this->GeneratePointCurvature();
    }

    void CTriSurface::NeedPointsColor(bool bColor)
    {
        data->savePointsColor = bColor;
    }

    void CTriSurface::NeedPointsLabel(bool bLabel)
    {
        data->savePointsLabel = bLabel;
    }

    float CTriSurface::GetPointComponent(int pointId, int component)
    {
        if(pointId < 0 || pointId > data->numOfPoints){
            cerr << "Bad point id, valid range 0~" << data->numOfPoints << endl;
            return -1;
        }
        if(0 == component)
            return data->points[pointId].x;

        else
            if(1 == component){
                return data->points[pointId].y;
            }else
                if(2 == component){
                    return data->points[pointId].z;
                }else{
                    cerr << "Bad component number, valid id is 0~2" << endl;
                    return -1;
                }


    }

    void CTriSurface::SetPointComponent(int pointId, int component, float componentValue)
    {
        if (pointId<0 || pointId>data->numOfPoints)
		{
			cerr<<"Bad point id, valid range 0~"<<data->numOfPoints<<endl;
			exit(EXIT_FAILURE);
		}
        if (0==component)
		{
			data->points[pointId].x=componentValue;
		}
		else if (1==component)
		{
			 data->points[pointId].y=componentValue;
		}
		else if (2==component)
		{
			 data->points[pointId].z=componentValue;
		}
		else
		{
			cerr<<"Bad component number, valid id is 0~2"<<endl;
			exit(EXIT_FAILURE);
		}
    }

    void CTriSurface::BuildNeighbors(void)
    {

      if(! data->NeighborsBuilt)
      {
// 	cout << "build neighbors...";
        //build neighbors;
        int cellIndex = 0;
        for(vector<CellType>::iterator itCells = data->cells.begin();itCells != data->cells.end();itCells++){
            //build neighboring points for each point;
            if(find(data->neighborPoints[itCells->x].begin(), data->neighborPoints[itCells->x].end(), itCells->y) == data->neighborPoints[itCells->x].end()){
                data->neighborPoints[itCells->x].push_back(itCells->y);
                data->neighborPoints[itCells->y].push_back(itCells->x);
            }
            if(find(data->neighborPoints[itCells->x].begin(), data->neighborPoints[itCells->x].end(), itCells->z) == data->neighborPoints[itCells->x].end()){
                data->neighborPoints[itCells->x].push_back(itCells->z);
                data->neighborPoints[itCells->z].push_back(itCells->x);
            }
            if(find(data->neighborPoints[itCells->y].begin(), data->neighborPoints[itCells->y].end(), itCells->z) == data->neighborPoints[itCells->y].end()){
                data->neighborPoints[itCells->y].push_back(itCells->z);
                data->neighborPoints[itCells->z].push_back(itCells->y);
            }
            //build neighboring cells for each point;
            if(find(data->neighborCells[itCells->x].begin(), data->neighborCells[itCells->x].end(), cellIndex) == data->neighborCells[itCells->x].end()){
                data->neighborCells[itCells->x].push_back(cellIndex);
            }
            if(find(data->neighborCells[itCells->y].begin(), data->neighborCells[itCells->y].end(), cellIndex) == data->neighborCells[itCells->y].end()){
                data->neighborCells[itCells->y].push_back(cellIndex);
            }
            if(find(data->neighborCells[itCells->z].begin(), data->neighborCells[itCells->z].end(), cellIndex) == data->neighborCells[itCells->z].end()){
                data->neighborCells[itCells->z].push_back(cellIndex);
            }
            ++cellIndex;
        };
	
	data->NeighborsBuilt=true; 
// 	cout<<"...done."<<endl;
	}
	else
	{
	  cout<<"neighbors have been built."<<endl; 
	}
	
    }

    void CTriSurface::GeneratePointNormal(void)
    {
        //////////////////////////////////////////////////////////////////////////
        //read file using vtk;
        vtkPolyDataReader *surfaceDataReader = vtkPolyDataReader::New();
        surfaceDataReader->SetFileName(data->inputFileName.c_str());
        surfaceDataReader->Update();
        //////////////////////////////////////////////////////////////////////////
        // generate norms
        vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
        normals->SetInput(surfaceDataReader->GetOutput());
        normals->SetNonManifoldTraversal(true);
        normals->SetSplitting(false);
        normals->SetComputePointNormals(true);
        normals->SetFeatureAngle(30);
        normals->Update(); //need this call
        vtkPolyData *afterNorm = vtkPolyData::New();
        afterNorm = normals->GetOutput();
        vtkPointData *normData = vtkPointData::New();
        normData = afterNorm->GetPointData();
        normData->Update();
        vtkDataArray* normArray=NULL;
        normArray = normData->GetNormals("Normals");
        if (NULL!=normArray)
		{
			for (int i=0; i<data->numOfPoints; ++i)
			{
				data->pointsNormals[i].x=normArray->GetComponent(i,0);
				data->pointsNormals[i].y=normArray->GetComponent(i,1);
				data->pointsNormals[i].z=normArray->GetComponent(i,2);
			}

		}
        surfaceDataReader->Delete();
        //end of norm generation;
    }

    void CTriSurface::BuildEdge(void)
    {       
      
      if( !data->EdgesBuilt)
      {
	
	cout<<"build edges...";
	data->EdgesBuilt=true; 
	
        data->pointEdges.resize(data->numOfPoints);
        data->cellEdges.resize(data->numOfCells);
        map<CEdgeInfoType,int> edgeVisitStatus;
        map<CEdgeInfoType,int> edgeId;
        int id = 0;
        for(int cellIndex = 0;cellIndex < this->GetNumOfCells();++cellIndex){
            vector<CEdgeInfoType> allEdges(3);
            int vertex1 = data->cells[cellIndex].x;
            int vertex2 = data->cells[cellIndex].y;
            int vertex3 = data->cells[cellIndex].z;
            CEdgeInfoType edge1(vertex1, vertex2);
            CEdgeInfoType edge2(vertex1, vertex3);
            CEdgeInfoType edge3(vertex2, vertex3);
            edge1.Sort();
            edge2.Sort();
            edge3.Sort();
            allEdges[0] = edge1;
            allEdges[1] = edge2;
            allEdges[2] = edge3;
            for(int edgeIndex = 0;edgeIndex < 3;++edgeIndex){
                if(0 == edgeVisitStatus[allEdges[edgeIndex]]){
                    edgeId[allEdges[edgeIndex]] = id;
                    edgeVisitStatus[allEdges[edgeIndex]] = 1;
                    //record edges;
                    data->edges.push_back(allEdges[edgeIndex]);
                    data->edges[id].cellId1 = cellIndex;
                    //record in cell;
                    data->cellEdges[cellIndex].push_back(id);
                    //record in points;
                    data->pointEdges[allEdges[edgeIndex].pointId1].insert(id);
                    data->pointEdges[allEdges[edgeIndex].pointId2].insert(id);
                    ++id;
                }else
                    if(1 == edgeVisitStatus[allEdges[edgeIndex]]){
                        //only need to update edge information ;
                        int oriEdgeId = edgeId[allEdges[edgeIndex]];
			data->edges[oriEdgeId].cellId2 = cellIndex;
                        data->cellEdges[cellIndex].push_back(oriEdgeId);
                        edgeVisitStatus[allEdges[edgeIndex]] = 2;
                    }else{
                        // no opeartion;
                    }

            } //end of for loop

        } //end of for loop

      cout << "done." << endl;
      }
      else
      {
	cout<<"edges have been built."<<endl;	

      }
    }

    int CTriSurface::GetNumOfEdges(void) const
    {
        return data->edges.size();
    }

    CEdgeInfoType & CTriSurface::GetEdge(int lineId)
    {
        return data->edges[lineId];
    }

    set<int> & CTriSurface::GetPointEdges(int pointId)
    {
        return data->pointEdges[pointId];
    }

    vector<int> & CTriSurface::GetCellEdges(int cellId)
    {
        return data->cellEdges[cellId];
    }

    void CTriSurface::WarpPointLabel2CellLabel(void)
    {
    	for (int cellIndex=0; cellIndex<data->numOfCells; ++cellIndex)
    	{
				CellType& currentCell=data->cells[cellIndex];

       			if (data->pointsLabel[currentCell.x]== data->pointsLabel[currentCell.y] || data->pointsLabel[currentCell.x]== data->pointsLabel[currentCell.z]  )	{
       				data->cellsLabel[cellIndex]=data->pointsLabel[currentCell.x];
    			}
    			else	{
    				data->cellsLabel[cellIndex]=data->pointsLabel[currentCell.y];
    			}
		}

    	data->saveCellsLabel=true;
    }

    void CTriSurface::GeneratePointCurvature(void)
    {
        //////////////////////////////////////////////////////////////////////////
        //read file using vtk;
        vtkPolyDataReader *surfaceDataReader = vtkPolyDataReader::New();
        surfaceDataReader->SetFileName(data->inputFileName.c_str());
        surfaceDataReader->Update();
        //vtkPolyData* surfaceData=vtkPolyData::New();
        //surfaceData=surfaceDataReader->GetOutput();
        //////////////////////////////////////////////////////////////////////////
        //generate curvature;
        vtkCurvatures *curv = vtkCurvatures::New();
        curv->SetInput(surfaceDataReader->GetOutput());
        curv->SetCurvatureTypeToMean();
        //curv->SetCurvatureTypeToMaximum();
        curv->Update();
        vtkPolyData *afterCurv = vtkPolyData::New();
        afterCurv = curv->GetOutput();
        afterCurv->Update();
        vtkPointData *curvData = vtkPointData::New();
        curvData = afterCurv->GetPointData();
        vtkDataArray* curvArray=NULL;
        curvArray = curvData->GetScalars("Mean_Curvature");
        if(NULL!= curvArray)
		{

			for (int i=0; i<data->numOfPoints; ++i)
			{
				data->pointsCurvature[i]=curvArray->GetComponent(i,0);
			}
			curvArray->Delete();
		}
        surfaceDataReader->Delete();
    }

    void CTriSurface::AddPointDataScalar(vector<float> & src, string dataName)
    {
        if (src.size()==data->numOfPoints)
		{
			data->scalarPointDataAdded=true;
			data->pointDataScalors[dataName]=src;
		}
		else
		{
			cerr<<"Number of elements are not consist"<<endl;
			cerr<<"Point number should be :"<<data->numOfPoints<<endl;
			exit(EXIT_FAILURE);
		}
    }

    void CTriSurface::AddPointDataVector(vector<VectorType> & src, string vectorName)
    {
        if (src.size()==data->numOfPoints)
		{
			data->vectorPointDataAdded=true;
			data->pointDataVectors[vectorName]=src;
		}
		else
		{
			cerr<<"Number of elements are not consist"<<endl;
			cerr<<"Point number should be :"<<data->numOfPoints<<endl;
			exit(EXIT_FAILURE);
		}
    }

    void CTriSurface::AddCellDataScalar(vector<float> & src, string dataName)
    {
        data->cellDataScalors[dataName] = src;
        data->scalarCellDataAdded = true;
    }

    void CTriSurface::AddCellDataVector(vector<VectorType> & src, string vectorName)
    {
        data->cellDataVectors[vectorName] = src;
        data->vectorCellDataAdded = true;
    }

    void CTriSurface::ReadPointDataScalar(string dataName, int skipLines)
    {
        data->scalarPointDataAdded = true;
        fstream oriSrcFile;
        oriSrcFile.open(data->inputFileName.c_str(), ios::in);
        if (NULL==oriSrcFile)
		{
			cerr<<"can not open file "<<data->inputFileName.c_str()<<endl;
			exit(EXIT_FAILURE);
		}
        // skip the lines of vertex and faces;
        if(-1 == skipLines){
            skipLines = data->numOfCells + data->numOfPoints;
        }
        string line;
        while(skipLines--){
            getline(oriSrcFile, line);
        }
        // read each line and try finding the matching property name;
        do{
            getline(oriSrcFile, line);
        } while(line.find(dataName) == string::npos && !oriSrcFile.eof());
        getline(oriSrcFile, line); // lookup tabel line;
        vector<float> scalarData(data->numOfPoints, 0);
        for(int pointIndex = 0;pointIndex < data->numOfPoints;pointIndex++){
            oriSrcFile >> scalarData[pointIndex];
        }
        oriSrcFile.close();
        data->pointDataScalors[dataName] = scalarData;
    }

    void CTriSurface::ReadCellDataScalar(string dataName, int skipLines)
    {
        data->scalarCellDataAdded = true;
        fstream oriSrcFile;
        oriSrcFile.open(data->inputFileName.c_str(), ios::in);
        if (NULL==oriSrcFile)
		{
			cerr<<"can not open file "<<data->inputFileName.c_str()<<endl;
			exit(EXIT_FAILURE);
		}
        // skip the lines of vertex and faces;
        if(-1 == skipLines){
            skipLines = data->numOfCells + data->numOfPoints;
        }
        string line;
        while(skipLines--){
            getline(oriSrcFile, line);
        }
        // read each line and try finding the matching property name;
        do{
            getline(oriSrcFile, line);
        } while(line.find(dataName) == string::npos && !oriSrcFile.eof());
        getline(oriSrcFile, line); // lookup tabel line;
        vector<float> scalarData(data->numOfCells, 0);
        for(int cellIndex = 0;cellIndex < data->numOfCells;cellIndex++){
            oriSrcFile >> scalarData[cellIndex];
        }
        oriSrcFile.close();
        data->cellDataScalors[dataName] = scalarData;
    }

    void CTriSurface::ReadPointDataVector(string vectorName, int skipLines)
    {
        data->vectorPointDataAdded = true;
        fstream oriSrcFile;
        oriSrcFile.open(data->inputFileName.c_str(), ios::in);
        if (NULL==oriSrcFile)
		{
			cerr<<"can not open file "<<data->inputFileName.c_str()<<endl;
			exit(EXIT_FAILURE);
		}
        // skip the lines of vertex and faces;
        if(-1 == skipLines){
            skipLines = data->numOfCells + data->numOfPoints;
        }
        string line;
        while(skipLines--){
            getline(oriSrcFile, line);
        }
        // read each line and try finding the matching property name;
        do{
            getline(oriSrcFile, line);
        } while(line.find(vectorName) == string::npos && !oriSrcFile.eof());
        /*getline(oriSrcFile,line); // lookup table line; */
        // no lookup table for vector data;
        vector<VectorType> vectorData(data->numOfPoints);
        for(int pointIndex = 0;pointIndex < data->numOfPoints;pointIndex++){
            oriSrcFile >> vectorData[pointIndex].x;
            oriSrcFile >> vectorData[pointIndex].y;
            oriSrcFile >> vectorData[pointIndex].z;
        }
        oriSrcFile.close();
        data->pointDataVectors[vectorName] = vectorData;
    }

    void CTriSurface::ReadCellDataVector(string vectorName, int skipLines)
	{
		data->vectorCellDataAdded=true;
		fstream oriSrcFile;
		oriSrcFile.open(data->inputFileName.c_str(),ios::in);
		if (NULL==oriSrcFile)
		{
			cerr<<"can not open file "<<data->inputFileName.c_str()<<endl;
			exit(EXIT_FAILURE);
		}

		// skip the lines of vertex and faces;
		if (-1==skipLines)
		{
			skipLines= data->numOfCells+data->numOfPoints;
		}

		string line;
		while(skipLines--)
		{
			getline(oriSrcFile,line);
		}

		// read each line and try finding the matching property name;
		do
		{
			getline(oriSrcFile,line);
		} while (line.find(vectorName)==string::npos && !oriSrcFile.eof());

		/*getline(oriSrcFile,line); // lookup table line; */ // no lookup table for vector data;


		vector<VectorType> vectorData(data->numOfCells);
		for (int cellIndex=0; cellIndex<data->numOfCells; cellIndex++)
		{
			oriSrcFile>>vectorData[cellIndex].x;
			oriSrcFile>>vectorData[cellIndex].y;
			oriSrcFile>>vectorData[cellIndex].z;
		}

		oriSrcFile.close();
		data->cellDataVectors[vectorName]=vectorData;
	}


    vector<int>& CTriSurface::GetPointDataLabel(void)
    {
    	return data->pointsLabel;
    }

    vector<int>& CTriSurface::GetCellDataLable(void)
    {
    	return data->cellsLabel;
    }

    vector<VectorType>& CTriSurface::GetPointDataColor(void)
    {
    	return data->pointsColor;
    }

	vector<float>& CTriSurface::GetPointDataScalar( string dataName )
	{
		if (data->pointDataScalors[dataName].size())
		{
			return data->pointDataScalors[dataName];
		}
		else
		{
			cerr<<"No such data"<<dataName.c_str()<<endl;
			exit(EXIT_FAILURE);
		}
	}

	vector<float>& CTriSurface::GetCellDataScalar( string dataName )
	{
		if (data->cellDataScalors[dataName].size())
		{
			return data->cellDataScalors[dataName];
		}
		else
		{

			cerr<<"No such data"<<dataName.c_str()<<endl;
			exit(EXIT_FAILURE);
		}

	}

	vector<VectorType>& CTriSurface::GetPointDataVector( string dataName )
	{
		if (data->pointDataVectors[dataName].size())
		{
			return data->pointDataVectors[dataName];
		}
		else
		{
			cerr<<"No such data"<<dataName.c_str()<<endl;
			exit(EXIT_FAILURE);
		}

	}

	vector<VectorType>& CTriSurface::GetCellDataVector( string dataName )
	{
		if (data->cellDataVectors[dataName].size())
		{
			return data->cellDataVectors[dataName];
		}
		else
		{
			cerr<<"No such data"<<dataName.c_str()<<endl;
			exit(EXIT_FAILURE);
		}
	}

	void CTriSurface::WriteLineBorders( vector<bool>& tag, string fileName )
	{
		fstream polydataLineStrm;
		polydataLineStrm.open(fileName.c_str(),ios::out);
		assert(polydataLineStrm);

		//write header
		polydataLineStrm<<"# vtk DataFile Version 2.0"<<endl;
		polydataLineStrm<<"Brain surface "<<endl;
		polydataLineStrm<<"ASCII"<<endl;
		polydataLineStrm<<"DATASET POLYDATA"<<endl;


		//write points;
		polydataLineStrm<<"POINTS "<<data->numOfPoints<<" FLOAT"<<endl;
		for (int i=0; i<data->numOfPoints; ++i)
		{
			polydataLineStrm<<data->points[i].x<<" "<<data->points[i].y<<" "<<data->points[i].z<<endl;
		}

		//write lines;
		int numOfLines=0;
		for(int tagIndex=0; tagIndex<tag.size(); ++tagIndex)
		{

				if (tag[tagIndex])
				{
					++numOfLines;
				}
		}

		polydataLineStrm<<"LINES "<<numOfLines<<" "<<numOfLines*3<<endl;

		for(int tagIndex=0; tagIndex<tag.size(); ++tagIndex)
		{

			if (tag[tagIndex])
			{
				polydataLineStrm<<"2 "<<data->edges[tagIndex].pointId1<<" "<<data->edges[tagIndex].pointId2<<endl;
			}
		}

		polydataLineStrm.close();

	}

	vector<VectorType>& CTriSurface::GetAllPointColor( void )
	{
		return data->pointsColor;
	}

	void CTriSurface::SmoothCellLabelBySimpleMajority(int itN)
	{
	    this->BuildEdge();
	    vector<CEdgeInfoType>& allEdges= data->edges; 
   
	    while (itN--)
	    {
	      vector<int> oldLabel= data->cellsLabel;
	      cout<<"SmoothCellLabel Iteration: "<<itN<<endl;
	      for (int cellIdx = 0; cellIdx < this->GetNumOfCells(); ++cellIdx) {
		 if(cellIdx%15000==0)
		   cout<<cellIdx<<endl;
		 int edge1= data->cellEdges[cellIdx][0];
		 int edge2= data->cellEdges[cellIdx][1];
		 int edge3= data->cellEdges[cellIdx][2];
		 int cell1=allEdges[edge1].cellId1==cellIdx ? allEdges[edge1].cellId2 :allEdges[edge1].cellId1; 
		 int cell2=allEdges[edge2].cellId1==cellIdx ? allEdges[edge2].cellId2 :allEdges[edge2].cellId1; 
		 int cell3=allEdges[edge3].cellId1==cellIdx ? allEdges[edge3].cellId2 :allEdges[edge3].cellId1; 
		 if( oldLabel[cell1] == oldLabel[cell2] || oldLabel[cell1] ==oldLabel[cell3])
		   data->cellsLabel[cellIdx] = oldLabel[cell1];
		 if( oldLabel[cell2] == oldLabel[cell3] )
		   data->cellsLabel[cellIdx] = oldLabel[cell2];
	      } //end for loop::cellIdx      
	      
	    } // end of while      
	}

	void CTriSurface::SmoothPointLabel( int range, int numBins,  int itN )
	{
		while (itN--)
		{
			cout<<"SmoothPointLabel Iteration: "<<itN<<endl;
			vector<int> oldLabel=data->pointsLabel;

			for(int pointId=0; pointId<data->numOfPoints; ++pointId)
			{
				if (pointId%15000==0)
				{
					cout<<pointId<<endl;
				}
				list<int> nbrs;
				this->GetNbrsOfPointOnRange(pointId,nbrs,range);

				//voting; 
				vector<int> hist(numBins,0);
				for(list<int>::iterator itNbrs= nbrs.begin(); itNbrs!=nbrs.end(); itNbrs++)
				{
					++hist[oldLabel[*itNbrs]];
				}

				int maxLabel=0; 
				int maxValue=-1;
				findMaxValueAndIndex(hist,maxValue,maxLabel);
				data->pointsLabel[pointId]=maxLabel;

			}
		} // end of while

	} 

	void CTriSurface::RemoveSmallIsolatedLabeledRegion( int sizeThresh/*=30 */ )
	{
		vector<bool> isVisited( this->GetNumOfPoints(),false);
		vector<int>& pointData= this->GetAllPointLabel();

		for (int pointIndex = 0; pointIndex <  this->GetNumOfPoints() ; ++pointIndex)
		{
			if (false==isVisited[pointIndex])
			{
				map<int, list<int> > nbrLabels;
				list<int> groupMembers;
				this->Propagate2All(pointIndex,groupMembers,nbrLabels,pointData);


				for(list<int>::iterator itGroupMembers= groupMembers.begin(); itGroupMembers!=groupMembers.end(); itGroupMembers++)
				{
					isVisited[*itGroupMembers]=true;
				}

				if (groupMembers.size()< sizeThresh )
				{
					cout<<"Delegate:"<<pointIndex<<" Label:"<<pointData[pointIndex]<<" GroupSize:"<<groupMembers.size()<<endl;
					cout<<"Remove this region and Relabel it using neighbor label: ";
					size_t maxNbr=0;
					int maxNbrLabel=0;
					for(map<int, list<int> >::iterator itNbrLabels= nbrLabels.begin(); itNbrLabels!=nbrLabels.end(); itNbrLabels++)
					{
						if (itNbrLabels->second.size() > maxNbr)
						{
							maxNbr=itNbrLabels->second.size();
							maxNbrLabel=itNbrLabels->first;
						}					
					}
					if (0==maxNbr)
					{
						cout<<"Warning! Isolated vertex: "<<pointIndex<<endl;
						continue;
					}				

					cout<<maxNbrLabel<<endl; 
					cout<<"Neighbors of this kind: "<<maxNbr<<endl;
					for(list<int>::iterator itGroupToRelabel= groupMembers.begin(); itGroupToRelabel!=groupMembers.end(); itGroupToRelabel++)
					{
						pointData[*itGroupToRelabel]=maxNbrLabel;
					}

				}		

			}


		} //end of for;
	}

	void CTriSurface::ReadPointLabel( string dataName, int skipLines )
	{
		data->savePointsLabel = true;
		fstream oriSrcFile;
		oriSrcFile.open(data->inputFileName.c_str(), ios::in);
		if (NULL==oriSrcFile)
		{
			cerr<<"can not open file "<<data->inputFileName.c_str()<<endl;
			exit(EXIT_FAILURE);
		}
		// skip the lines of vertex and faces;
		if(-1 == skipLines){
			skipLines = data->numOfCells + data->numOfPoints;
		}
		string line;
		while(skipLines--){
			getline(oriSrcFile, line);
		}
		// read each line and try finding the matching property name;
		do{
			getline(oriSrcFile, line);
		} while(line.find(dataName) == string::npos && !oriSrcFile.eof());
		getline(oriSrcFile, line); // lookup tabel line;
		for(int pointIndex = 0;pointIndex < data->numOfPoints;pointIndex++){
			oriSrcFile >> data->pointsLabel[pointIndex];
		}
		oriSrcFile.close();
	}

	void CTriSurface::ReadCellLabel( string dataName, int skipLines/*=-1*/ )
	{
		data->saveCellsLabel = true;
		fstream oriSrcFile;
		oriSrcFile.open(data->inputFileName.c_str(), ios::in);
		if (NULL==oriSrcFile)
		{
			cerr<<"can not open file "<<data->inputFileName.c_str()<<endl;
			exit(EXIT_FAILURE);
		}
		// skip the lines of vertex and faces;
		if(-1 == skipLines){
			skipLines = data->numOfCells + data->numOfPoints;
		}
		string line;
		while(skipLines--){
			getline(oriSrcFile, line);
		}
		// read each line and try finding the matching property name;
		do{
			getline(oriSrcFile, line);
		} while(line.find(dataName) == string::npos && !oriSrcFile.eof());
		getline(oriSrcFile, line); // lookup tabel line;
		for(int cellIndex = 0;cellIndex < data->numOfCells;cellIndex++){
			oriSrcFile >> data->cellsLabel[cellIndex];
		}
		oriSrcFile.close();
	}

	void CTriSurface::ReadPointColor( string vectorName, int skipLines/*=-1*/ )
	{
		data->savePointsColor = true;
		fstream oriSrcFile;
		oriSrcFile.open(data->inputFileName.c_str(), ios::in);
		if (NULL==oriSrcFile)
		{
			cerr<<"can not open file "<<data->inputFileName.c_str()<<endl;
			exit(EXIT_FAILURE);
		}
		// skip the lines of vertex and faces;
		if(-1 == skipLines){
			skipLines = data->numOfCells + data->numOfPoints;
		}
		string line;
		while(skipLines--){
			getline(oriSrcFile, line);
		}
		// read each line and try finding the matching property name;
		do{
			getline(oriSrcFile, line);
		} while(line.find(vectorName) == string::npos && !oriSrcFile.eof());
		/*getline(oriSrcFile,line); // lookup table line; */
		// no lookup table for vector data;

		for(int pointIndex = 0;pointIndex < data->numOfPoints;pointIndex++){
			oriSrcFile >> data->pointsColor[pointIndex].x;
			oriSrcFile >>  data->pointsColor[pointIndex].y;
			oriSrcFile >>  data->pointsColor[pointIndex].z;
		}
		oriSrcFile.close();
	}

void CTriSurface::SmoothPointDataScalar(std::vector< float >& pattern, int range, SmoothType mode, int iterN)
{ 
  
  while(iterN--)
  {
    cout<<"Iterations Left"<<": "<<iterN+1<<endl;
    vector<float> tmpPattern=pattern;
    if (mode== CTriSurface::Voting ) // vote; 
    {

      int min=10000000; 
      int max=-1000000; 
      for (int i=0; i<this->GetNumOfPoints(); ++i)
      {	
	if(min > tmpPattern[i])
	  min= tmpPattern[i]; 
	if(max < tmpPattern[i])
	  max= tmpPattern[i]; 
      }
      
      //set number of bins; 
      int binSize= max-min+1;       
   
      for (int i=0; i< this->GetNumOfPoints(); ++i)
      {
	vector<int> histogram( binSize,0 ); 
	list<int> allNbrs; 
	this->GetNeighboringPointsOfPoint(i,allNbrs,range);

	//find the voting results; 	
	for( list<int>::const_iterator it= allNbrs.begin(); it!= allNbrs.end(); it++)
	{
	  int pointID= *it; 
	  ++histogram[tmpPattern[pointID]-min]; 
	}

	  int index=-1; 
	  int maxValue=0; 
	  findMaxValueAndIndex(histogram,maxValue,index);
	  pattern[i]=index+ min; 
      }    
      
    }
    else if (mode== CTriSurface::Averaging ) // average;
    {
      for (int i=0; i<this->GetNumOfPoints(); ++i)
      {
	      double allValue=0;
	      KML::PointIDsType nbrPoints; 
	      this->GetNeighboringPointsOfPoint(i,nbrPoints,range);
	      for ( KML::PointIDsType::const_iterator itNbrPoints=nbrPoints.begin(); itNbrPoints!= nbrPoints.end(); itNbrPoints++ )
	      {
		      allValue+=tmpPattern[*itNbrPoints];
	      }

	      pattern[i]=(float) (allValue/nbrPoints.size() );
      }

    }
    else
    {
      cerr<<"Unsupported smooth type! "<<endl;
      exit(EXIT_FAILURE);
    }
 
  
  }//end of while; 


}

int CTriSurface::FindClosestSurfaceID(VectorType& coord)
{
  int theId =-1; 
  float minDis =1000000; 
  float currentDis=0;
  for(int idxSurf = 0; idxSurf < data->numOfPoints; ++idxSurf )
  {  
    VectorType dis = coord - data->points[idxSurf]; 
    currentDis = dis.Norm1stOrder();
    if(minDis> currentDis)
    {
      minDis = currentDis;
      theId= idxSurf;
    }
  } // end for::idxSurf;
  
  return theId;
}

} //end of KML ; 


