/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#pragma once

#include "Vector3D.h"
#include "triSurfEdge.h"
#include <map>
#include <string>
#include <set>
#include <vector>
#include <list>
#include "kaimingCommon.h"
using namespace std;

namespace KML
{

	typedef Vector3D<float> PointType;
	typedef Vector3D<float> VectorType;
	typedef Vector3D<int> CellType;
	typedef Vector3D<float> PointDataType;
	typedef Vector3D<float> CellDataType;
	typedef list<int> PointIDsType;
	typedef list<int> CellIDsType;

	class CTriSurfacePrivacy;
	class CTriSurface
	{
	public:
		//////////////////////////////////////////////////////////////////////////
		enum SmoothType{ Averaging, Voting}; 
		CTriSurface( string filename);
		CTriSurface(CTriSurface& mySurface);
		~CTriSurface(void);

		void SaveAs( string filename);
		void Save(void);
		const string& GetFileName(void) const;


		int GetNumOfPoints(void)const;
		int GetNumOfCells(void)const;
		int GetNumOfEdges(void) const;


		// return edge if a edge id is given; 
		CEdgeInfoType& GetEdge(int lineId);
		set<int>& GetPointEdges(int pointId);
		vector<int>& GetCellEdges(int cellId);

		const VectorType& GetPointCoords(int pointId);
		const VectorType& GetPointNormal(int pointId) const;


		VectorType GetCellNormal(int cellId);
		float GetPointCurvature (int pointId);
		const list<int>&  GetNeighboringCellsOfPoint(int pointId );
		void GetNeighboringPointsOfPoint(int pointId, list<int>& nbrs, int range=1);
                
                ///backward compatibility, deprecated ;
		void GetNeighboringCellsOfPoint(int pointId, set<int>& nbrCells,int range=1);
                void GetNeighboringCellsOfPoints(list<int>& nbrs, set<int>& nbrCells);
                
                
                
                void GetNeighboringCellsOfPoint(int pointId, vector<int>& nbrCells,int range=1);
                void GetNeighboringCellsOfPoints(list<int>& nbrs, vector<int>& nbrCells);
                
                
		void GetNeighboringPointsOfPoint(int pointId, list<int>& nbrs, vector<float>&  data, float interestedValue, int range=1);
		void GetNeighboringPointsOfPoint(int pointId, list<int>& nbrs, vector<int>&  data, int interestedValue, int range=1);
		void GetNbrsOfPointOnRange(int pointId,list<int>& outerNbrs, int range);
		
		//isROI: false if you want the point to be considered. 
		void GetNeighboringPoints(std::list< int >& points, std::list< int >& nbrs, std::vector< bool >& isROI); 
		// propagate from a single point, until find all its group members(includes itself) or a maximum member is found. return its group members and the labels of the nbring classes; maxSize=-1 mean no limit to the groupsize;; 
		void Propagate2All(int pointId, list<int>& groupMembers, map<int,list<int> >&  nbrLabels,vector<int>& dataProperty, int maxSize=-1);
		void Propagate2All(int pointId, list<int>& groupMembers, map<float, list<int> >& nbrLabels,vector<float>& dataProperty,int maxSize=-1);


		vector<VectorType>& GetAllPointColor(void); 
		const VectorType& GetPointColor(int pointId)const;
		void SetPointColor(int pointId, const VectorType& floatColor);

		float GetPointComponent(int pointId, int component);
		void SetPointComponent(int pointId, int component, float value);

		// label points and cells;
		int GetPointLabel(int pointId);
		vector<int>& GetAllPointLabel(void);
		void SetPointLabel(int pointId, int label);
		void  GetPointsWithLabel(int label, list<int>& target);
		int GetCellLabel(int cellId) ;
		vector<int>& GetAllCellLabel(void);
		void SetCellLabel(int cellId, int label);
		void  GetCellsWithLabel(int label, list<int>& target) ;

		float GetCellArea(int cellId);
		float GetCellsTotalArea(const CellIDsType& cellIds);
		float GetCellsTotalArea(int startId, int numOfCells);
		const CellType& GetCell(int cellId) ;

		double GetCurvatureThreshold(void) ;
		void SetCurvatureThreshold(double thresh);

		void NeedPointsColor(bool bColor=true);
		void NeedPointsNormal(bool bNormal=true);
		void NeedPointsCurvature(bool bCurvature=true);
		void NeedPointsLabel(bool bLabel=true);
		void NeedCellsLabel(bool bLabel=true);

		void AddPointDataScalar(vector<float>& src, string dataName);
		void AddPointDataVector(vector<VectorType>& src, string vectorName);
		void AddCellDataScalar(vector<float>& src, string dataName);
		void AddCellDataVector(vector< VectorType>& src, string vectorName);

		void ReadPointDataScalar(string dataName, int skipLines=-1);
		void ReadPointDataVector(string vectorName, int skipLines=-1);
		void ReadCellDataScalar(string dataName, int skipLines=-1);
		void ReadCellDataVector(string vectorName, int skipLines=-1);
		void ReadPointLabel(string dataName, int skipLines=-1);
		void ReadPointColor(string vectorName, int skipLines=-1);
		void ReadCellLabel(string dataName, int skipLines=-1);
		void ReadCellColor(string vectorName, int skipLines=-1);

		vector<float>& GetPointDataScalar(string dataName);
		vector<int>& GetPointDataLabel(void);
		vector<int>& GetCellDataLable(void);
		vector<VectorType>& GetPointDataColor(void);
		vector<float>& GetCellDataScalar(string dataName);
		vector<VectorType>& GetPointDataVector(string dataName);
		vector<VectorType>& GetCellDataVector(string dataName);


		void WriteLineBorders(vector<bool>& tag, string fileName);

		//get neighbor information for points and cells;
		void BuildNeighbors(void);
		//get edge information;
		void BuildEdge(void);

		//warp point label to cell label;
		void WarpPointLabel2CellLabel(void);
		void SmoothPointLabel(int range, int numBins, int itN=3);
		void SmoothCellLabelBySimpleMajority(int itN = 3);
		void SmoothPointDataScalar(vector<float>& data, int range, SmoothType mode, int itN=3);
		void SmoothCellDataScalar(vector<float>& data, int range, int numBins, int itN=3);
		void RemoveSmallIsolatedLabeledRegion( int sizeThresh=30 );
                
                //find the closed sid for a coords; 
                int FindClosestSurfaceID(VectorType& coord);

	private:

		CTriSurfacePrivacy* data;
		void InitData(void);
		void GeneratePointNormal(void);
		void GeneratePointCurvature(void);
		CTriSurface(void);

	};

};
