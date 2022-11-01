/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#pragma once

#include <vector>
#include <set>
//#include <omp.h>
#include "triSurface.h"
#include "fibers.h"
#include "Vector3D.h"
#include "newimageall.h"
#include "newmat.h"
#include "indexer.h"
#include "kaimingCommon.h"

using namespace std;
using namespace KML;
using namespace NEWMAT;
using namespace NEWIMAGE;

namespace KML
{

	class CJointModelPrivacy;
	class CJointModel;

	// for intersection of line with cells;
	class CIntersectionType
	{
		friend class CJointModelPrivacy;
		friend class CJointModel;
		bool isValid; // false means this intersection is not used or means nothing;
		int cellId;
		float coord[3];
		Vector3D<size_t> index;
	public:
		CIntersectionType(){isValid=false;}
		bool IsValid(void){ return isValid; }
		int GetCellId(void) const {return cellId;}
		Vector3D<float> GetCoord(void) { return Vector3D<float> (coord[0],coord[1],coord[2]) ; }
		float* GetCoordFloat(void) {return coord; }
		const Vector3D<size_t>& GetIndex(void) { return index; }
		~CIntersectionType(){};
	};
	// for joint modeling using dti, fmri. 
	class CJointModel
	{
	private:
		CJointModelPrivacy* data;
		CJointModel(const CJointModel&);
		void operator = (const CJointModel&);
		void CenteringBolds(void); 
		int  GetGMBoldsByNormal(int pointId, int nbrSize, vector<float>& timeSeries, float serachStep=1., int maxStep=15);
		int  GetGMBoldsByFiber(int pointId, int fiberId, vector<float>& timeSeries, int nbrSize=1 , int firstPointsNum=15);
		void GetImageCoordDti( Vector3D<float>& phisicalCoord, Vector3D<int>& imageCoord); 
		void ReLabelClusters(int tobeKill, int dst) ;
		void AssignColor2Cluster(Vector3D<float>& color, int clusterLabel);
		bool IsVertexRightMode(int pointId) const;
		float GetGaussianKenerlValue(float x, float sigma);
                void InitData(string surfaceName);


	public:
		CJointModel(string surfaceName, string fibersName, string boldsName, string dtiSegName );
		CJointModel (string surfaceName, string fibersName );
		CJointModel(CTriSurface& mySurface, CFibers& myFiber);
		CJointModel(CTriSurface& mySurface, CFibers& myFiber, volume4D<float>& myBolds, volume<float>& myGM); 
		~CJointModel();


		//find the front and rear intersection between fiber and surface;
		void FindFiberSurfIntersections(void);
		void SaveALLIntersectionPoints(string fileName) ;
		void SavePartialIntersectionPoints(const vector<int>& fidLists, string fileName);
		void SavePartialIntersectionPoints(const set<int>& fidList, string fileName) ;

		//for each cell , get a list of intersecting fibers; note FindFiberSurfIntersections must be called
		// to find the fiber ids on a certain cell;
		void FindFibersOnCells(void);
		void WarpingBoldsOnSurface(void);
		void GeometricSmoothing(void);

		// a group of gets and sets;
		vector<CIntersectionType>& GetFrontIntersectingCellOfFibers(void);
		vector<CIntersectionType>& GetRearIntersectingCellOfFibers(void);
		vector<vector<int> >& GetIntersectingFibersOfCells(void);

		CIntersectionType& GetFrontIntersectingCellOfFiber(int fiberId);
		CIntersectionType& GetRearIntersectingCellOfFiber(int fiberId);
		vector<int>& GetIntersectingFibersOfCell(int cellId);

		vector<float>& GetFoldingPattern(void);
		int GetNumSurfacePoints(void) const; 
		int GetNumSurfaceCells(void) const; 
		int GetNumSurfaceGyriPoints(void) const; 
		int GetNumSurfaceSulciPoints(void) const; 
		int GetNumTimeSeries(void) const; 
		int GetNumFibers(void) const; 	
		int GetVolume2Del(void) const; 
		vector<float>& GetVertexTimeSeris(int vertexId) const;
		float GetAbsCorrelThresh(void) const; 
		Vector3D<float>& GetOffSet(void) const; 
		float GetBoldDimX(void) const; 
		float GetBoldDimY(void) const; 
		float GetBoldDimZ(void) const; 
		int GetNumSamples(void) const;
		float GetSmoothingSigma(void) const;
		int GetSmoothingRange(void) const;
		
		
		CTriSurface& GetSurface(void);
		CFibers& GetFibers(void);
		volume4D<float>& GetBOLDs(void); 
		volume<float>& GetGM(void); 
		
		// get the fiber density on the cortical surface;
		void CalculateFiberDensity(vector<float>& dst,int range=1);
		
		void SetOffSet(Vector3D<float>& offset); 
		void SetOffSet(float x, float y, float z); 
		void SetVolume2Del(int vol2del); 
		void SetAbsCorrelThresh(float thresh); 
		void SetNumSamples(int num);
		void SetMode(int m);
		void SetTimePointRange4SimCalc(int start, int end);
		void SetSmoothingSigma(float sigma);
		void SetSmoothingRange(int range);
                //deprecated, for compatability only; 
//                 void SetIndexer(CIndexer* pIndexer);
		void SetIndexer(CIndexer theIndex);
		CIndexer& GetIndexer(void);
		//index the intersections, FindFiberSurfIntersections must be called before;
		void IndexingIntersectionsOnCell(void);

		void CalculateFibersInsideCubics(void);
		void SaveFibersInsideCubic(Vector3D<size_t>&  index, string fileName);


		void GetCorrelation(list<int>& v1,list<int>& v2, list<float>& correl,vector<unsigned int>& indexMapID);
		void ClassifyUnselectedVertices(vector<unsigned int>& clusterIndexMapSurfacePointID, vector<unsigned int>& centersId);
		void MakeUpMisClassification(string cingulateFile);
		// merge clusters that are highly correlated;
		void MergeCluters(vector<unsigned int>& centersId, vector<unsigned int>& newCenterIds, float simThreshold=0.8);
		// smoothing clustering result by voting
		void SmoothingByVoting(int iterN=1);
		//coloring top number of clusters;
		void ColoringTopClusters(int num);
		//set background color if necessary;
		void SetBackGroundColor(Vector3D<float>& backgroundColor);
		//prune small clusters;
		void PruneSmallClusters(vector<unsigned int>& centersId,int sizeThresh=50);
		// remove small isolated buble
		void RemoveSmallIsolatedRegion(int sizeThresh=30);

		void Test(void);

		//roi fibers;
		void GetFibersInsideROI(volume<float>& roiIndtiSpace, list<int>& resutls);
		void GetFibersInsideROI(volume<float>& roiIndtiSpace, list<int>& resutls, int roiLabelValue);
		void GetFibersInsideROI(volume<float>& roiIndtiSpace, list<int>& fiberConsidered, list<int>& results);
		void GetFibersInsideROI(volume<float>& roiIndtiSpace, list<int>& fiberConsidered, list<int>& results, int roiLabelValue);
                void GetFibersByVerticeList(list<int>& vertices, vector<int>& result);
                
                ///deprecated method, backward compatability;
		void GetFibersBySidNeighborhood( int fid, int ring, set<int>& results);
                
                void GetFibersBySidNeighborhood( int fid, int ring, vector<int>& results);
                
                ///deprecated method, backward compatability;
                void GetDajiangFibersBySidNeighborhood( int fid, int ring, set<int>& results);
                
                void GetDajiangFibersBySidNeighborhood( int fid, int ring, vector<int>& results);

		const set<int>& GetFibersInsideCubeFront( int fid);
		const set<int>& GetFibersInsideCubicByCenterCoord(const Vector3D<float>& centerCoord);
                const set<int>& GetFibersInsideCubicByIndexGrid(const Vector3D<size_t>& grid);
		//find invalid fibers, and save them into a file;
		void GetUnIntersectedFibers(string fileName);
                
                //GetBOLDs for bad mapping, fiber/norm doesnt work; 
                int GetBolds4Failures(int pointId, std::vector< float >& timeSeries, int nbrSize, int fiberId = -1, int firstPointsNum=15);

                //find the intersection of surface and grid, save in a volume file; note: a indexer must be set; 
                void FindSurfaceVolumeIntersection(void); 
		//check whether or not a intersection grid; 
		bool IsAIntersectionGrid(Vector3D<size_t>& grid);
                
                //build index for fibers and voxels at the same time, a fiber will have a series of grids that are along the fiber; and edge voxel will have a series of fibers; note: indexer must be set and FindSurfaceVolumeIntersection must be called before; ( will call automatically)
                void BuildFiberEdgeIndiceWIndexer(void);
                
                //save the volume and surface intersection; 
                void SaveVolSurfX(string fileName); 
                


	};
};
