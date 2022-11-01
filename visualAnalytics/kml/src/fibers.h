/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#ifndef  __FIBERS__H
#define  __FIBERS__H

#include "Vector3D.h"
#include <map>
#include <string>
#include <set>
#include <vector>
#include <list>

using namespace std;


namespace KML
{
	typedef Vector3D<float> PointCoordType;
	
	class CIntersectionType;
	class CFibersPrivacy; 
	class CFibers
	{
		private:
			CFibersPrivacy* m_data; 
			int ShortenCurrentFiber(int fid, CIntersectionType& front, CIntersectionType& rear, float dis);
		public:
			CFibers(string fileName );
			CFibers(CFibers& fibers);
			~CFibers();
			void SwapEndians();
	
			vector<int>& GetFiber(int fiberId);
			float GetFiberLength(int fiberId) const;
			int GetFiberSize(int fiberId) const; 
			const string& GetFiberName()const;
			PointCoordType& GetPointCoords(int pointId);
			void SetPointCoords(int pointId, PointCoordType& newCoords);
			
			PointCoordType& GetPointColorAscii( int pointId); 
			KML::Vector3D< unsigned char >& GetPointColorBinary(int pointId); 
			
			int GetNumFibers(void);
			int GetNumPoints(void);
			void SaveFibers(list<int>& fiberIds,string fileName,bool toAsciiMode=false);
			void SaveAllFibers(string fileName,bool save2AsciiMode=false);
			void SaveFibersWithLabel(string fileName, int labelValue );
			void ReadFiberLabel(string fileName);
			vector<int>& GetALLFiberLabel(void);
			int& GetFiberLabel(int fiberId);
			void SetFiberLabel(int fiberId, int fiberLabel);
			void SetFiberLabel(vector<int>& allLables);
			void IsSwapNeeded(bool needSwap);
			void Translate(float x, float y, float z);
			void SaveFibers(set<int>& fiberIds, string fileName);
			void SaveFibers( const set<int>& fiberIds, string fileName);
			void SaveFibers( const vector<int>& fiberIds, string fileName);
// 			void ShortenFibers(CTriSurface& mySurface, float maxDis=0.6);
			void ShortenFibers(int numPoints=3);
			//the max length to extend, default is 4mm;
// 			void ExtendFibers(CTriSurface& mySurface,float max2Extent= 4);
// 			void RemoveInvalidFibers(CTriSurface& mySurface);
			void NeedLineColor(void);
			void SetLineColor(int pid, PointCoordType& theColor);
			
			//resample points, 1/resampelDegree; 
			void ResamplePoints(int resampelDegree=4); 
			void ResampleLines(int resampelDegree=3);
			
			//remove short fibers; 
			void RemoveShortFibers(float minFiberLength=20);
			//get the resolution of fiber, i.e. the distance between two points; 
			float GetResolution(void); 
	};


};

#endif
