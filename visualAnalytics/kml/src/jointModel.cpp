/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#include "jointModel.h"
#include "ColorSchemeRGB.h"
#include "vtkCellLocator.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkGenericCell.h"
#include <cstdlib>
#include <ctime>
#include <map>
#include <set>
#include "algorithm"
#include "newimageio.h"
#include "newimage.h"
#include "newimagefns.h"
#include "kaimingCommon.h"

using namespace KML;
using namespace NEWIMAGE; 

namespace KML
{
enum Mode
{
    GYRI,
    SULCI,
    ALL
};


class CJointModelPrivacy
{
    friend class CJointModel;

private:

    string surfaceName;
    string fibersName;
    string boldsName;
    string dtiSegName;

    CTriSurface surface;
    CFibers fibers;
    volume4D<float> bolds;
    volume<float> dtiSeg;

    vector<vector<float> > vertexTimeSeries;
    vector<vector<vector<set<int> > > > fibersInsideCubics;
    vector<float> foldingPattern;
    Vector3D<float> offSet;

    int numSurfacePoints;
    int numSurfaceGyriPoints;
    int numSurfaceSulciPoints;
    int numSurfaceCells;
    int numFibers;
    int numTimeSeries;
    int numSamples;
    int vol2Del;
    int numClusters;
    int timeSeriesStart;
    int timeSeriesEnd;
    float xdimfMRI;
    float ydimfMRI;
    float zdimfMRI;
    float xdimDti;
    float ydimDti;
    float zdimDti;
    Mode mode;
    float correlThresh;
    bool needSampling;
    float sigma;
    int smoothRange;

    vector<CIntersectionType> fiberIntersectionCellFront;// keep record the intersection cells of one end;
    vector<CIntersectionType> fiberIntersectionCellRear; // keep record the intersection cells of the other end;
    vector<vector<int> > cellIntersectingFibers;  // keep all the fibers that intersect with the cell;

    //for intersection locator;
    vtkPolyDataReader* surfaceReader;
    vtkPolyData* surfaceData;
    vtkCellLocator* locator;

    //for indexer
    CIndexer indexer;

    //this is a volume file has the same space with the indexer, and keeep records of the intersections of volume and surface, x for intersection;
    volume<short> SurfVolX;
    //for each fiber, record a series of grids/voxels along the fiber; 
    vector<vector<Vector3D<size_t> > > grids4EachFiber;
    //the each intersecting voxel, keep record of the fibers it contains;
    vector<vector<vector<set<int> > > > fibersInEdgeGrids;  
    //for each intersecting voxel, keep record of the surface vertices it contains; 
    vector<vector<vector<list<int> > > > verticesInEdgeGrids; 
    //for staus;
    bool  isDoneFiberSurfIntersection;
    bool isDoneIntersectingFibers4EachCell;
    bool isDoneCalculateFibersInsideCubic;
    bool isDoneFindSurfaceVolumeIntersection;
    bool isDoneBuildFiberEdgeIndiceWIndexer; 


public:
    CJointModelPrivacy( string surfaceName, string fibersName ): surface(surfaceName), fibers(fibersName), isDoneFindSurfaceVolumeIntersection(false),isDoneFiberSurfIntersection(false),isDoneIntersectingFibers4EachCell(false),isDoneCalculateFibersInsideCubic(false),isDoneBuildFiberEdgeIndiceWIndexer(false) {};
    CJointModelPrivacy(CTriSurface& mySurface, CFibers& myFiber): surface(mySurface), fibers(myFiber), isDoneFindSurfaceVolumeIntersection(false),isDoneFiberSurfIntersection(false),isDoneIntersectingFibers4EachCell(false),isDoneCalculateFibersInsideCubic(false),isDoneBuildFiberEdgeIndiceWIndexer(false) {};

};


CJointModel::CJointModel( string surfaceName, string fibersName, string boldsName, string dtiGMName )
{
    data = new CJointModelPrivacy(surfaceName,fibersName );

    data->surfaceName=string(surfaceName);
    data->fibersName=string(fibersName);
    data->boldsName=string(boldsName);
    data->dtiSegName=string(dtiGMName);
    
    this->InitData(surfaceName);
 
    data->boldsName=string(boldsName);
    data->dtiSegName=string(dtiGMName);
    read_volume(data->dtiSeg,data->dtiSegName);
    data->xdimDti= data->dtiSeg.xdim();
    data->ydimDti= data->dtiSeg.ydim();
    data->zdimDti= data->dtiSeg.zdim();


    read_volume4D(data->bolds,data->boldsName);
    this->CenteringBolds();

    data->numTimeSeries=data->bolds.tsize();
    data->timeSeriesStart=0;
    data->timeSeriesEnd=data->numTimeSeries;
    data->xdimfMRI=data->bolds.xdim();
    data->ydimfMRI=data->bolds.ydim();
    data->zdimfMRI=data->bolds.zdim();

}


void CJointModel::InitData(string surfaceName)
{
    data->surface.BuildNeighbors();
    data->surface.NeedPointsLabel(true);
    data->surface.NeedPointsColor(true);
    data->numFibers=data->fibers.GetNumFibers();
    data->numSurfacePoints=data->surface.GetNumOfPoints();
    data->numSurfaceCells=data->surface.GetNumOfCells();
    data->numSamples=20000;
    data->vol2Del=0;
    data->numClusters=0;
    data->correlThresh=0.3;
    data->mode=ALL;
    data->needSampling= true;
    data->sigma=4;
    data->smoothRange=3;

    data->fiberIntersectionCellFront.resize(data->numFibers);
    data->fiberIntersectionCellRear.resize(data->numFibers);
    data->cellIntersectingFibers.resize(data->numSurfaceCells);
    data->vertexTimeSeries.resize(data->numSurfacePoints);
    data->foldingPattern.resize(data->numSurfacePoints);

    data->surface.ReadPointDataScalar("SulciORGyri");
    data->foldingPattern=data->surface.GetPointDataScalar("SulciORGyri");
    data->numSurfaceGyriPoints=0;
    for (int pointIndex=0; pointIndex<data->numSurfacePoints; ++pointIndex)
    {
        if (-1==data->foldingPattern[pointIndex])
            ++data->numSurfaceGyriPoints;
    }
    data->numSurfaceSulciPoints=data->numSurfacePoints-data->numSurfaceGyriPoints;

    data->surface.NeedPointsNormal(true);


    data->surfaceReader=vtkPolyDataReader::New();
    data->surfaceReader->SetFileName(surfaceName.c_str());
    data->surfaceReader->Update();

    data->surfaceData=vtkPolyData::New();
    data->surfaceData=data->surfaceReader->GetOutput();

    data->locator = vtkCellLocator::New();
    data->locator->SetDataSet(data->surfaceData);
    data->locator->CacheCellBoundsOn();

    data->locator->AutomaticOn();
    data->locator->BuildLocator();

    this->FindFiberSurfIntersections();
    this->FindFibersOnCells();
}

CJointModel::CJointModel( string surfaceName, string fibersName )
{
    data = new CJointModelPrivacy(surfaceName,fibersName );

    data->surfaceName=string(surfaceName);
    data->fibersName=string(fibersName);
    
    this->InitData(surfaceName);
}


CJointModel::CJointModel(CTriSurface& mySurface, CFibers& myFiber) {

    data = new CJointModelPrivacy(mySurface,myFiber );
    data->surface.BuildNeighbors();
    data->surface.NeedPointsLabel(true);
    data->surface.NeedPointsColor(true);
    data->surfaceName=mySurface.GetFileName();
    data->fibersName=myFiber.GetFiberName();

    data->numFibers=data->fibers.GetNumFibers();
    data->numSurfacePoints=data->surface.GetNumOfPoints();
    data->numSurfaceCells=data->surface.GetNumOfCells();
    data->numSamples=20000;
    data->vol2Del=0;
    data->numClusters=0;
    data->correlThresh=0.3;
    data->mode=ALL;
    data->needSampling= true;
    data->sigma=4;
    data->smoothRange=3;

    data->fiberIntersectionCellFront.resize(data->numFibers);
    data->fiberIntersectionCellRear.resize(data->numFibers);
    data->cellIntersectingFibers.resize(data->numSurfaceCells);
    data->vertexTimeSeries.resize(data->numSurfacePoints);
    data->foldingPattern.resize(data->numSurfacePoints);

    data->surface.ReadPointDataScalar("SulciORGyri");
    data->foldingPattern=data->surface.GetPointDataScalar("SulciORGyri");
    data->numSurfaceGyriPoints=0;
    for (int pointIndex=0; pointIndex<data->numSurfacePoints; ++pointIndex)
    {
        if (-1==data->foldingPattern[pointIndex])
            ++data->numSurfaceGyriPoints;
    }
    data->numSurfaceSulciPoints=data->numSurfacePoints-data->numSurfaceGyriPoints;

    data->surface.NeedPointsNormal(true);


    data->surfaceReader=vtkPolyDataReader::New();
    data->surfaceReader->SetFileName(mySurface.GetFileName().c_str());
    data->surfaceReader->Update();

    data->surfaceData=vtkPolyData::New();
    data->surfaceData=data->surfaceReader->GetOutput();

    data->locator = vtkCellLocator::New();
    data->locator->SetDataSet(data->surfaceData);
    data->locator->CacheCellBoundsOn();

    data->locator->AutomaticOn();
    data->locator->BuildLocator();

 
    this->FindFiberSurfIntersections();
    this->FindFibersOnCells();
}

CJointModel::CJointModel(CTriSurface& mySurface, CFibers& myFiber, volume4D< float >& myBolds, volume< float >& myGM)
{
    data = new CJointModelPrivacy(mySurface,myFiber );
    data->surface.BuildNeighbors();
    data->surface.NeedPointsLabel(true);
    data->surface.NeedPointsColor(true);
    data->surfaceName=mySurface.GetFileName();
    data->fibersName=myFiber.GetFiberName();

    data->numFibers=data->fibers.GetNumFibers();
    data->numSurfacePoints=data->surface.GetNumOfPoints();
    data->numSurfaceCells=data->surface.GetNumOfCells();
    data->numSamples=20000;
    data->vol2Del=0;
    data->numClusters=0;
    data->correlThresh=0.3;
    data->mode=ALL;
    data->needSampling= true;
    data->sigma=4;
    data->smoothRange=3;

    data->fiberIntersectionCellFront.resize(data->numFibers);
    data->fiberIntersectionCellRear.resize(data->numFibers);
    data->cellIntersectingFibers.resize(data->numSurfaceCells);
    data->vertexTimeSeries.resize(data->numSurfacePoints);
    data->foldingPattern.resize(data->numSurfacePoints);

    data->surface.ReadPointDataScalar("SulciORGyri");
    data->foldingPattern=data->surface.GetPointDataScalar("SulciORGyri");
    data->numSurfaceGyriPoints=0;
    for (int pointIndex=0; pointIndex<data->numSurfacePoints; ++pointIndex)
    {
        if (-1==data->foldingPattern[pointIndex])
            ++data->numSurfaceGyriPoints;
    }
    data->numSurfaceSulciPoints=data->numSurfacePoints-data->numSurfaceGyriPoints;

    data->surface.NeedPointsNormal(true);

    data->dtiSeg= myGM;
    data->bolds=myBolds;
    this->CenteringBolds();

    data->xdimDti= data->dtiSeg.xdim();
    data->ydimDti= data->dtiSeg.ydim();
    data->zdimDti= data->dtiSeg.zdim();

    data->numTimeSeries=data->bolds.tsize();
    data->timeSeriesStart=0;
    data->timeSeriesEnd=data->numTimeSeries;
    data->xdimfMRI=data->bolds.xdim();
    data->ydimfMRI=data->bolds.ydim();
    data->zdimfMRI=data->bolds.zdim();


    data->surfaceReader=vtkPolyDataReader::New();
    data->surfaceReader->SetFileName(mySurface.GetFileName().c_str());
    data->surfaceReader->Update();

    data->surfaceData=vtkPolyData::New();
    data->surfaceData=data->surfaceReader->GetOutput();

    data->locator = vtkCellLocator::New();
    data->locator->SetDataSet(data->surfaceData);
    data->locator->CacheCellBoundsOn();

    data->locator->AutomaticOn();
    data->locator->BuildLocator();
 
    this->FindFiberSurfIntersections();
    this->FindFibersOnCells();

}

CJointModel::~CJointModel()
{
    data->locator->Delete();
    data->surfaceReader->Delete();
    delete data;
}

void CJointModel::CenteringBolds( void )
{
    int sizeX=data->bolds.xsize();
    int sizeY=data->bolds.ysize();
    int sizeZ=data->bolds.zsize();
    int sizeT=data->bolds.tsize();

    for (int gmZ=0; gmZ<sizeZ; ++gmZ) {
        for (int gmY=0; gmY<sizeY; ++gmY) {
            for (int gmX=0; gmX<sizeX; ++gmX) {
                if (data->bolds(gmX,gmY,gmZ,0)) {
                    float meanValue=0;
                    for (int index=0; index<sizeT; ++index)
                    {
                        meanValue+=data->bolds(gmX,gmY,gmZ,index);
                    }
                    meanValue/=sizeT;

                    for (int index=0; index<sizeT; ++index)
                    {
                        data->bolds(gmX,gmY,gmZ,index)-=meanValue;
                    }

                }
            }
        }
    }//end of all fors;


}


void CJointModel::FindFiberSurfIntersections( void )
{
    if (! data->isDoneFiberSurfIntersection) {


        //#pragma omp  parallel  for default(shared)  schedule(dynamic) num_threads(4)
        for (int fiberIndex=0; fiberIndex<data->numFibers; ++fiberIndex)
        {
            //cout<<omp_get_thread_num()<<endl;
            int sub_id;
            vtkIdType cell_id;
            double param_t, intersect[3], paraCoord[3];
            double sourcePnt[3], destinPnt[3];


            vector<int>& currentFiber=data->fibers.GetFiber(fiberIndex);

            bool isIntersectionFound=false;
            int fiberPointStart=0;
            int fiberPointEnd=fiberPointStart+1;

            while (false==isIntersectionFound   && fiberPointEnd< currentFiber.size() )
            {
                int startId= currentFiber[fiberPointStart];
                int endId= currentFiber[fiberPointEnd];

                sourcePnt[0]= data->fibers.GetPointCoords(startId).x;
                sourcePnt[1]= data->fibers.GetPointCoords(startId).y;
                sourcePnt[2]= data->fibers.GetPointCoords(startId).z;

                destinPnt[0]= data->fibers.GetPointCoords(endId).x;
                destinPnt[1]= data->fibers.GetPointCoords(endId).y;
                destinPnt[2]= data->fibers.GetPointCoords(endId).z;

                int result=data->locator->IntersectWithLine(sourcePnt,destinPnt,(double)(0.001), param_t, intersect, paraCoord, sub_id, cell_id);
                if (1 == result)// has intersection
                {
                    isIntersectionFound = true;
                    CIntersectionType theIntersection;
                    theIntersection.isValid = true;
                    theIntersection.cellId = cell_id;
                    theIntersection.coord[0] = intersect[0];
                    theIntersection.coord[1] = intersect[1];
                    theIntersection.coord[2] = intersect[2];
                    data->fiberIntersectionCellFront[fiberIndex] = theIntersection;
                } else {
                    fiberPointStart += 1;
                    fiberPointEnd += 1;
                }
            } //end of while loop ;

        } //end of for 4 front;

        //#pragma omp  parallel  for default(shared)  schedule(dynamic) num_threads(4)
        for (int fiberIndex = 0;fiberIndex < data->numFibers;++fiberIndex) {
            int sub_id;
            vtkIdType cell_id;
            double param_t, intersect[3], paraCoord[3];
            double sourcePnt[3], destinPnt[3];
            vector<int> & currentFiber = data->fibers.GetFiber(fiberIndex);
            bool isIntersectionFound = false;
            int fiberPointStart = currentFiber.size() - 1;
            int fiberPointEnd = fiberPointStart - 1;
            while (false == isIntersectionFound && fiberPointEnd >= 0) {
                int startId = currentFiber[fiberPointStart];
                int endId = currentFiber[fiberPointEnd];
                sourcePnt[0] = data->fibers.GetPointCoords(startId).x;
                sourcePnt[1] = data->fibers.GetPointCoords(startId).y;
                sourcePnt[2] = data->fibers.GetPointCoords(startId).z;
                destinPnt[0] = data->fibers.GetPointCoords(endId).x;
                destinPnt[1] = data->fibers.GetPointCoords(endId).y;
                destinPnt[2] = data->fibers.GetPointCoords(endId).z;
                int result = data->locator->IntersectWithLine(sourcePnt, destinPnt, 0.001, param_t, intersect, paraCoord, sub_id, cell_id);
                if (1 == result)// has intersection
                {
                    isIntersectionFound = true;
                    CIntersectionType theIntersection;
                    theIntersection.isValid = true;
                    theIntersection.cellId = cell_id;
                    theIntersection.coord[0] = intersect[0];
                    theIntersection.coord[1] = intersect[1];
                    theIntersection.coord[2] = intersect[2];
                    data->fiberIntersectionCellRear[fiberIndex] = theIntersection;
                } else {
                    fiberPointStart -= 1;
                    fiberPointEnd -= 1;
                }
            } //end of while loop ;

        } //end of for 4 rear;

        data->isDoneFiberSurfIntersection=true;
    }

}

//This method saves the intersection points of fibers with surface;
void CJointModel::SaveALLIntersectionPoints(string fileName)
{
    if ( ! data->isDoneFiberSurfIntersection)
        this->FindFiberSurfIntersections();

    fstream outFile;
    outFile.open(fileName.c_str(), ios::out);
    if (NULL==outFile)
    {
        cerr<<"Can not create file :"<<fileName<<endl;
        exit(EXIT_FAILURE);
    }
    //write header
    outFile << "# vtk DataFile Version 2.0" << endl;
    outFile << "Brain surface " << endl;
    outFile << "ASCII" << endl;
    outFile << "DATASET POLYDATA" << endl;
    //write points;
    // calculate the number of intersections first;
    int intersectionPointsCount = 0;
    for (int fiberIndex = 0;fiberIndex < data->numFibers;++fiberIndex) {
        if (data->fiberIntersectionCellFront[fiberIndex].isValid == true) {
            ++intersectionPointsCount;
        }
        if (data->fiberIntersectionCellRear[fiberIndex].isValid == true) {
            ++intersectionPointsCount;
        }
    }

    outFile << "POINTS " << intersectionPointsCount << " FLOAT" << endl;
    for (int fiberIndex = 0;fiberIndex < data->numFibers;++fiberIndex) {
        if (data->fiberIntersectionCellFront[fiberIndex].isValid == true) {
            outFile << data->fiberIntersectionCellFront[fiberIndex].coord[0] << " " << data->fiberIntersectionCellFront[fiberIndex].coord[1] << " " << data->fiberIntersectionCellFront[fiberIndex].coord[2] << endl;
        }
        if (data->fiberIntersectionCellRear[fiberIndex].isValid == true) {
            outFile << data->fiberIntersectionCellRear[fiberIndex].coord[0] << " " << data->fiberIntersectionCellRear[fiberIndex].coord[1] << " " << data->fiberIntersectionCellRear[fiberIndex].coord[2] << endl;
        }
    }

    //write vertices;
    outFile << "VERTICES " << intersectionPointsCount << " " << intersectionPointsCount * 2 << endl;
    for (int pointIndex = 0;pointIndex < intersectionPointsCount;++pointIndex) {
        outFile << 1 << " " << pointIndex << endl;
    }
    outFile.close();
}

void CJointModel::SavePartialIntersectionPoints(const vector<int>& fidList, string fileName)
{
    if ( ! data->isDoneFiberSurfIntersection)
        this->FindFiberSurfIntersections();

    fstream outFile;
    outFile.open(fileName.c_str(), ios::out);
    if (NULL==outFile)
    {
        cerr<<"Can not create file :"<<fileName<<endl;
        exit(EXIT_FAILURE);
    }
    //write header
    outFile << "# vtk DataFile Version 2.0" << endl;
    outFile << "Brain surface " << endl;
    outFile << "ASCII" << endl;
    outFile << "DATASET POLYDATA" << endl;
    //write points;
    // calculate the number of intersections first;
    int intersectionPointsCount = 0;
    for (int fiberIndex = 0;fiberIndex < fidList.size();++fiberIndex) {
        if (data->fiberIntersectionCellFront[fidList[fiberIndex]].isValid == true) {
            ++intersectionPointsCount;
        }
        if (data->fiberIntersectionCellRear[fidList[fiberIndex]].isValid == true) {
            ++intersectionPointsCount;
        }
    }

    outFile << "POINTS " << intersectionPointsCount << " FLOAT" << endl;
    for (int idx = 0;idx < fidList.size();++idx) {

        int fiberIndex= fidList[idx];
        if (data->fiberIntersectionCellFront[fiberIndex].isValid == true) {
            outFile << data->fiberIntersectionCellFront[fiberIndex].coord[0] << " " << data->fiberIntersectionCellFront[fiberIndex].coord[1] << " " << data->fiberIntersectionCellFront[fiberIndex].coord[2] << endl;
        }
        if (data->fiberIntersectionCellRear[fiberIndex].isValid == true) {
            outFile << data->fiberIntersectionCellRear[fiberIndex].coord[0] << " " << data->fiberIntersectionCellRear[fiberIndex].coord[1] << " " << data->fiberIntersectionCellRear[fiberIndex].coord[2] << endl;
        }
    }

    //write vertices;
    outFile << "VERTICES " << intersectionPointsCount << " " << intersectionPointsCount * 2 << endl;
    for (int pointIndex = 0;pointIndex < intersectionPointsCount;++pointIndex) {
        outFile << 1 << " " << pointIndex << endl;
    }
    outFile.close();
}

void CJointModel::SavePartialIntersectionPoints(const set<int>& fidList, string fileName) {
    vector<int> list( fidList.size());
    copy( fidList.begin(), fidList.end(), list.begin());
    this->SavePartialIntersectionPoints( list, fileName);
}
// to find the fiber ids on a certain cell;
void CJointModel::FindFibersOnCells(void)
{
    if ( ! data->isDoneFiberSurfIntersection)
        this->FindFiberSurfIntersections();

    if ( ! data->isDoneIntersectingFibers4EachCell) {

        for (int fiberId = 0;fiberId < data->numFibers;++fiberId) {
            if (data->fiberIntersectionCellFront[fiberId].isValid) {
                int cellId = data->fiberIntersectionCellFront[fiberId].cellId;
                data->cellIntersectingFibers[cellId].push_back(fiberId);
            }
            if (data->fiberIntersectionCellRear[fiberId].isValid) {
                int cellId = data->fiberIntersectionCellRear[fiberId].cellId;
                data->cellIntersectingFibers[cellId].push_back(fiberId);
            }
        }

        data->isDoneIntersectingFibers4EachCell=true;
    }
}

void CJointModel::WarpingBoldsOnSurface(void)
{
    float badFiberWarping = 0;
    float allFiberWarping = 0;
    float badNormalWarping = 0;
    float allNormalWarping = 0;
    float totalFailure = 0;
    vector<int> failureROIs;
    //#pragma omp parallel for schedule( dynamic ) default (shared) num_threads(4)
    for (int pointIndex = 0;pointIndex < data->numSurfacePoints;++pointIndex) {
        //for each vertex, find all its one ring nbr,
        //if there are fibers intersecting with its nbrs, then find closest one, and use its direction information for gm time series;
        //else, find the closest gray matters in normal direction; define a nbrhood, averaging the wm.
        data->vertexTimeSeries[pointIndex].resize(data->numTimeSeries - data->vol2Del, 0);
        const list<int> & nbrCells = data->surface.GetNeighboringCellsOfPoint(pointIndex);
        if (nbrCells.size() == 0)// should not happen, just in case;
        {
            cerr << "Warning: isolated vertex found, id: " << pointIndex << endl;
            continue;
        }
        // pointIndex coords;
        float x = data->surface.GetPointComponent(pointIndex, 0);
        float y = data->surface.GetPointComponent(pointIndex, 1);
        float z = data->surface.GetPointComponent(pointIndex, 2);
        float minDis = 10000;
        int closestFiberId = -1;
        for (list<int>::const_iterator itNbrCells = nbrCells.begin();itNbrCells != nbrCells.end();itNbrCells++) {
            int cellId = *itNbrCells;
            if (data->cellIntersectingFibers[cellId].size())//this cell has intersecting fibers;
            {
                for (int fibersIndexOnThisCell = 0;fibersIndexOnThisCell < data->cellIntersectingFibers[cellId].size();++fibersIndexOnThisCell) {
                    int fiberId = data->cellIntersectingFibers[cellId][fibersIndexOnThisCell];
                    if (data->fiberIntersectionCellFront[fiberId].isValid) {
                        float intersectX = data->fiberIntersectionCellFront[fiberId].coord[0];
                        float intersectY = data->fiberIntersectionCellFront[fiberId].coord[1];
                        float intersectZ = data->fiberIntersectionCellFront[fiberId].coord[2];
                        float disTmp = abs(x - intersectX) + abs(y - intersectY) + abs(z - intersectZ);
                        if (disTmp < minDis) {
                            minDis = disTmp;
                            closestFiberId = fiberId;
                        }
                    }
                    else {
                        float intersectX = data->fiberIntersectionCellRear[fiberId].coord[0];
                        float intersectY = data->fiberIntersectionCellRear[fiberId].coord[1];
                        float intersectZ = data->fiberIntersectionCellRear[fiberId].coord[2];
                        float disTmp = abs(x - intersectX) + abs(y - intersectY) + abs(z - intersectZ);
                        if (disTmp < minDis) {
                            minDis = disTmp;
                            closestFiberId = fiberId;
                        }
                    }

                }

            }

        } // end of all the point's nbr cells;

        int result;
        if (-1 == closestFiberId)// fiber not found; using normal direction for searching;
        {
            ++allNormalWarping;
            result = this->GetGMBoldsByNormal(pointIndex, 0, data->vertexTimeSeries[pointIndex]);
            if (-1 == result)
                result = this->GetGMBoldsByNormal(pointIndex, 1, data->vertexTimeSeries[pointIndex]); // search a larger nbrhood.

            if (-1 == result)
                ++badNormalWarping;

        } else// using the closed fiber to search;
        {
            ++allFiberWarping;
            result = this->GetGMBoldsByFiber(pointIndex, closestFiberId, data->vertexTimeSeries[pointIndex], 0);
            if (-1 == result)
                result = this->GetGMBoldsByFiber(pointIndex, closestFiberId, data->vertexTimeSeries[pointIndex], 1);
            if (-1 == result)
                ++badFiberWarping;
        }

//             if(-1==result) //still no bolds. might be caused by bad
//             {
//                 result=this->GetBolds4Failures(pointIndex,m_data->vertexTimeSeries[pointIndex],1,closestFiberId);
//                 if(-1==result)
//                   result=this->GetBolds4Failures(pointIndex,m_data->vertexTimeSeries[pointIndex],2,closestFiberId);
//                 if(-1==result)
//                   result=this->GetBolds4Failures(pointIndex,m_data->vertexTimeSeries[pointIndex],3,closestFiberId);
//             }
        if (-1==result)
        {
            ++totalFailure;
            failureROIs.push_back(pointIndex);
        }

    } //end for all vertex;

    cout << "Bad Fiber Mapping: " << badFiberWarping / allFiberWarping << endl;
    cout << "Bad Normal Directionnn Mapping: " << badNormalWarping / allNormalWarping << endl;
    cout << "Total Mapping Failure: "<< totalFailure/data->numSurfacePoints <<endl;

    //final correction using neighbor time series;
//         int nbrSize =0;
//         int finalCorrected=0;
//         vector<bool> isFound(failureROIs.size(),false);
//         while(nbrSize++<3)
//         {
//           for(int idx = 0; idx < isFound.size(); ++idx )
//           {
//             if(false == isFound[idx])
//             {
//               list<int> nbrs;
//               m_data->surface.GetNeighboringPointsOfPoint(failureROIs[idx],nbrs,nbrSize);
//               for(list<int>::iterator it = nbrs.begin(); it!= nbrs.end(); ++it)
//               {
//                 if(m_data->vertexTimeSeries[*it][0]!=0 || m_data->vertexTimeSeries[*it][1]!=0)
//                 {
//                   m_data->vertexTimeSeries[failureROIs[idx]] = m_data->vertexTimeSeries[*it];
//                   ++finalCorrected;
//                   isFound[idx] = true;
//                   break;
//                 }
//               }
//             }
//           } // end for::idx;
//         }
//
//         cout << "CompleteFailure: "<< (totalFailure-finalCorrected)/m_data->numSurfacePoints <<endl;
}

void CJointModel::GeometricSmoothing()
{
    vector < vector<float> > oriSurfaceTimeSeries = data->vertexTimeSeries;
    for (int pointIndex = 0;pointIndex < data->numSurfacePoints;++pointIndex) {
        if (pointIndex % 2000 == 0)
            cout << "smoothing..." << pointIndex << endl;

        list<int> allNbrs;
        data->surface.GetNeighboringPointsOfPoint(pointIndex, allNbrs, data->smoothRange);
        // now we have all the nbrs within range;
        float weight = 1;
        const VectorType & srcPoint = data->surface.GetPointCoords(pointIndex);
        for (list<int>::iterator itAllNbrPoints = allNbrs.begin();itAllNbrPoints != allNbrs.end();++itAllNbrPoints) {
            int currentNbr = *itAllNbrPoints;
            const VectorType & dstPoint = data->surface.GetPointCoords(currentNbr);
            float dis = DistanceVector3D(srcPoint, dstPoint);
            float valueGaussian = GetGaussianKenerlValue(dis, data->sigma);
            weight += valueGaussian;
            data->vertexTimeSeries[pointIndex] += oriSurfaceTimeSeries[currentNbr] * valueGaussian;
        }
        data->vertexTimeSeries[pointIndex] *= 1. / weight;
    }

}

int CJointModel::GetNumSurfacePoints(void) const
{
    return data->numSurfacePoints;
}

int CJointModel::GetNumSurfaceCells(void) const
{
    return data->numSurfaceCells;
}

int CJointModel::GetNumSurfaceGyriPoints(void) const
{
    if (0 == data->numSurfaceGyriPoints) {
        for (int pointIndex = 0;pointIndex < data->numSurfacePoints;++pointIndex) {
            if (data->foldingPattern[pointIndex] == -1) {
                ++data->numSurfaceGyriPoints;
            }
        }

    }

    return data->numSurfaceGyriPoints;
}

int CJointModel::GetNumSurfaceSulciPoints(void) const
{
    if (data->numSurfaceSulciPoints == data->numSurfacePoints) {
        int gyriPoints = this->GetNumSurfaceGyriPoints();
        data->numSurfaceSulciPoints = data->numSurfacePoints - gyriPoints;
    }
    return data->numSurfaceSulciPoints;
}

int CJointModel::GetNumTimeSeries(void) const
{
    return data->numTimeSeries;
}

int CJointModel::GetNumFibers(void) const
{
    return data->numFibers;
}

int CJointModel::GetVolume2Del(void) const
{
    return data->vol2Del;
}

float CJointModel::GetAbsCorrelThresh(void) const
{
    return data->correlThresh;
}

float CJointModel::GetBoldDimX(void) const
{
    return data->xdimfMRI;
}

float CJointModel::GetBoldDimY(void) const
{
    return data->ydimfMRI;
}

float CJointModel::GetBoldDimZ(void) const
{
    return data->zdimfMRI;
}

Vector3D<float> & CJointModel::GetOffSet(void) const
{
    return data->offSet;
}

int CJointModel::GetNumSamples(void) const
{
    return data->numSamples;
}

vector<float> & CJointModel::GetVertexTimeSeris(int vertexId) const
{
    if (vertexId<0 || vertexId >= data->numSurfacePoints)
    {
        cerr<<"Invalid vertex id: "<<vertexId<<endl;
        cerr<<"Should be within range 0~"<<data->numSurfacePoints<<endl;
        exit(EXIT_FAILURE);
    }
    return data->vertexTimeSeries[vertexId];
}

CTriSurface & CJointModel::GetSurface(void)
{
    return data->surface;
}

void CJointModel::SetOffSet(Vector3D<float> & offset)
{
    data->offSet = offset;
}

void CJointModel::SetOffSet(float x, float y, float z)
{
    data->offSet.x = x;
    data->offSet.y = y;
    data->offSet.z = z;
}

void CJointModel::SetVolume2Del(int vol2del)
{
    if (vol2del<0 || vol2del> data->numTimeSeries)
    {
        cerr<<"Invalid number to ignore: "<<vol2del<<endl;
        exit(EXIT_FAILURE);
    }
    data->vol2Del = vol2del;
}

void CJointModel::SetAbsCorrelThresh(float thresh)
{
    if (thresh<0 || thresh>1)
    {
        cerr<<"Invalid correlation threshold: "<<thresh<<endl;
        exit(EXIT_FAILURE);
    }
    data->correlThresh = thresh;
}

void CJointModel::SetNumSamples(int num)
{
    if (num > data->numSurfacePoints) {
        cerr << "warning: too large sample number." << endl;
        cerr << "we will use all the surface points" << endl;
        num = data->numSurfacePoints - 1000;
    }
    if (GYRI == data->mode && num > data->numSurfaceGyriPoints) {
        cout << "we have only " << data->numSurfaceGyriPoints << " gyral vertices" << endl;
        cout << "all gyral vertices are selected as samples" << endl;
        num = data->numSurfaceGyriPoints;
        data->needSampling = false;
    }
    if (SULCI == data->mode && num > data->numSurfaceSulciPoints) {
        cout << "we have only " << data->numSurfaceSulciPoints << " sulcal vertices" << endl;
        cout << "all sulcal vertices are selected as samples" << endl;
        num = data->numSurfaceSulciPoints;
        data->needSampling = false;
    }
    data->numSamples = num;
}

void CJointModel::SetIndexer(CIndexer theIndexer) {
  if(data->indexer != theIndexer);
  {
    if(true==data->indexer.IsInitialized())
      cout<<"warning: indexer was reset!, please recall relavant functions, e.g.,CalculateFibersInsideCubic,FindSurfaceVolumeIntersection,BuildFiberEdgeIndiceWIndexer"<<endl;
    data->indexer = theIndexer; 
    data->offSet= data->indexer.GetOffset();
    data->isDoneCalculateFibersInsideCubic=false;
    data->isDoneFindSurfaceVolumeIntersection=false;
    data->isDoneBuildFiberEdgeIndiceWIndexer=false; 
    data->isDoneCalculateFibersInsideCubic=false; 
  }
} 

// void CJointModel::SetIndexer(CIndexer* pIndexer)
// {
//   if(data->indexer != *pIndexer);
//   {
//     if(true==data->indexer.IsInitialized())
//       cout<<"warning: indexer was reset!, please recall relavant functions, e.g.,CalculateFibersInsideCubic,FindSurfaceVolumeIntersection,BuildFiberEdgeIndiceWIndexer"<<endl;
//     data->indexer = *pIndexer; 
//     data->offSet= data->indexer.GetOffset();
//     data->isDoneCalculateFibersInsideCubic=false;
//     data->isDoneFindSurfaceVolumeIntersection=false;
//     data->isDoneBuildFiberEdgeIndiceWIndexer=false; 
//   }
// }

void CJointModel::IndexingIntersectionsOnCell(void) {

    for (int fid = 0; fid < data->numFibers; ++fid) {
        data->fiberIntersectionCellFront[fid].index= data->indexer.GetIndex( data->fiberIntersectionCellFront[fid].GetCoordFloat());
        data->fiberIntersectionCellRear[fid].index=  data->indexer.GetIndex( data->fiberIntersectionCellRear[fid].GetCoordFloat());
    } //end loop:: fid

}

int CJointModel::GetBolds4Failures(int pointId,vector<float> & timeSeries,int nbrSize, int fiberId,int firstPointsNum)
{

    if (fiberId != -1)
    {
        const Vector3D<float> pointCoord = data->surface.GetPointCoords(pointId);
        //determine which end to start;
        vector<int> & currentFiber = data->fibers.GetFiber(fiberId);
        Vector3D<float> & headPoint = data->fibers.GetPointCoords(currentFiber[0]);
        Vector3D<float> & endPoint = data->fibers.GetPointCoords(currentFiber[currentFiber.size() - 1]);
        Vector3D<float> startPoint;
        float headDis = AbsDistanceVector3D(headPoint, pointCoord);
        float endDis = AbsDistanceVector3D(endPoint, pointCoord);
        if (headDis < endDis) {
            startPoint = headPoint;
        } else {
            startPoint = endPoint;
        }
        int gmCount = 0;
        int fiberPointIndex = 4;

        while (fiberPointIndex < firstPointsNum) {
            Vector3D<int> startPointImgCoord;
            this->GetImageCoordDti(startPoint, startPointImgCoord);
            int boldx = startPointImgCoord.x * data->dtiSeg.xdim() / data->bolds.xdim() + 0.5;
            int boldy = startPointImgCoord.y * data->dtiSeg.ydim() / data->bolds.ydim() + 0.5;
            int boldz = startPointImgCoord.z * data->dtiSeg.zdim() / data->bolds.zdim() + 0.5;
            //check whether there is gm voxel in nbr hood;
            if (data->bolds(boldx, boldy, boldz,0)!= 0)
            {
                for (int t = 0;t < timeSeries.size();++t) {
                    timeSeries[t] = data->bolds(boldx, boldy, boldz, t + data->vol2Del);
                }
                return 1;
            }

            //update startpoint;
            ++fiberPointIndex;
            if (fiberPointIndex >= currentFiber.size())
                break;

            startPoint = data->fibers.GetPointCoords(currentFiber[fiberPointIndex]);
        } //end of while;

        //increase nbrsize;
        fiberPointIndex =4;
        while (fiberPointIndex < firstPointsNum) {
            Vector3D<int> startPointImgCoord;
            this->GetImageCoordDti(startPoint, startPointImgCoord);
            //check whether there is gm voxel in nbr hood;
            for (int zz = -nbrSize;zz < nbrSize + 1;++zz) {
                for (int yy = -nbrSize;yy < nbrSize + 1;++yy) {
                    for (int xx = -nbrSize;xx < nbrSize + 1;++xx) {
                        int nbrVoxelX = startPointImgCoord.x + xx;
                        int nbrVoxelY = startPointImgCoord.y + yy;
                        int nbrVoxelZ = startPointImgCoord.z + zz;
                        int boldx = nbrVoxelX * data->dtiSeg.xdim() / data->bolds.xdim() + 0.5;
                        int boldy = nbrVoxelY * data->dtiSeg.ydim() / data->bolds.ydim() + 0.5;
                        int boldz = nbrVoxelZ * data->dtiSeg.zdim() / data->bolds.zdim() + 0.5;
                        if (data->bolds(boldx, boldy, boldz,0)!= 0)//// gm voxel;
                        {
                            for (int t = 0;t < timeSeries.size();++t) {
                                timeSeries[t] = data->bolds(boldx, boldy, boldz, t + data->vol2Del);
                            }
                            return 1; //once we find one, quit.
                        }
                    }
                }
            } //end of check;

            //update startpoint;
            ++fiberPointIndex;
            if (fiberPointIndex >= currentFiber.size())
                break;
            startPoint = data->fibers.GetPointCoords(currentFiber[fiberPointIndex]);
        } //end of while;

    }
    else //using normal direction;
    {
        const Vector3D<float> & normal = data->surface.GetPointNormal(pointId);
        const Vector3D<float> & startPos = data->surface.GetPointCoords(pointId);
        Vector3D<float> currentPos = startPos;
        int step = 0;
        Vector3D<int> imageCoord;
        int gmCount = 0;
        float searchStep = 1;
        while (step < firstPointsNum) {
            Vector3D<float> nextPos = startPos + searchStep * normal; // for this pos, if there's no gm voxel within its nbr; go on searching; if has, then using this voxel's nbr average as the gm;
            this->GetImageCoordDti(nextPos, imageCoord);
            int boldx = imageCoord.x * data->dtiSeg.xdim() / data->bolds.xdim() + 0.5;
            int boldy = imageCoord.y * data->dtiSeg.ydim() / data->bolds.ydim() + 0.5;
            int boldz = imageCoord.z * data->dtiSeg.zdim() / data->bolds.zdim() + 0.5;
            if (data->bolds(boldx,boldy,boldz,0)!= 0)
            {
                for (int t = 0;t < timeSeries.size();++t) {
                    timeSeries[t] = data->bolds(boldx, boldy, boldz, t + data->vol2Del);
                }
                return 1;
            }
            //update step ;
            ++step;
        } // end of while;

        step = 0;
        while (step < firstPointsNum) {
            Vector3D<float> nextPos = startPos + searchStep * normal; // for this pos, if there's no gm voxel within its nbr; go on searching; if has, then using this voxel's nbr average as the gm;
            this->GetImageCoordDti(nextPos, imageCoord);
            //check whether there is gm voxel in nbr hood;
            for (int zz = -nbrSize;zz < nbrSize + 1;++zz) {
                for (int yy = -nbrSize;yy < nbrSize + 1;++yy) {
                    for (int xx = -nbrSize;xx < nbrSize + 1;++xx) {
                        int nbrVoxelX = imageCoord.x + xx;
                        int nbrVoxelY = imageCoord.y + yy;
                        int nbrVoxelZ = imageCoord.z + zz;
                        int boldx = nbrVoxelX * data->dtiSeg.xdim() / data->bolds.xdim() + 0.5;
                        int boldy = nbrVoxelY * data->dtiSeg.ydim() / data->bolds.ydim() + 0.5;
                        int boldz = nbrVoxelZ * data->dtiSeg.zdim() / data->bolds.zdim() + 0.5;
                        if ( data->bolds(boldx,boldy,boldz,0)!= 0 )//// gm voxel;
                        {
                            for (int t = 0;t < timeSeries.size();++t) {
                                timeSeries[t] = data->bolds(boldx, boldy, boldz, t + data->vol2Del);
                            }
                            return 1; //once we find one voxel, we quit.
                        }
                    }
                }
            }
            //update step ;
            ++step;
        } // end of while;
    } //end of else
    // if not return, use the current voxels bolds;
    Vector3D<float>  startPos = data->surface.GetPointCoords(pointId);
    Vector3D<int> imageCoord;
    this->GetImageCoordDti(startPos, imageCoord);
    for (int zz = -nbrSize;zz < nbrSize + 1;++zz) {
        for (int yy = -nbrSize;yy < nbrSize + 1;++yy) {
            for (int xx = -nbrSize;xx < nbrSize + 1;++xx) {
                int nbrVoxelX = imageCoord.x + xx;
                int nbrVoxelY = imageCoord.y + yy;
                int nbrVoxelZ = imageCoord.z + zz;
                int boldx = nbrVoxelX * data->dtiSeg.xdim() / data->bolds.xdim() + 0.5;
                int boldy = nbrVoxelY * data->dtiSeg.ydim() / data->bolds.ydim() + 0.5;
                int boldz = nbrVoxelZ * data->dtiSeg.zdim() / data->bolds.zdim() + 0.5;
                if (0!=data->bolds(boldx, boldy, boldz, 0) || 0!=data->bolds(boldx, boldy, boldz, 1) ||0!=data->bolds(boldx, boldy, boldz,3))
                {
                    for (int t = 0;t < timeSeries.size();++t) {
                        timeSeries[t] = data->bolds(boldx, boldy, boldz, t + data->vol2Del);
                    }
                    return 1; //once we find one voxel, we quit.
                }
            }
        }
    }

    return -1;
}
int CJointModel::GetGMBoldsByNormal(int pointId, int nbrSize, vector<float> & timeSeries, float searchStep, int maxStep)
{
    const Vector3D<float> & normal = data->surface.GetPointNormal(pointId);
    const Vector3D<float> & startPos = data->surface.GetPointCoords(pointId);
    Vector3D<float> currentPos = startPos;
    int step = 0;
    Vector3D<int> imageCoord;
    int gmCount = 0;
    if (0 == nbrSize) {
        while (step < maxStep) {
            Vector3D<float> nextPos = startPos + searchStep * normal; // for this pos, if there's no gm voxel within its nbr; go on searching; if has, then using this voxel's nbr average as the gm;
            this->GetImageCoordDti(nextPos, imageCoord);
            //check whether there is gm voxel in nbr hood;
            if (data->dtiSeg(imageCoord.x, imageCoord.y, imageCoord.z))//// gm voxel;
            {
                int boldx = imageCoord.x * data->dtiSeg.xdim() / data->bolds.xdim() + 0.5;
                int boldy = imageCoord.y * data->dtiSeg.ydim() / data->bolds.ydim() + 0.5;
                int boldz = imageCoord.z * data->dtiSeg.zdim() / data->bolds.zdim() + 0.5;
                if (data->bolds(boldx,boldy,boldz,0)!= 0)
                {
                    for (int t = 0;t < timeSeries.size();++t) {
                        timeSeries[t] = data->bolds(boldx, boldy, boldz, t + data->vol2Del);
                    }
                    return 1;
                }
            }

            //update step ;
            ++step;
        } // end of while;

    }
    else {
        while (step < maxStep) {
            Vector3D<float> nextPos = startPos + searchStep * normal; // for this pos, if there's no gm voxel within its nbr; go on searching; if has, then using this voxel's nbr average as the gm;
            this->GetImageCoordDti(nextPos, imageCoord);
            //check whether there is gm voxel in nbr hood;
            for (int zz = -nbrSize;zz < nbrSize + 1;++zz) {
                for (int yy = -nbrSize;yy < nbrSize + 1;++yy) {
                    for (int xx = -nbrSize;xx < nbrSize + 1;++xx) {
                        int nbrVoxelX = imageCoord.x + xx;
                        int nbrVoxelY = imageCoord.y + yy;
                        int nbrVoxelZ = imageCoord.z + zz;
                        int boldx = nbrVoxelX * data->dtiSeg.xdim() / data->bolds.xdim() + 0.5;
                        int boldy = nbrVoxelY * data->dtiSeg.ydim() / data->bolds.ydim() + 0.5;
                        int boldz = nbrVoxelZ * data->dtiSeg.zdim() / data->bolds.zdim() + 0.5;
                        if (data->dtiSeg(nbrVoxelX, nbrVoxelY, nbrVoxelZ)  && data->bolds(boldx,boldy,boldz,0)!= 0 )//// gm voxel;
                        {
                            for (int t = 0;t < timeSeries.size();++t) {
                                timeSeries[t] = data->bolds(boldx, boldy, boldz, t + data->vol2Del);
                            }
                            return 1; //once we find one voxel, we quit.
                        }

                    }

                }

            }

            //update step ;
            ++step;
        } // end of while;

    } //end of if else;

    return -1;
}

int CJointModel::GetGMBoldsByFiber(int pointId, int fiberId, vector<float> & timeSeries, int nbrSize, int firstPointsNum)
{
    const Vector3D<float> pointCoord = data->surface.GetPointCoords(pointId);
    //determine which end to start;
    vector<int> & currentFiber = data->fibers.GetFiber(fiberId);
    Vector3D<float> & headPoint = data->fibers.GetPointCoords(currentFiber[0]);
    Vector3D<float> & endPoint = data->fibers.GetPointCoords(currentFiber[currentFiber.size() - 1]);
    Vector3D<float> startPoint;
    float headDis = AbsDistanceVector3D(headPoint, pointCoord);
    float endDis = AbsDistanceVector3D(endPoint, pointCoord);
    if (headDis < endDis) {
        startPoint = headPoint;
    } else {
        startPoint = endPoint;
    }
    int gmCount = 0;
    int fiberPointIndex = 0;
    if (0 == nbrSize) {
        while (fiberPointIndex < firstPointsNum) {
            Vector3D<int> startPointImgCoord;
            this->GetImageCoordDti(startPoint, startPointImgCoord);
            //check whether there is gm voxel in nbr hood;
            if (data->dtiSeg(startPointImgCoord.x, startPointImgCoord.y, startPointImgCoord.z))//// gm voxel;
            {
                int boldx = startPointImgCoord.x * data->dtiSeg.xdim() / data->bolds.xdim() + 0.5;
                int boldy = startPointImgCoord.y * data->dtiSeg.ydim() / data->bolds.ydim() + 0.5;
                int boldz = startPointImgCoord.z * data->dtiSeg.zdim() / data->bolds.zdim() + 0.5;
                if (data->bolds(boldx, boldy, boldz,0)!= 0)
                {
                    for (int t = 0;t < timeSeries.size();++t) {
                        timeSeries[t] = data->bolds(boldx, boldy, boldz, t + data->vol2Del);
                    }
                    return 1;
                }
            }

            //update startpoint;
            ++fiberPointIndex;
            if (fiberPointIndex >= currentFiber.size())
                break;

            startPoint = data->fibers.GetPointCoords(currentFiber[fiberPointIndex]);
        } //end of while;
    }
    else {
        while (fiberPointIndex < firstPointsNum) {
            Vector3D<int> startPointImgCoord;
            this->GetImageCoordDti(startPoint, startPointImgCoord);
            //check whether there is gm voxel in nbr hood;
            for (int zz = -nbrSize;zz < nbrSize + 1;++zz) {
                for (int yy = -nbrSize;yy < nbrSize + 1;++yy) {
                    for (int xx = -nbrSize;xx < nbrSize + 1;++xx) {
                        int nbrVoxelX = startPointImgCoord.x + xx;
                        int nbrVoxelY = startPointImgCoord.y + yy;
                        int nbrVoxelZ = startPointImgCoord.z + zz;
                        int boldx = nbrVoxelX * data->dtiSeg.xdim() / data->bolds.xdim() + 0.5;
                        int boldy = nbrVoxelY * data->dtiSeg.ydim() / data->bolds.ydim() + 0.5;
                        int boldz = nbrVoxelZ * data->dtiSeg.zdim() / data->bolds.zdim() + 0.5;
                        if (data->dtiSeg(nbrVoxelX, nbrVoxelY, nbrVoxelZ) && data->bolds(boldx, boldy, boldz,0)!= 0)//// gm voxel;
                        {
                            for (int t = 0;t < timeSeries.size();++t) {
                                timeSeries[t] = data->bolds(boldx, boldy, boldz, t + data->vol2Del);
                            }
                            return 1; //once we find one, quit.
                        }

                    }

                }

            } //end of check;

            //update startpoint;
            ++fiberPointIndex;
            if (fiberPointIndex >= currentFiber.size())
                break;

            startPoint = data->fibers.GetPointCoords(currentFiber[fiberPointIndex]);
        } //end of while;
    }

    return -1;
}

void CJointModel::GetImageCoordDti(Vector3D<float> & phisicalCoord, Vector3D<int> & imageCoord)
{
    imageCoord.x = int((phisicalCoord.x - data->offSet.x) / data->xdimDti + 0.5);
    imageCoord.y = int((phisicalCoord.y - data->offSet.y) / data->ydimDti + 0.5);
    imageCoord.z = int((phisicalCoord.z - data->offSet.z) / data->zdimDti + 0.5);
}

int CJointModel::GetSmoothingRange(void) const
{
    return data->smoothRange;
}

float CJointModel::GetSmoothingSigma() const
{
    return data->sigma;
}

float CJointModel::GetGaussianKenerlValue(float x, float sigma)
{
    float val;
    val = exp(-(x * x) / (2.0 * sigma * sigma));
    return val;
}

void CJointModel::MakeUpMisClassification(string fileName)
{
    fstream cellsCingulateStrm;
    cellsCingulateStrm.open(fileName.c_str(), ios::in);
    vector<float> & gyriOrSulci = data->surface.GetPointDataScalar("SulciORGyri");
    while (!cellsCingulateStrm.eof()) {
        int cellId = 0;
        cellsCingulateStrm >> cellId;
        const CellType & currentPoints = data->surface.GetCell(cellId);
        gyriOrSulci[currentPoints.x] = -1; // -1 for gyri; and 1 for sulci;
        gyriOrSulci[currentPoints.y] = -1;
        gyriOrSulci[currentPoints.z] = -1;
        data->foldingPattern[currentPoints.x] = -1;
        data->foldingPattern[currentPoints.y] = -1;
        data->foldingPattern[currentPoints.z] = -1;
    }
    //recalculate the number of gyri and sulci;
    data->numSurfaceGyriPoints = 0;
    data->numSurfaceSulciPoints = 0;
    for (int pointIndex = 0;pointIndex < data->numSurfacePoints;++pointIndex) {
        if (1 == data->foldingPattern[pointIndex])
            ++data->numSurfaceSulciPoints;

        if (-1 == data->foldingPattern[pointIndex])
            ++data->numSurfaceGyriPoints;

    }
}

void CJointModel::GetCorrelation(list<int> & v1, list<int> & v2, list<float> & corr, vector<unsigned int> & index2SurfaceId)
{
    vector<bool> pointNotSelected(data->numSurfacePoints, true);
    vector<unsigned int> selectedSamples;
    int pointCount = 0;
    if (data->needSampling) {
        srand((unsigned int)(time(NULL)));
        while (pointCount<data->numSamples) {
            int newPointId= rand()%data->numSurfacePoints;
            if ( pointNotSelected[newPointId] && IsVertexRightMode(newPointId))// never got it, and in right mode; ;
            {
                pointNotSelected[newPointId]=false;
                selectedSamples.push_back(newPointId);
                ++pointCount;
            }
        }
    } else {
        for (int vertexId = 0 ; vertexId < data->numSurfacePoints; ++vertexId) {
            if (IsVertexRightMode(vertexId))
                selectedSamples.push_back(vertexId);
        } //end of for loop:: vertexId
    }
    assert(data->numSamples==selectedSamples.size());
    //          pointCount=0;
    //#pragma  omp parallel for default(shared) schedule(dynamic) num_threads(4)
    //          for (int pointIndex1=0; pointIndex1<m_data->numSamples; ++pointIndex1)
    //          {
    //#pragma omp critical
    //                  {
    //                          if (0==pointCount%500)
    //                          {
    //                                  cout<<pointCount<<" "<<corr.size()<<endl;
    //                          }
    //                          ++pointCount;
    //                  }
    //
    //                  int pointId1=selectedSamples[pointIndex1];
    //                  vector<float>& s1=m_data->vertexTimeSeries[pointId1];
    //
    //                  for (int pointIndex2=pointIndex1+1; pointIndex2<m_data->numSamples;++pointIndex2)
    //                  {
    //                          int pointId2= selectedSamples[pointIndex2];
    //                          vector<float>& s2=m_data->vertexTimeSeries[pointId2];
    //                          float absCorrel=AbsCoorelOfTwoCenteredSeries(s1,s2);
    //                          if (absCorrel> m_data->correlThresh)
    //                          {
    //#pragma omp critical
    //                                  {
    //                                          v1.push_back(pointIndex1);
    //                                          v2.push_back(pointIndex2);
    //                                          corr.push_back(absCorrel);
    //
    //                                  }
    //                          }
    //                  } //end of point 2 search;
    //
    //          } // all point1 search end here;
    int newPointId = 0;
    pointCount = 0;
    vector<int> idMapIndex(data->numSurfacePoints, -1);
    for (int pointIndex1 = 0;pointIndex1 < data->numSamples;++pointIndex1) {
        {
            if (0 == pointCount % 1000) {
                cout << pointCount << " " << corr.size() << endl;
            }
            ++pointCount;
        }
        vector<float> & s1 = data->vertexTimeSeries[selectedSamples[pointIndex1]];
        //cout<<s1<<endl;
        for (int pointIndex2 = pointIndex1 + 1;pointIndex2 < data->numSamples;++pointIndex2) {
            vector<float> & s2 = data->vertexTimeSeries[selectedSamples[pointIndex2]];
            //cout<<s2<<endl;
            float absCorrel = AbsCoorelOfTwoCenteredSeries(s1, s2, data->timeSeriesStart, data->timeSeriesEnd);
            //cout<<absCorrel<<endl;
            if (absCorrel > data->correlThresh) {
                {
                    int index1 = idMapIndex[selectedSamples[pointIndex1]];
                    int index2 = idMapIndex[selectedSamples[pointIndex2]];
                    if (-1 == index1) {
                        index1 = newPointId++;
                        idMapIndex[selectedSamples[pointIndex1]] = index1;
                        index2SurfaceId.push_back(selectedSamples[pointIndex1]);
                    }
                    if (-1 == index2) {
                        index2 = newPointId++;
                        idMapIndex[selectedSamples[pointIndex2]] = index2;
                        index2SurfaceId.push_back(selectedSamples[pointIndex2]);
                    }
                    v1.push_back(index1);
                    v2.push_back(index2);
                    corr.push_back(absCorrel);
                }
            }

        } //end of point 2 search;

    } // all point1 search end here;

}

//      void CJointModel::GetWholeCorrelation( list<int>& v1,list<int>& v2, list<float>& corr,vector<unsigned int>& selectedSamples)
//      {
//              vector<bool>  isPointAlreadySelected(m_data->numSurfacePoints,false);
//
//              int pointCount=0;
//              srand((unsigned int)(time(NULL)));
//              while(pointCount<m_data->numSamples){
//                      int newPointId= rand()%m_data->numSurfacePoints;
//                      if( false==isPointAlreadySelected[newPointId])// never got it;
//                      {
//                              isPointAlreadySelected[newPointId]=true;
//                              selectedSamples.push_back(newPointId);
//                              ++pointCount;
//                      }
//              }
//              assert(m_data->numSamples==selectedSamples.size());
//
//              pointCount=0;
//#pragma  omp parallel for default(shared) schedule(static) num_threads(4)
//              for (int pointIndex1=0; pointIndex1<m_data->numSamples; ++pointIndex1)
//              {
//#pragma omp critical
//                      {
//                              if (0==pointCount%500)
//                              {
//                                      cout<<pointCount<<" "<<corr.size()<<endl;
//                              }
//                              ++pointCount;
//                      }
//
//                      int pointId1=selectedSamples[pointIndex1];
//                      vector<float>& s1=m_data->vertexTimeSeries[pointId1];
//
//                      for (int pointIndex2=pointIndex1+1; pointIndex2<m_data->numSamples;++pointIndex2)
//                      {
//                              int pointId2= selectedSamples[pointIndex2];
//                              vector<float>& s2=m_data->vertexTimeSeries[pointId2];
//                              float absCorrel=AbsCoorelOfTwoCenteredSeries(s1,s2);
//                              if (absCorrel> m_data->correlThresh)
//                              {
//#pragma omp critical
//                                      {
//                                              v1.push_back(pointIndex1);
//                                              v2.push_back(pointIndex2);
//                                              corr.push_back(absCorrel);
//
//                                      }
//                              }
//                      } //end of point 2 search;
//
//              } // all point1 search end here;
//      }
void CJointModel::ClassifyUnselectedVertices(vector<unsigned int> & clusterIndexMapSurfacePointID, vector<unsigned int> & centersId)
{
    //cout<<clusterIndexMapSurfacePointID<<endl;
    cout << centersId << endl;
    vector<bool> isVertexSelected(data->numSurfacePoints, false);
    for (size_t vertexIndex = 0;vertexIndex < clusterIndexMapSurfacePointID.size();++vertexIndex) {
        isVertexSelected[clusterIndexMapSurfacePointID[vertexIndex]] = true;
    } //end of for loop:: vertexId
    vector<float> absSim2Centers(centersId.size(), 0);
    for (int vertexId = 0;vertexId < data->numSurfacePoints;++vertexId) {
        if (isVertexSelected[vertexId] || !IsVertexRightMode(vertexId))
            continue;

        //for unselected vertex, caculate the sim to centers, and return the max and max position;
        for (int centerIndex = 0;centerIndex < centersId.size();++centerIndex) {
            int tmpCenterid = centersId[centerIndex];
            float tmpAbsSim = AbsCoorelOfTwoCenteredSeries(data->vertexTimeSeries[vertexId], data->vertexTimeSeries[tmpCenterid], data->timeSeriesStart, data->timeSeriesEnd);
            absSim2Centers[centerIndex] = tmpAbsSim;
        } //end of for loop:: centerIndex
        float maxSim = -11111111;
        int maxPos = -1;
        findMaxValueAndIndex(absSim2Centers, maxSim, maxPos);
        data->surface.SetPointLabel(vertexId, centersId[maxPos]);
    } //end of for loop:: vertexId

}

void CJointModel::MergeCluters(vector<unsigned int> & centersId, vector<unsigned int> & newCenterIds, float simThreshold)
{
    bool stopFlag = false;
    vector<int> mapMergeIndex(centersId.size(), -1); // if not -1, the cluster will be merged ;
    vector<int> mergeOrder;
    while (false == stopFlag) {
        //find the correlations of the series of clusters;
        int maxOuterPos = 0;
        int maxInnerPos = 0;
        float maxCorrel = -100;
        for (int outerIt = 0;outerIt < centersId.size();++outerIt) {
            if (-1 == mapMergeIndex[outerIt]) {
                for (int innerIt = outerIt + 1;innerIt < centersId.size();++innerIt) {
                    if (-1 == mapMergeIndex[innerIt]) {
                        float absCorrel = AbsCoorelOfTwoCenteredSeries(data->vertexTimeSeries[centersId[outerIt]], data->vertexTimeSeries[centersId[innerIt]], data->timeSeriesStart, data->timeSeriesEnd);
                        if (absCorrel > maxCorrel) {
                            maxOuterPos = outerIt;
                            maxInnerPos = innerIt;
                            maxCorrel = absCorrel;
                        }
                    }

                } // end of for loop

            }

        } //end of for loop

        if (maxCorrel > simThreshold) {
            //this->RelabelClusteringRst(maxInnerPos+1,maxOuterPos+1); // repleace maxInner (small size) by outer; +1 for label;
            mapMergeIndex[maxInnerPos] = maxOuterPos;
            mergeOrder.push_back(maxInnerPos);
            cout << "merge " << centersId[maxInnerPos] << " to " << centersId[maxOuterPos] << " (Correl=" << maxCorrel << ")" << endl;
            // relabel the clusters;
        } else {
            stopFlag = true;
        }
    } //end of while;

    //relabel the surface;
    vector<int> mapId2Label(data->numSurfacePoints, -1);
    int label = 0;
    newCenterIds.push_back(-1); // label 0 has no correspondence
    for (int centerIndex = 0;centerIndex < centersId.size();++centerIndex) {
        int isReplaced = mapMergeIndex[centerIndex];
        if (-1 == isReplaced) {
            int centerId = centersId[centerIndex];
            mapId2Label[centerId] = label + 1;
            ++label;
            newCenterIds.push_back(centerId);
        }
    }

    data->numClusters = label;
    cout << "number of clusters after merging: " << label << endl;
    for (int mergeIndex = 0;mergeIndex < mergeOrder.size();++mergeIndex) {
        int killPos = mergeOrder[mergeIndex];
        int newClusterPos = mapMergeIndex[killPos];
        int oriCenterId = centersId[killPos];
        int newCenterId = centersId[newClusterPos];
        this->ReLabelClusters(oriCenterId, newCenterId);
    }
    for (int mergeIndex = 0;mergeIndex < mergeOrder.size();++mergeIndex) {
        int killPos = mergeOrder[mergeIndex];
        centersId[killPos] = -1;
    }
    for (int index = 0;index < mapId2Label.size();++index) {
        int finalLabel = mapId2Label[index];
        if (-1 != finalLabel)
            this->ReLabelClusters(index, finalLabel);

    } //end of for loop:: index
    cout << "New cluster Ids: " << endl;
    cout << newCenterIds << endl;
}

void CJointModel::ReLabelClusters(int tobeKill, int dst)
{
    int pointCount = 0;
    for (int vertexId = 0;vertexId < data->numSurfacePoints;++vertexId) {
        if (tobeKill == data->surface.GetPointLabel(vertexId)) {
            ++pointCount;
            data->surface.SetPointLabel(vertexId, dst);
        }
    } //end of for loop:: vertexId

    cout << pointCount << " vertices with label " << tobeKill << " were replaced by " << dst << endl;
}

void CJointModel::SmoothingByVoting(int iterN)
{
    int iter = 0;
    if ( ALL==data->mode )//whole surface;
    {
        while (iter <iterN) {
            cout<<"Smoothing iteration # :"<<++iter<<endl;
            vector<int> oriLabeledSurf=data->surface.GetPointDataLabel();

            for (int vertexId = 0 ; vertexId < data->numSurfacePoints; ++vertexId) {
                // get the nbr of current vertexId;
                list<int> nbrs;
                data->surface.GetNeighboringPointsOfPoint(vertexId,nbrs, data->smoothRange);

                vector<int> hist(data->numClusters+1, 0);
                for (list<int>::iterator itNbrs=nbrs.begin(); itNbrs!=nbrs.end(); itNbrs++)
                {
                    int nbrLabel= oriLabeledSurf[*itNbrs];
                    ++hist[nbrLabel];
                }
                int max=-10000;
                int maxPos=-1;
                findMaxValueAndIndex(hist,max,maxPos);
                assert(-1!=maxPos);
                data->surface.SetPointLabel(vertexId,maxPos);
            } //end of for loop:: vertexId

        }// end of while;

    }
    else
    {
        float mode=0;
        if (GYRI==data->mode)
            mode=-1;
        else
            mode=1;

        while (iter <iterN) {
            cout<<"Smoothing iteration # :"<<++iter<<endl;
            vector<int> oriLabeledSurf=data->surface.GetPointDataLabel();

            for (int vertexId = 0 ; vertexId < data->numSurfacePoints; ++vertexId) {

                if (data->foldingPattern[vertexId]!= float(mode)) // no interested in this vertex;
                    continue;
                // get the nbr of current vertexId;
                list<int> nbrs;
                data->surface.GetNeighboringPointsOfPoint(vertexId,nbrs,data->smoothRange);

                vector<int> hist(data->numClusters+1, 0);
                for (list<int>::iterator itNbrs=nbrs.begin(); itNbrs!=nbrs.end(); itNbrs++)
                {
                    if (data->foldingPattern[*itNbrs]!=float(mode)) // no interested in this vertex;
                        continue;
                    int nbrLabel= oriLabeledSurf[*itNbrs];
                    ++hist[nbrLabel];
                }
                int max=-10000;
                int maxPos=-1;
                findMaxValueAndIndex(hist,max,maxPos);
                assert(-1!=maxPos);
                data->surface.SetPointLabel(vertexId,maxPos);
            } //end of for loop:: vertexId

        }// end of while;
    }
}

void CJointModel::ColoringTopClusters(int num)
{
    CColorSchemeRGB myColor;
    if (num > myColor.GetMaxNumOfColors()) {
        cerr << "Warning: support maximum clusters : " << myColor.GetMaxNumOfColors() << endl;
        cerr << "I will select only the top  " << myColor.GetMaxNumOfColors() << " clusters." << endl;
        num = myColor.GetMaxNumOfColors();
    }
    if (num > data->numClusters) {
        cerr << "Not enough clusters. I will select all of the clusters; " << endl;
        num = data->numClusters;
    }
    //calculate the size of each clusters;
    vector<int> clusterSize(data->numClusters + 1, 0);
    vector<int> & surfaceLabel = data->surface.GetPointDataLabel();
    for (int vertexId = 0;vertexId < data->numSurfacePoints;++vertexId) {
        int label = surfaceLabel[vertexId];
        ++clusterSize[label];
    } //end of for loop:: vertexId
    clusterSize[0] = 0;
    string baseName;
    if (ALL == data->mode)
        baseName = "all";

    else
        if (GYRI == data->mode)
            baseName = "gyri";

        else
            baseName = "sulci";


    string clusterSizeName = baseName + "ClusterSize.txt";
    fstream clusterSizeStr;
    clusterSizeStr.open(clusterSizeName.c_str(), ios::out);
    for (int label = 1;label < clusterSize.size();++label) {
        if (clusterSize[label])
            clusterSizeStr << label << " " << clusterSize[label] << endl;

    }
    clusterSizeStr.close();
    for (int clusterIndex = 0;clusterIndex < num;++clusterIndex) {
        int maxValue = -10000;
        int maxPos = -1;
        findMaxValueAndIndex(clusterSize, maxValue, maxPos);
        Vector3D<float> currentColor = myColor.GetColorByIndex(clusterIndex);
        //color the vertex using this color;
        this->AssignColor2Cluster(currentColor, maxPos);
        clusterSize[maxPos] = -1;
    } //end of for loop:: clusterIndex
    // whether or not to set background;
}

void CJointModel::AssignColor2Cluster(Vector3D<float> & color, int clusterlabel)
{
    for (int vertexId = 0;vertexId < data->numSurfacePoints;++vertexId) {
        if (clusterlabel == data->surface.GetPointLabel(vertexId))
            data->surface.SetPointColor(vertexId, color);

    } //end of for loop:: vertexId
}

void CJointModel::SetBackGroundColor(Vector3D<float> & color)
{
    this->AssignColor2Cluster(color, 0);
}

void CJointModel::SetMode(int m)
{
    if (1 == m) {
        data->mode = SULCI;
        cout << "Find network on sulci vertices" << endl;
    } else
        if (-1 == m) {
            data->mode = GYRI;
            cout << "Find network on gyri vertices" << endl;
        } else {
            data->mode = ALL;
            cout << "Find network on whole surface vertices" << endl;
        }

}

CIndexer& CJointModel::GetIndexer(void) {
    return data->indexer;
}

void CJointModel::SetTimePointRange4SimCalc(int start, int end)
{
    assert(start<end);
    assert(start>=0);
    assert(end<=data->numTimeSeries);
    data->timeSeriesStart = start;
    data->timeSeriesEnd = end;
}

void CJointModel::SetSmoothingRange(int range)
{
    data->smoothRange = range;
}

void CJointModel::SetSmoothingSigma(float sigma)
{
    data->sigma = sigma;
}

bool CJointModel::IsVertexRightMode(int pointId) const
{
    if (ALL == data->mode || (-1 == data->foldingPattern[pointId] && GYRI == data->mode) || (1 == data->foldingPattern[pointId] && SULCI == data->mode))
        return true;

    return false;
}
void CJointModel::PruneSmallClusters(vector<unsigned int> & centersId, int sizeThresh)
{
    vector<int> clusterSize(centersId.size(), 0); // label start from 1;
    for (int index = 0;index < data->numSurfacePoints;++index) {
        ++clusterSize[data->surface.GetPointLabel(index)];
    } //end of for loop:: index
    cout << "clusterID   clusterSize" << endl;
    for (int i = 1;i < clusterSize.size();++i) {
        cout << i << " " << clusterSize[i] << endl;
    }
    clusterSize[0] = 1000000; // it will not be the min one;
    vector<bool> cluster2Kill(data->numClusters + 1, false); //false means not to kill;
    while (true) {
        int minValue = 1000000;
        int minPos = -1;
        findMinValueAndIndex(clusterSize, minValue, minPos);
        if (minValue > sizeThresh)
            break;

        cluster2Kill[minPos] = true;
        cout << "cluster " << minPos << " will be pruned, size: " << clusterSize[minPos] << endl;
        clusterSize[minPos] = 1000000;
    } // end of while;
    // for each vertex of the small clusters, relabel them;
    for (int clusterLabel = 1;clusterLabel < centersId.size();++clusterLabel) {
        if (cluster2Kill[clusterLabel]) {
            for (int vertexId = 0;vertexId < data->numSurfacePoints;++vertexId) {
                int currentLabel = data->surface.GetPointLabel(vertexId);
                if (currentLabel == clusterLabel) {
                    // relabel this vertex with the closest cluster ;
                    vector<float> correl(centersId.size(), -100);
                    for (int cid = 1;cid < centersId.size();++cid) {
                        if (cluster2Kill[cid])
                            continue;

                        float absCorrel = AbsCoorelOfTwoCenteredSeries(data->vertexTimeSeries[vertexId], data->vertexTimeSeries[centersId[cid]], data->timeSeriesStart, data->timeSeriesEnd);
                        correl[cid] = absCorrel;
                    } //end of for loop:: cid
                    float maxValue = -11111111;
                    int maxPos = -1;
                    findMaxValueAndIndex(correl, maxValue, maxPos);
                    // change current label to maxpos;
                    data->surface.SetPointLabel(vertexId, maxPos);
                }

            } //end of for loop:: vertexId

        }

    } //end of for loop:: clusterIndex

}

void CJointModel::RemoveSmallIsolatedRegion(int sizeThresh)
{
    vector<bool> isVisited(data->surface.GetNumOfPoints(), false);
    vector<int> & pointData = data->surface.GetAllPointLabel();
    for (int pointIndex = 0;pointIndex < data->surface.GetNumOfPoints();++pointIndex) {
        if (false == isVisited[pointIndex]) {
            map<int,list<int> > nbrLabels;
            list<int> groupMembers;
            data->surface.Propagate2All(pointIndex, groupMembers, nbrLabels, pointData);
            for (list<int>::iterator itGroupMembers = groupMembers.begin();itGroupMembers != groupMembers.end();itGroupMembers++) {
                isVisited[*itGroupMembers] = true;
            }
            if (groupMembers.size() < sizeThresh) {
                cout << "Delegate:" << pointIndex << " Label:" << pointData[pointIndex] << " GroupSize:" << groupMembers.size() << endl;
                cout << "Remove this region and Relabel it using neighbor label: ";
                size_t maxNbr = 0;
                int maxNbrLabel = 0;
                for (map<int,list<int> >::iterator itNbrLabels = nbrLabels.begin();itNbrLabels != nbrLabels.end();itNbrLabels++) {
                    if (itNbrLabels->second.size() > maxNbr) {
                        maxNbr = itNbrLabels->second.size();
                        maxNbrLabel = itNbrLabels->first;
                    }
                }

                if (0 == maxNbr) {
                    cout << "Warning! Isolated vertex: " << pointIndex << endl;
                    continue;
                }
                cout << maxNbrLabel << endl;
                cout << "Neighbors of this kind: " << maxNbr << endl;
                for (list<int>::iterator itGroupToRelabel = groupMembers.begin();itGroupToRelabel != groupMembers.end();itGroupToRelabel++) {
                    pointData[*itGroupToRelabel] = maxNbrLabel;
                }
            }

        }

    } //end of for;

}

vector<CIntersectionType> & CJointModel::GetFrontIntersectingCellOfFibers(void)
{
    if (! data->isDoneFiberSurfIntersection)
        this->FindFiberSurfIntersections();
    return data->fiberIntersectionCellFront;
}

vector<CIntersectionType> & CJointModel::GetRearIntersectingCellOfFibers(void)
{
    if (! data->isDoneFiberSurfIntersection)
        this->FindFiberSurfIntersections();
    return data->fiberIntersectionCellRear;
}

vector<vector<int> > & CJointModel::GetIntersectingFibersOfCells(void)
{
    if (! data->isDoneIntersectingFibers4EachCell)
        this->FindFibersOnCells();
    return data->cellIntersectingFibers;
}

CIntersectionType & CJointModel::GetFrontIntersectingCellOfFiber(int fiberId)
{
    if (! data->isDoneFiberSurfIntersection)
        this->FindFiberSurfIntersections();
    if (fiberId<0 || fiberId>data->numFibers-1)
    {
        cerr<<"Invalid fiber id: "<<fiberId<<endl;
        exit(EXIT_FAILURE);
    }
    return data->fiberIntersectionCellFront[fiberId];
}

CIntersectionType & CJointModel::GetRearIntersectingCellOfFiber(int fiberId)
{
    if (! data->isDoneFiberSurfIntersection)
        this->FindFiberSurfIntersections();
    if (fiberId<0 || fiberId>data->numFibers-1)
    {
        cerr<<"Invalid fiber id: "<<fiberId<<endl;
        exit(EXIT_FAILURE);
    }
    return data->fiberIntersectionCellRear[fiberId];
}

vector<int> & CJointModel::GetIntersectingFibersOfCell(int cellId)
{
    if (! data->isDoneIntersectingFibers4EachCell)
        this->FindFibersOnCells();

    if (cellId<0 || cellId> data->numSurfaceCells-1)
    {
        cerr<<"Invalid cell id: "<<cellId<<endl;
        exit(EXIT_FAILURE);
    }
    return data->cellIntersectingFibers[cellId];
}

CFibers & CJointModel::GetFibers(void)
{
    return data->fibers;
}

vector<float> & CJointModel::GetFoldingPattern(void)
{
    return data->foldingPattern;
}

volume4D<float> & CJointModel::GetBOLDs(void)
{
    return data->bolds;
}

volume<float> & CJointModel::GetGM(void)
{
    return data->dtiSeg;
}

void CJointModel::CalculateFiberDensity(std::vector<float> & dst, int range)
{
}

CJointModel::CJointModel(const KML::CJointModel& )
{
}

void CJointModel::operator =(const KML::CJointModel& )
{
}


void CJointModel::GetFibersInsideROI( volume<float>& roi, list<int>& resutls) {

    data->xdimDti = roi.xdim();
    data->ydimDti = roi.ydim();
    data->zdimDti= roi.zdim();

    for (int fiberIndex = 0; fiberIndex < data->numFibers; ++fiberIndex) {
        vector<int>& currentFiber = this->data->fibers.GetFiber(fiberIndex);
        for (int eleIndex = 0; eleIndex < currentFiber.size(); ++eleIndex) {
            Vector3D<float>& currentElem= data->fibers.GetPointCoords( currentFiber[eleIndex]);
            Vector3D<int> imgGrid;
            this->GetImageCoordDti(currentElem, imgGrid);
            if ( roi( imgGrid.x, imgGrid.y, imgGrid.z))
            {
                resutls.push_back(fiberIndex);
                break;
            }
        } //end loop:: eleIndex
    } //end loop:: fiberIndex
}

void CJointModel::GetFibersInsideROI( volume<float>& roi, list<int>& resutls , int label) {

    data->xdimDti = roi.xdim();
    data->ydimDti = roi.ydim();
    data->zdimDti= roi.zdim();

    for (int fiberIndex = 0; fiberIndex < data->numFibers; ++fiberIndex) {
        vector<int>& currentFiber = this->data->fibers.GetFiber(fiberIndex);
        for (int eleIndex = 0; eleIndex < currentFiber.size(); ++eleIndex) {
            Vector3D<float>& currentElem= data->fibers.GetPointCoords( currentFiber[eleIndex]);
            Vector3D<int> imgGrid;
            this->GetImageCoordDti(currentElem, imgGrid);
            if ( roi( imgGrid.x, imgGrid.y, imgGrid.z)==label)
            {
                resutls.push_back(fiberIndex);
                break;
            }
        } //end loop:: eleIndex
    } //end loop:: fiberIndex
}


void CJointModel::GetFibersInsideROI(volume<float>& roi, list<int>& fibersCondisered,  list<int>& results, int label) {

    data->xdimDti = roi.xdim();
    data->ydimDti = roi.ydim();
    data->zdimDti= roi.zdim();

    for (list<int>::iterator it= fibersCondisered.begin(); it!= fibersCondisered.end(); ++it) {
        int fiberId= *it;
        vector<int>& currentFiber = this->data->fibers.GetFiber(fiberId);
        for (int eleIndex = 0; eleIndex < currentFiber.size(); ++eleIndex) {
            Vector3D<float>& currentElem=data->fibers.GetPointCoords( currentFiber[eleIndex]);;
            Vector3D<int> imgGrid;
            this->GetImageCoordDti(currentElem, imgGrid);
            if ( roi( imgGrid.x, imgGrid.y, imgGrid.z) ==label) {
                results.push_back(fiberId);
                break;
            }
        } //end loop:: eleIndex
    } //end loop:: fiberIndex
}


void CJointModel::GetFibersInsideROI(volume<float>& roi, list<int>& fibersCondisered,  list<int>& results) {

    data->xdimDti = roi.xdim();
    data->ydimDti = roi.ydim();
    data->zdimDti= roi.zdim();

    for (list<int>::iterator it= fibersCondisered.begin(); it!= fibersCondisered.end(); ++it) {
        int fiberId= *it;
        vector<int>& currentFiber = this->data->fibers.GetFiber(fiberId);
        for (int eleIndex = 0; eleIndex < currentFiber.size(); ++eleIndex) {
            Vector3D<float>& currentElem=data->fibers.GetPointCoords( currentFiber[eleIndex]);;
            Vector3D<int> imgGrid;
            this->GetImageCoordDti(currentElem, imgGrid);
            if ( roi( imgGrid.x, imgGrid.y, imgGrid.z)) {
                results.push_back(fiberId);
                break;
            }
        } //end loop:: eleIndex
    } //end loop:: fiberIndex
}
void CJointModel::GetFibersByVerticeList(std::list< int >& vertices, std::vector< int >& result)
{
    vector<int> nbrCells;
    data->surface.GetNeighboringCellsOfPoints(vertices,nbrCells);

    for (int idx = 0; idx < nbrCells.size(); ++idx )
    {
        int cell= nbrCells[idx];
        vector<int>& fibersOnCell = data->cellIntersectingFibers[cell];
//           result.insert(result.end(),fibersOnCell.begin(),fibersOnCell.end());
        for (int idxFiber = 0; idxFiber < fibersOnCell.size(); ++idxFiber )
        {
            result.push_back(fibersOnCell[idxFiber]);
        } // end for::idxFiber;
    } // end for::idx;

}

void CJointModel::GetFibersBySidNeighborhood( int fid, int ring, set<int>& results) {
    set<int> nbrcells;
    data->surface.GetNeighboringCellsOfPoint(fid, nbrcells, ring);

    for ( set<int>::iterator it= nbrcells.begin(); it!= nbrcells.end(); ++it) {
        int cid=*it;
        vector<int>& currentFibers=data->cellIntersectingFibers[cid];
        for (int index = 0; index < currentFibers.size(); ++index) {
            results.insert(currentFibers[index]);
        } //end loop:: index
    }
}
void CJointModel::GetFibersBySidNeighborhood(int fid, int ring, std::vector< int >& results)
{
    vector<int> nbrcells;
    data->surface.GetNeighboringCellsOfPoint(fid, nbrcells, ring);

    vector<int> fiberStaus(data->numFibers,0);
    for (int idx = 0; idx < nbrcells.size(); ++idx )
    {
        vector<int>& currentFibers=data->cellIntersectingFibers[nbrcells[idx]];
        for (int index = 0; index < currentFibers.size(); ++index) {
            if (fiberStaus[currentFibers[index]]==0)
            {
                results.push_back(currentFibers[index]);
                fiberStaus[currentFibers[index]]=1;
            }
        } //end loop:: index
    } // end for::idx;
}

void CJointModel::GetDajiangFibersBySidNeighborhood(int fid, int ring, std::set< int >& results)
{
    if (false==data->indexer.IsInitialized())
    {
        cout<< "please specify indexer, i will quit"<<endl;
        exit(1);
    }

    this->CalculateFibersInsideCubics();
    list<int> nbrPoints;
    data->surface.GetNeighboringPointsOfPoint(fid,nbrPoints,ring);

    for (list<int>::iterator it = nbrPoints.begin(); it!= nbrPoints.end(); it++)
    {
        const Vector3D<float>& coord = data->surface.GetPointCoords(*it);
        const set<int>& frbIt=this->GetFibersInsideCubicByCenterCoord(coord);
        results.insert(frbIt.begin(),frbIt.end());
    }
}

void CJointModel::GetUnIntersectedFibers(string fileName) {
    list<int> allInvalidFibers;

    for (int fiberIndex = 0; fiberIndex < data->numFibers; ++fiberIndex) {
        if (! data->fiberIntersectionCellFront[fiberIndex].IsValid() || ! data->fiberIntersectionCellRear[fiberIndex].IsValid())
            allInvalidFibers.push_back(fiberIndex);
    } //end loop:: fiberIndex

    data->fibers.SaveFibers(allInvalidFibers, fileName);
}

void CJointModel::CalculateFibersInsideCubics(void) {

    if (! data->isDoneCalculateFibersInsideCubic) {

        cout<< "calulating the fibers inside cubics...";

        const Vector3D<size_t>& size= data->indexer.GetSize();
        data->fibersInsideCubics.resize(size.x+1);
        for (int x = 0; x <= size.x; ++x) {
            data->fibersInsideCubics[x].resize(size.y+1);
            for (int y = 0; y <= size.y; ++y) {
                data->fibersInsideCubics[x][y].resize(size.z+1);
            } //end loop:: y
        } //end loop:: x

        for (int fid = 0; fid < data->numFibers; ++fid) {
            vector<int>& currentFiber= data->fibers.GetFiber(fid);
            for (int currentElement = 0; currentElement < currentFiber.size(); ++currentElement) {
                Vector3D<float>& coord= data->fibers.GetPointCoords( currentFiber[currentElement]);
                if (coord.x < data->indexer.GetOffset().x || coord.y < data->indexer.GetOffset().y ||coord.z < data->indexer.GetOffset().z )
                {
                    //cout<<"warning: fiber element is beyond the outline of index. please check. fid:"<< fid<<" ,element:"<<currentElement<<endl;
                    continue;
                }
                Vector3D<size_t> index= data->indexer.GetIndex(coord);
                data->fibersInsideCubics[index.x][index.y][index.z].insert(fid);
            } //end loop:: currentElement
        } //end loop:: fid

        cout<<"done."<<endl;


        data->isDoneCalculateFibersInsideCubic=true;
    }
}

void CJointModel::SaveFibersInsideCubic( Vector3D<size_t>&  index, string fileName) {

    if (! data->isDoneCalculateFibersInsideCubic)
        this->CalculateFibersInsideCubics();

    set<int>& fibersSample= data->fibersInsideCubics[index.x][index.y][index.z];
    cout<<"fiber size inside cubic for index "<< index.x <<" "<<index.y<<" "<<index.z <<" :"  <<fibersSample.size()<<endl;
    data->fibers.SaveFibers( fibersSample, fileName);
}


const set<int>& CJointModel::GetFibersInsideCubeFront(int fid) {
    CIntersectionType& oneSample= data->fiberIntersectionCellFront[fid];
    const Vector3D<size_t>&  index = oneSample.GetIndex();
    return data->fibersInsideCubics[index.x][index.y][index.z];
}


const set<int>& CJointModel::GetFibersInsideCubicByCenterCoord(const Vector3D<float>& centerCoord) {
  this->CalculateFibersInsideCubics();
    Vector3D<size_t> index = data->indexer.GetIndex(centerCoord );
    return data->fibersInsideCubics[index.x][index.y][index.z];
}
const std::set< int >& CJointModel::GetFibersInsideCubicByIndexGrid(const KML::Vector3D<size_t>& grid)
{
    return data->fibersInsideCubics[grid.x][grid.y][grid.z];
}


void CJointModel::BuildFiberEdgeIndiceWIndexer(void )
{
  if(false==data->isDoneBuildFiberEdgeIndiceWIndexer){
    
    //set indexer first;
    if(false==data->indexer.IsInitialized())
    {
      cout<<"Please set the indexer first!"<<endl;
      exit(1);
    }
    
    //find surf vol x first; x for intersection; 
    this->FindSurfaceVolumeIntersection();
    
    //some cleaning work incase index has been built previously; 
    data->grids4EachFiber.clear();
    data->fibersInEdgeGrids.clear();
    data->verticesInEdgeGrids.clear();
    
    int graceSize = 5; // fiver more voxels than the indexer cubic; 
    data->fibersInEdgeGrids.resize(data->indexer.GetSize().x+graceSize);
    data->verticesInEdgeGrids.resize(data->indexer.GetSize().x+graceSize);

    for( int idx=0; idx < data->indexer.GetSize().x+graceSize; ++idx ){
        data->fibersInEdgeGrids[idx].resize(data->indexer.GetSize().y+graceSize); 
        data->verticesInEdgeGrids[idx].resize(data->indexer.GetSize().y+graceSize); 
        for( int idy=0; idy < data->indexer.GetSize().y+graceSize; ++idy ){
            data->fibersInEdgeGrids[idx][idy].resize(data->indexer.GetSize().z+graceSize); 
            data->verticesInEdgeGrids[idx][idy].resize(data->indexer.GetSize().z+graceSize); 
        }//end of for::idy
    }//end of for::idx
     
    //for each fiber, find the chain of grids; 
    int nFiberLines= data->numFibers; 
    for( int fid=0; fid < nFiberLines ; ++fid ){
        vector<Vector3D<size_t> > theChain;
        vector<int>& theLine= data->fibers.GetFiber(fid); 
        Vector3D<size_t> rearGrid; 
        for( int idxElem=0; idxElem < theLine.size(); ++idxElem ){
            int pid= theLine[idxElem]; 
            Vector3D<float>& theCoord = data->fibers.GetPointCoords(pid); 
            Vector3D<size_t> theGrid= data->indexer.GetIndex(theCoord);
            if(rearGrid!=theGrid)
            {
                theChain.push_back(theGrid);
                rearGrid=theGrid; 
            }
            if(data->SurfVolX(theGrid.x,theGrid.y,theGrid.z))
            {
               data->fibersInEdgeGrids[theGrid.x][theGrid.y][theGrid.z].insert(fid);  
            }
        }//end of for::idxElem
        
        data->grids4EachFiber.push_back(theChain);
    }//end of for::fid 
    
    
    //for each vertex, find the containing boundary voxel; 
    
    for(int idx = 0; idx < data->surface.GetNumOfPoints(); ++idx )
    {  
      //get the coord for current vertex; 
      //get the grid for the vertex; 
      //push back ; 
      const VectorType& coord = data->surface.GetPointCoords(idx);
      Vector3D<size_t> theGrid = data->indexer.GetIndex(coord);
      data->verticesInEdgeGrids[theGrid.x][theGrid.y][theGrid.z].push_back(idx);
    } // end for::idx;
   
    data->isDoneBuildFiberEdgeIndiceWIndexer=true; 
  }
}
void CJointModel::FindSurfaceVolumeIntersection(void )
{
  if(false==data->isDoneFindSurfaceVolumeIntersection)
  {
   
    data->SurfVolX.reinitialize(data->indexer.GetSize().x,data->indexer.GetSize().y,data->indexer.GetSize().z);
    data->SurfVolX=0;
    data->SurfVolX.setdims(data->indexer.GetDims().x,data->indexer.GetDims().y,data->indexer.GetDims().z);
    data->SurfVolX.setDisplayMaximumMinimum(1,0);
    
    
    for( int sid=0; sid < data->surface.GetNumOfPoints(); ++sid ){
        const Vector3D<float>& theCoord= data->surface.GetPointCoords(sid);
        Vector3D<size_t> theGrid=data->indexer.GetIndex(theCoord);
        data->SurfVolX(theGrid.x,theGrid.y,theGrid.z)=1; 
    }//end of for::sid
    
    int count=0; 
    for( int z=0; z < data->SurfVolX.zsize(); ++z ){
        for( int y=0; y < data->SurfVolX.ysize(); ++y ){
            for( int x=0; x < data->SurfVolX.xsize(); ++x ){
                if(data->SurfVolX(x,y,z))
                    ++count;
            }//end of for::x
        }//end of for::y
    }//end of for::z
    cout<<"number of intersection voxels: "<<count<<endl;
    data->isDoneFindSurfaceVolumeIntersection=true; 
  }  
}

bool CJointModel::IsAIntersectionGrid(Vector3D<size_t>& grid)
{
  this->FindSurfaceVolumeIntersection();
  if(data->SurfVolX(grid.x,grid.y,grid.z))
    return true;
  else
    return false;
}

void CJointModel::GetDajiangFibersBySidNeighborhood(int fid, int ring, std::vector< int >& results)
{
    if (false==data->indexer.IsInitialized())
    {
        cout<< "please specify indexer, i will quit"<<endl;
        exit(1);
    }
    
 
    this->CalculateFibersInsideCubics();
    list<int> nbrPoints;
    data->surface.GetNeighboringPointsOfPoint(fid,nbrPoints,ring);

    for (list<int>::iterator it = nbrPoints.begin(); it!= nbrPoints.end(); it++)
    {
 
        const Vector3D<float>& coord = data->surface.GetPointCoords(*it);
        const set<int>& frbIt=this->GetFibersInsideCubicByCenterCoord(coord);
        if(frbIt.size())
        {
          std::set< int >::iterator iter = frbIt.begin(); 
          while(iter!=frbIt.end())
          {
            results.push_back(*iter);
            iter++;
          }
        }
    }
}

void CJointModel::SaveVolSurfX(string fileName)
{
  NEWIMAGE::save_volume<>(data->SurfVolX,fileName); 
}

}//end of namespace KML;
