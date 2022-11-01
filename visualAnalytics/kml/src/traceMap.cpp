/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#include "kmPCA.h"
#include "traceMap.h"
#include "fibers.h"
#include "vtkIdList.h"
#include "kaimingCommon.h"
#include "algorithm"

namespace KML
{
class CPrivacyTraceMap
{
public:
    vector<float> feature;
    vector<int> bundle;
    vector<vector<Vector3D<float> > > allTracePoints;
    Vector3D<float> seed;
    vector<Vector3D<float> > polars;
    int minLength;
    int interval;
    float angleStep;
    float densityRadius;
    double nAllTracePoints;
    CPCA pca;
    CFibers* pFibers;
    bool isDoneSetSeed,isDoneSetBundle,isDoneSetFibers,isDoneDoMapCalculation;

    CPrivacyTraceMap()
    {
        pFibers=NULL;
        interval=8;
        angleStep= PAI/6;
        densityRadius=0.3;
        nAllTracePoints=0;
        isDoneDoMapCalculation=false;
        isDoneSetFibers=false;
        isDoneSetBundle=false;
        isDoneSetSeed=false;
    }
    ~CPrivacyTraceMap() {};
};

CTraceMap::CTraceMap()
{
    data = new CPrivacyTraceMap();
    this->CalculatePolars();
}

CTraceMap::CTraceMap(const KML::CTraceMap& oth)
{
    data = new CPrivacyTraceMap();
    this->operator=(oth);
//   this->CalculatePolars();
}
CTraceMap& CTraceMap::operator=(const KML::CTraceMap& oth)
{
    data->pFibers=oth.data->pFibers;
    data->interval=oth.data->interval;
    data->angleStep=oth.data->angleStep;
    data->densityRadius=oth.data->densityRadius;
    data->nAllTracePoints=oth.data->nAllTracePoints;
    data->seed=oth.data->seed;
    data->bundle=oth.data->bundle;
    data->feature=oth.data->feature;
    data->allTracePoints=oth.data->allTracePoints;
    data->polars=oth.data->polars;
    data->isDoneDoMapCalculation=oth.data->isDoneDoMapCalculation;
    data->isDoneSetFibers=oth.data->isDoneSetFibers;
    data->isDoneSetBundle=oth.data->isDoneSetBundle;
    data->isDoneSetSeed=oth.data->isDoneSetSeed;
}

CTraceMap::~CTraceMap()
{
    delete data;
}
void CTraceMap::Reinitialized(void )
{
    data->feature.clear();
    data->bundle.clear();
    data->allTracePoints.clear();
    data->nAllTracePoints=0;
}

const std::vector< float >& CTraceMap::GetTraceMapFeatures(void )
{
    this->DoTraceMap();
    return data->feature;
}
const vector<vector< KML::Vector3D< float > > >& CTraceMap::GetTraceMapsPoints(void )
{
    return data->allTracePoints;
}
void CTraceMap::SetBundle(const std::vector< int >& list)
{
    {
        this->Reinitialized();
        data->bundle=list;
        data->isDoneSetBundle=true;
        data->isDoneDoMapCalculation=false;
    }
}

void CTraceMap::SetBundle(const std::list< int >& list)
{
    {
        this->Reinitialized();
        data->bundle.resize(list.size());
        copy(list.begin(),list.end(),data->bundle.begin());
        data->isDoneSetBundle=true;
        data->isDoneDoMapCalculation=false;
    }
}

void CTraceMap::SetBundle(const std::set< int >& list)
{
    {
        this->Reinitialized();
        data->bundle.resize(list.size());
        copy(list.begin(),list.end(),data->bundle.begin());
        data->isDoneSetBundle=true;
        data->isDoneDoMapCalculation=false;
    }
}

void CTraceMap::SetBundle(vtkSmartPointer<vtkIdList>& list)
{
    {
        this->Reinitialized();
        for ( int idx=0; idx < list->GetNumberOfIds() ; ++idx ) {
            data->bundle.push_back(list->GetId(idx));
        }//end of for::idx
        data->isDoneSetBundle=true;
        data->isDoneDoMapCalculation=false;
    }
}

void CTraceMap::SetFiber(CFibers* theFiber)
{
    if (NULL!=theFiber)
    {
        data->pFibers=theFiber;
        data->isDoneSetFibers = true;
    }
    else
    {
        cerr<<"Invalid Fiber Pointer"<<endl;
        exit(1);
    }
}

void KML::CTraceMap::SetSeedPoint(const KML::Vector3D< float >& seed)
{
    if (false==data->isDoneSetSeed || seed!=data->seed)
    {
        data->seed = seed;
        data->isDoneSetSeed=true;
        data->isDoneDoMapCalculation=false;
    }
}

void CTraceMap::DoTraceMap(void )
{
    if (false ==data->isDoneDoMapCalculation)
    {
        //for each fibe in the list;
        //calculate its trace points;
        // after the loop
        // calculate the features;

        if ( data->isDoneSetBundle && data->isDoneSetFibers && data->isDoneSetSeed)
        {
            for ( int index=0; index < data->bundle.size() ; ++index ) {
                int fiberID= data->bundle[index];
                vector<Vector3D<float> > theTrace4CurrentFiber;
                this->TracingSingleFiber(fiberID,theTrace4CurrentFiber);
                if (theTrace4CurrentFiber.size())
                {
                    data->nAllTracePoints+=theTrace4CurrentFiber.size();
                    data->allTracePoints.push_back(theTrace4CurrentFiber);
                }

            }//end of for::index
            this->CalculateFeatures();
            //     cout<<"done trace map"<<endl;
            data->isDoneDoMapCalculation=true;
        }
        else
        {
            cout<<"Something not set. I will quit calculation of tracemaps."<<endl;
            exit(1);
        }
    }
}

int CTraceMap::TracingSingleFiber(int id, vector<Vector3D<float> >& result)
{
    ///
    vector<int>& currentFiber = data->pFibers->GetFiber(id);
//   cout<<"id size:"<<id<<" "<<currentFiber.size()<<endl;
    if (currentFiber.size() < data->interval)
        return 0;

    //determine the fiber orientation;
    // if not start from the seed side; reverse the fiber;
    Vector3D<float>& headPoint= data->pFibers->GetPointCoords(currentFiber[0]);
    Vector3D<float>& rearPoint= data->pFibers->GetPointCoords(currentFiber[currentFiber.size()-1]);
    float dis1 = DistanceVector3D<float>(headPoint,data->seed);
    float dis2 = DistanceVector3D<float>(rearPoint,data->seed);
    if (dis1>dis2)
    {
        std::reverse(currentFiber.begin(), currentFiber.end());
    }

    for ( int j = data->interval; j < ( currentFiber.size()-data->interval) ; j +=data->interval ) {

        int tmpStart=j-data->interval;
        int tmpEnd=j+data->interval;
        int segSize=tmpEnd-tmpStart;

        Vector3D<float>& startPoint=data->pFibers->GetPointCoords(currentFiber[tmpStart]);
        Vector3D<float>& endPoint=data->pFibers->GetPointCoords(currentFiber[tmpEnd]);
        Vector3D<float> refDirection=endPoint-startPoint;
        refDirection.Normalize();

        vector<vector<double> > rawData(segSize, vector<double>(3,0));
        for ( int segIndex=0; segIndex < segSize; ++segIndex ) {
            int id=currentFiber[ segIndex+tmpStart];
            rawData[segIndex][0]= data->pFibers->GetPointCoords(id).x;
            rawData[segIndex][1]= data->pFibers->GetPointCoords(id).y;
            rawData[segIndex][2]= data->pFibers->GetPointCoords(id).z;
        }//end of for::segIndex


        data->pca.SetData(rawData);
        vector<vector<double> > basis;
        data->pca.DoPCA(basis);

        Vector3D<float> tracePointCurrentSeg (basis[0][0],basis[1][0],basis[2][0]);
        if (DotProduct(refDirection,tracePointCurrentSeg)<0)
        {
            tracePointCurrentSeg*=-1;
//        cout<<"reverse..."<<endl;
        }
        result.push_back(tracePointCurrentSeg);

    }//end of for::j

    return 1;

}

void CTraceMap::CalculatePolars()
{
    int nStep= int(2*PAI/data->angleStep+0.5);

    double angleTheta = 0.0;
    double anglePhi = 0.0;
    data->polars.clear();
    for (int i = 0; i < nStep; i++) {
        angleTheta = data->angleStep * i;
        for (int j = 0; j < nStep; j++) {
            anglePhi = data->angleStep * j;
            Vector3D<float> samplePoint;
            samplePoint.x = (float) (cos(angleTheta) * cos(anglePhi));
            samplePoint.y = (float) (sin(angleTheta) * cos(anglePhi));
            samplePoint.z = (float) (sin(anglePhi));
            data->polars.push_back(samplePoint);
        }
    }
}
void CTraceMap::CalculateFeatures(void )
{

//     int nStep= int(2*PAI/data->angleStep+0.5);
    data->feature.clear();
    data->feature.resize(data->polars.size(),0);

    if (data->nAllTracePoints)
    {
        for (int idxPolar = 0; idxPolar < data->polars.size(); ++idxPolar )
        {
            int pointCount=0;
            for ( int idxLine=0; idxLine < data->allTracePoints.size() ; ++idxLine ) {
                for ( int idxCol=0; idxCol < data->allTracePoints[idxLine].size() ; ++idxCol ) {
                    Vector3D<float>& theTracePoint= data->allTracePoints[idxLine][idxCol];
                    Vector3D<float> diff=data->polars[idxPolar]-theTracePoint;
                    if (diff.Norm()< data->densityRadius)
                        ++pointCount;
                }//end of for::idxCol
            }//end of for::idxLine

            float theDensity=pointCount/data->nAllTracePoints;
            data->feature[idxPolar] = theDensity;
        } // end for::idxPolar;
        
        KML::StatsNormalizeHist(data->feature);
        /*      float tmpSum = 0.0f;
                double angleTheta = 0.0;
                double anglePhi = 0.0;
                for (int i = 0; i < nStep; i++) {
                  angleTheta = data->angleStep * i;
                  for (int j = 0; j < nStep; j++) {
                      anglePhi = data->angleStep * j;

                      Vector3D<float> samplePoint;
                      samplePoint.x = (float) (cos(angleTheta) * cos(anglePhi));
                      samplePoint.y = (float) (sin(angleTheta) * cos(anglePhi));
                      samplePoint.z = (float) (sin(anglePhi));
                      Vector3D<float> diff;

                      int pointCount=0;


                      for( int idxLine=0; idxLine < data->allTracePoints.size() ; ++idxLine ){
                        for( int idxCol=0; idxCol < data->allTracePoints[idxLine].size() ; ++idxCol ){
                          Vector3D<float>& theTracePoint= data->allTracePoints[idxLine][idxCol];
                          diff=samplePoint-theTracePoint;
                          if(diff.Norm()< data->densityRadius)
                            ++pointCount;
                        }//end of for::idxCol
                      }//end of for::idxLine

                      float theDensity=pointCount/data->nAllTracePoints;
                      data->feature.push_back(theDensity);
                  }
                }
              */
    }
}

void CTraceMap::SetAngleStep(float step)
{
    if (step> PAI )
    {
        cout<<"warning: will transfer to ringe [0 Pai]"<<endl;
        data->angleStep=step*PAI/180;
    }
    data->angleStep=step;
    this->CalculatePolars();
}

void CTraceMap::SetSegLength(int length)
{
    data->interval=length;
}

void CTraceMap::SetDensityRadius(float radius)
{
    data->densityRadius=radius;
}



}//endof KML;
