#include "Vector3D.h"
#include "kaimingCommon.h"
#include "indexer.h"
using namespace KML;
namespace KML{
class CIndexerPrivacy{
	friend class CIndexer;
public:
	CIndexerPrivacy(){};
	Vector3D<size_t> gridSize;
	Vector3D<float> dims;
	Vector3D<float> offSets;
        size_t zOrder; 
        size_t yOrder; 
        bool isInitialzed; 
};

CIndexer::CIndexer()
{
	data= new CIndexerPrivacy;
	data->dims=1;
	data->offSets=1;
	data->gridSize=1;
        data->zOrder=0; 
        data->yOrder=0; 
        data->isInitialzed=false; 
}
CIndexer::CIndexer(const char* mhdFile)
{
	data= new CIndexerPrivacy;
	KML::GetDimsFromMHD(mhdFile, data->dims);
	KML::GetOffsetFromMHD(mhdFile, data->offSets);
	KML::GetSizesFromMHD(mhdFile,data->gridSize);

        data->zOrder = data->gridSize.y*data->gridSize.x; 
        data->yOrder = data->gridSize.x; 
	cout<<"index size: "<< data->gridSize<<endl;
	cout<<"index dims: "<<data->dims<<endl;
        data->isInitialzed=true; 

}
CIndexer::CIndexer(string mhdFile)
{
        data= new CIndexerPrivacy;
        KML::GetDimsFromMHD(mhdFile.c_str(), data->dims);
        KML::GetOffsetFromMHD(mhdFile.c_str(), data->offSets);
        KML::GetSizesFromMHD(mhdFile.c_str(),data->gridSize);

        data->zOrder = data->gridSize.y*data->gridSize.x; 
        data->yOrder = data->gridSize.x; 
        cout<<"index size: "<< data->gridSize<<endl;
        cout<<"index dims: "<<data->dims<<endl;
        data->isInitialzed=true; 

}
CIndexer::~CIndexer()
{
	delete data;
}

void CIndexer::SetDims(const Vector3D<float>& dims){
	data->dims= dims;
}

void CIndexer::SetOffset( const Vector3D<float>& offSets)
{
	data->offSets= offSets;
}

void CIndexer::SetSize( const Vector3D<size_t>& size){
	data->gridSize= size;
        data->zOrder = data->gridSize.y*data->gridSize.x; 
        data->yOrder = data->gridSize.x; 
}

const Vector3D<size_t>& CIndexer::GetSize(void) const {
	return data->gridSize;
}

const Vector3D<float>& CIndexer::GetOffset(void) const{
	return data->offSets;
}

const Vector3D<float>& CIndexer::GetDims(void) const{
	return data->dims;
}

Vector3D<size_t> CIndexer::GetIndex(const Vector3D<float>& pos)const{
	Vector3D<size_t> theGridIndex;
//         theGridIndex.x= size_t((pos.x-data->offSets.x)/ data->dims.x /*+0.5*/); // get the floor ; 
// 	   theGridIndex.y= size_t((pos.y-data->offSets.y)/ data->dims.y /*+0.5*/);
//         theGridIndex.z= size_t((pos.z-data->offSets.z)/ data->dims.z /*+0.5*/);
        
        theGridIndex.x = floor((pos.x-data->offSets.x)/ data->dims.x);
        theGridIndex.y = floor((pos.y-data->offSets.y)/ data->dims.y);
        theGridIndex.z = floor((pos.z-data->offSets.z)/ data->dims.z);
        
        if(pos.x<=data->offSets.x)
          theGridIndex.x=0;
        if(pos.y<=data->offSets.y)
          theGridIndex.y=0;
        if(pos.z<=data->offSets.z)
          theGridIndex.z=0; 
        theGridIndex.x= theGridIndex.x >= data->gridSize.x ? data->gridSize.x-1 : theGridIndex.x;
        theGridIndex.y= theGridIndex.y >= data->gridSize.y ? data->gridSize.y-1 : theGridIndex.y; 
        theGridIndex.z= theGridIndex.z >= data->gridSize.z ? data->gridSize.z-1 : theGridIndex.z; 
        
//         if(theGridIndex.x >= data->gridSize.x  || theGridIndex.y >= data->gridSize.y || theGridIndex.z >= data->gridSize.z )
//         {
//           cout<<"warning: invalid index."<<theGridIndex<<endl;
//           theGridIndex.x=theGridIndex.x >= data->gridSize.x ? data->gridSize.x-1 : theGridIndex.x;
//           theGridIndex.y=theGridIndex.y >= data->gridSize.y ? data->gridSize.y-1 : theGridIndex.y; 
//           theGridIndex.z=theGridIndex.z >= data->gridSize.z ? data->gridSize.z-1 : theGridIndex.z ;
//         }
//         if(theGridIndex.x<0 || theGridIndex.y<0 || theGridIndex.z<0 )
//         {
//           cout<<"warning: invalid index."<<theGridIndex<<endl;
//           theGridIndex.x=theGridIndex.x <0 ? 0 : theGridIndex.x;
//           theGridIndex.y=theGridIndex.y <0 ? 0 : theGridIndex.y; 
//           theGridIndex.z=theGridIndex.z <0 ? 0 : theGridIndex.z ;
//         }
	return theGridIndex;
}

Vector3D<size_t> CIndexer::GetIndex(float pos[3])const{
	Vector3D<float> newPos ( pos[0], pos[1], pos[2]);
	return this->GetIndex( newPos);
}

Vector3D<size_t> CIndexer::GetIndex(double pos[3])const{
	Vector3D<float> newPos ( pos[0], pos[1], pos[2]);
	return this->GetIndex( newPos);
}

void CIndexer::ResetDims(const Vector3D<float>& newDims){
	Vector3D<size_t> newSize;
	newSize.x = (size_t) ( 0.5+ data->gridSize.x* data->dims.x/ newDims.x);
	newSize.y= (size_t) ( 0.5+ data->gridSize.y* data->dims.y/ newDims.y);
	newSize.z= (size_t) ( 0.5+ data->gridSize.z* data->dims.z/ newDims.z);
	data->gridSize= newSize;
	data->dims= newDims;
        data->zOrder = data->gridSize.y*data->gridSize.x; 
        data->yOrder = data->gridSize.x; 
// 	cout<<"new size: "<< data->gridSize<<endl;
// 	cout<<"new dims: "<<data->dims<<endl;
}

void CIndexer::ResetDims(float newRes)
{
  Vector3D<float> newDims(newRes,newRes,newRes);
  this->ResetDims(newDims);
}


bool CIndexer::operator == (const KML::CIndexer& oth)const{
  
  return ( data->gridSize == oth.data->gridSize) && (data->dims==oth.data->dims)  && (data->offSets==oth.data->offSets);
  
}

CIndexer& CIndexer::operator = (const CIndexer& oth){
  if(this == &oth)
    return *this; 
  
  data->gridSize = oth.data->gridSize;
  data->dims= oth.data->dims; 
  data->offSets= oth.data->offSets; 
  data->zOrder= oth.data->zOrder; 
  data->yOrder= oth.data->yOrder; 
  data->isInitialzed=oth.data->isInitialzed;
  
  return *this; 
}

CIndexer::CIndexer(const KML::CIndexer& oth)
{
  data= new CIndexerPrivacy;
  data->gridSize = oth.data->gridSize;
  data->dims= oth.data->dims; 
  data->offSets= oth.data->offSets; 
  data->zOrder= oth.data->zOrder; 
  data->yOrder= oth.data->yOrder; 
  data->isInitialzed=oth.data->isInitialzed;
    
}

size_t CIndexer::GetHistIndex(Vector3D< ::size_t >& idx)const
{
    return idx.z*data->zOrder+ idx.y*data->yOrder+idx.x; 
}


size_t CIndexer::GetNumGrids(void)const
{
    return data->gridSize.x*data->gridSize.y*data->gridSize.z;
}

size_t CIndexer::GetHistIndex(const KML::Vector3D< float >& pos)const
{
  Vector3D<size_t> index=this->GetIndex(pos);
  return this->GetHistIndex(index); 
}

bool CIndexer::IsInitialized(void )const
{
  return data->isInitialzed;
}

bool CIndexer::operator!=(const KML::CIndexer& oth)const
{
  return !(this->operator==(oth)); 
}
KML::Vector3D< float > CIndexer::GetWorldCood(const KML::Vector3D< int >& grid)
{

  Vector3D<float> coord = this->GetWorldCood(Vector3D<size_t> (grid.x,grid.y,grid.z));  
  return coord; 
}
KML::Vector3D< float > CIndexer::GetWorldCood(const KML::Vector3D< size_t>& grid)
{
  
  Vector3D< float > coord=data->offSets;
  coord.x +=(grid.x+0.5)*data->dims.x; 
  coord.y +=(grid.y+0.5)*data->dims.y; 
  coord.z +=(grid.z+0.5)*data->dims.z; 
  return coord; 

}
  void CIndexer::OutputIndexerAppend(string fileName)const
  {
    fstream outStream;
    outStream.open(fileName.c_str(),ios::out|ios::app);
    if(NULL==outStream)
    {
      cout<<"can not open output file :"<<fileName<<endl;
      exit(1);
    }
    
    outStream<<"################################"<<endl;
    outStream<<"dims: "<<data->dims<<endl;
    outStream<<"sizes: "<<data->gridSize<<endl;
    outStream<<"offsets: "<<data->offSets<<endl;
    outStream.close();    
  }
  
}// end of KM:


