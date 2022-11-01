/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#ifndef ___TraceMap__H
#define ___TraceMap__H
#include <vector> 
#include <Vector3D.h>
#include <list>
#include "vtkSmartPointer.h"
#include "vtkIdList.h"
#include <set> 

using namespace std; 

namespace KML{
class CPrivacyTraceMap; 
class CFibers; 

class CTraceMap
{
private:
  CPrivacyTraceMap* data; 
  int TracingSingleFiber(int id, vector< KML::Vector3D< float > >& result); 
  void CalculateFeatures(void);
  void Reinitialized(void);
  void CalculatePolars(void); 
public:
  CTraceMap();
  CTraceMap(const CTraceMap& oth); //copy structor, for push_back or =; 
  ~CTraceMap(); 
  void SetFiber(CFibers* theFiber); 
  void SetBundle(const std::vector< int >& list);  
  void SetBundle(const std::set< int >& list);  
  void SetBundle(const std::list< int >& list); 
  void SetBundle(vtkSmartPointer<vtkIdList>& list); 
  void SetAngleStep(float step=30);
  void SetSegLength(int length=8);
  void SetDensityRadius(float radius=0.3);
  void SetSeedPoint(const Vector3D<float>& seed); 
  const std::vector< vector< KML::Vector3D< float > > >& GetTraceMapsPoints(void); 
  const std::vector< float >& GetTraceMapFeatures(void);  
  //deprecated, no need to call this function now; will be declared as private soon; 
  void DoTraceMap(void); 
  CTraceMap& operator=(const CTraceMap& oth);
};

}; //end of namespace KML;


#endif
