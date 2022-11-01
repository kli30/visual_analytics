/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#include "fibers.h"
#include <iostream>
#include "jointModel.h"
#include "tictoc.hpp"
#include <list> 

using namespace std;
using namespace KML;

int main(int argc, char** argv)
{
  if( 5!=argc)
  {
    cout<<"usage: "<<argv[0]<<" surf fiber mhd sid"<<endl;
    exit(1);
  }
  CJointModel myMdl(argv[1],argv[2]); 
  CIndexer myIdx(argv[3]); 
  Vector3D<float> newDims;
  newDims=2.5; 
  myIdx.ResetDims(newDims);

  myMdl.FindFiberSurfIntersections();
  myMdl.SaveALLIntersectionPoints("intersection.vtk");
  myMdl.FindFibersOnCells();
  
  myMdl.SetIndexer(myIdx);
  myMdl.CalculateFibersInsideCubics();
  
  set<int> bundles;
  int id = atoi(argv[4]);
  tictoc timer;
  timer.tic();
  myMdl.GetFibersBySidNeighborhood(id,3,bundles);
  timer.toc();
  
  cout<<"my time :"<<timer.totalTimeSec()<<endl;
  timer.clear();;
  CTriSurface& mySurf = myMdl.GetSurface(); 
  CFibers& myFiber = myMdl.GetFibers(); 
  myFiber.SaveFibers(bundles,"myBundle.vtk");  
  cout<<"number of fibers: "<<bundles.size()<<endl;
  bundles.clear();
  timer.tic();
  myMdl.GetDajiangFibersBySidNeighborhood(id,3,bundles);
  timer.toc();
  cout<<"dajiang time: "<<timer.totalTimeSec()<<endl;
  myFiber.SaveFibers(bundles,"djBundle.vtk");
  cout<<"number of fibers: "<<bundles.size()<<endl;
  
  
}
