/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/


#include "roidrag.h"
#include "vtkPropPicker.h"
#include "vtkPointPicker.h"



void ROIDrag::MyMouseMove(void)
{
  vtkPropPicker* pker=  static_cast< vtkPropPicker*> ( iren->GetPicker());
  
  if(pker->GetActor()&& iren->GetControlKey()){
	pker->GetActor()->SetPickable(0);
	iren->GetMousePosition(&cx,&cy);
	vtkPointPicker* pPicker= vtkPointPicker::New(); 
	iren->SetPicker(pPicker);

	if(pPicker->Pick(cx,cy,0,ren)){
	    worldCoordinate=pPicker->GetPickPosition();   
//	    cout<<worldCoordinate[0]<<" "<<worldCoordinate[1]<<" "<<worldCoordinate[2]<<endl;
	    roiSphere->SetCenter(worldCoordinate);
	    iren->GetRenderWindow()->Render();
	}
	pPicker->Delete();
	iren->SetPicker(pker);
  }
  else
  {
    cout<<"not moving..."<<endl;
  }
 
}
void ROIDrag::SetInteractor(vtkRenderWindowInteractor* interactor)
{
  iren= interactor; 
}
void ROIDrag::SetRenderer(vtkRenderer* theRen)
{
  ren= theRen; 
}
void ROIDrag::SetROISource(vtkSphereSource* theROI)
{
  roiSphere= theROI; 
}
void ROIDrag::SetRenderWindow(vtkRenderWindow* renderWindow)
{
  renw= renderWindow; 
}
ROIDrag::ROIDrag()
{

}
ROIDrag::~ROIDrag()
{

}


