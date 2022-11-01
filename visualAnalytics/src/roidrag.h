/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/


#ifndef ROIDRAG_H
#define ROIDRAG_H


 
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkCommand.h"
#include "vtkSphereSource.h" 
#include "vtkRenderWindow.h"
#include <QObject>
#include <QMainWindow>

class ROIDrag : public QMainWindow
{
  Q_OBJECT
  
  public:
	ROIDrag(); 
	~ROIDrag(); 
        void SetInteractor(vtkRenderWindowInteractor* interactor); 
	void SetRenderer(vtkRenderer* theRen);
	void SetROISource(vtkSphereSource* theROI); 
	void SetRenderWindow(vtkRenderWindow* renderWindow); 
	
  public slots:
  void MyMouseMove(void); 
  
  
  private:

        vtkRenderWindowInteractor *iren;
	vtkRenderer* ren; 
	vtkRenderWindow* renw; 
	vtkSphereSource* roiSphere;
	int cx,cy; 
	double* worldCoordinate; 	
};

#endif // ROIDRAG_H
