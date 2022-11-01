/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/
#include "unit.h"
#include "ui_unit.h"
#include "kaimingCommon.h"
#include "sstream"
#include "vtkTextProperty.h"
#include <algorithm>
#include <QMenu>
#include <QtGui>
#include "QDataStream"
#include "vtkInteractorStyle.h"
#include "vtkTDxInteractorStyleCamera.h"
#include "vtkTDxInteractorStyleSettings.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkCommand.h"
#include "vtkConeSource.h"
#include "vtkCamera.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataReader.h"
#include "vtkActor.h"
#include "vtkMetaImageReader.h"
#include "vtkImageViewer2.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRegressionTestImage.h"
#include "vtkImageData.h"
#include "vtkEventQtSlotConnect.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkChartXY.h"
#include "vtkPlot.h"
#include "vtkTable.h"
#include "vtkFloatArray.h"
#include "vtkContextView.h"
#include "vtkContextScene.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRegressionTestImage.h"
#include <vtkActorCollection.h>
#include "vtkCylinderSource.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "newimageall.h"
#include "kaimingCommon.h"
#include <vtkSphereSource.h>
#include "Vector3D.h"
#include "triSurface.h"
#include "fibers.h"
#include "vtkCellArray.h"
#include "vtkProperty.h"
#include "vtkInteractorStyleRubberBandPick.h"
#include "newmat.h"
#include "vtkLookupTable.h"
#include "vtkPropPicker.h"
#include "vtkPointPicker.h"
#include "vtkCellPicker.h"
#include "vtkArrowSource.h"
#include "grangerCausality.h"
#include "vtkTextSource.h"
#include "vtkTextActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "kaimingCommon.h"
#include "vtkAxis.h"

using namespace KML;
using namespace NEWIMAGE;
using namespace NEWMAT;
using namespace Ui;

#ifndef VTK_CREATE
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()
#endif

#ifndef MyDestroy
#define MyDestroy(name) \
  if(NULL != name ) { delete name; }
#endif

#ifndef VTKDestroy
#define VTKDestroy(name) \
  if(NULL != name ) { name->Delete(); }
#endif

CUnit::CUnit()
{
    this->setupUi(this);

    //init values;
    surfRenderer=NULL;
    fiberRenderer=NULL;
    linkCameraConnection=NULL;
    mouseMoveConnection=NULL;
    releaseMouseMoveConnection=NULL;
    actionAdjust_ROIs=NULL;
//   actionFunctional_Network=NULL;
//   actionEffective_Network=NULL;
    volViewX=NULL;
    volViewY=NULL;
    volViewZ=NULL;
    boldView=NULL;
    linePlot=NULL;
    roiDefDial=NULL;
    myFiber=NULL;
    mySurface=NULL;
    myModel=NULL;
    boldSignalTable=NULL;
    boldSignalChart=NULL;
    pPicker=NULL;
    pTemplateWnd=NULL;
    textSource=NULL;
    textActor=NULL;
    textMapper=NULL;
    roiPicked=-1;


    b0Spacing[0]=0;
    b0Spacing[1]=0;
    b0Spacing[2]=0;

    surf_KEY= string("SURFACE");
    dtiB0_KEY= string("DTIB0");
    bolds_KEY= string("BOLDS");
    fiber_KEY= string("FIBERS");
    network_KEY= string("NETWORK");
    dtiseg_KEY=string("DTIGM");
    t1_KEY=string("T1");
    t1Label_KEY=string("T1LABEL");
    atlas_KEY=string("ATLAS");


    roiListArea = NULL;
    isSingleROIMode=true;
    isVolT1=false;

    linkCameraConnection = vtkEventQtSlotConnect::New();
    mouseMoveConnection= vtkEventQtSlotConnect::New();
    releaseMouseMoveConnection= vtkEventQtSlotConnect::New();
    extractor=vtkSmartPointer<vtkExtractPolyDataGeometry>::New();
    extractorSphere=vtkSmartPointer<vtkSphere>::New();
    extractor->SetImplicitFunction(extractorSphere);
    extractor->ExtractInsideOn();

    modeType=EffectiveConnectivity;

    actionAdjust_ROIs = new QAction(this);
    actionAdjust_ROIs->setObjectName(QString::fromUtf8("actionAdjust_ROIs"));
    actionAdjust_ROIs->setEnabled(false);

    surfWidget->GetRenderWindow()->AddRenderer(vtkSmartPointer<vtkRenderer>::New());
    BoldWidget->GetRenderWindow()->AddRenderer(vtkSmartPointer<vtkRenderer>::New());
    fiberWidget->GetRenderWindow()->AddRenderer(vtkSmartPointer<vtkRenderer>::New());
    volXWidget->GetRenderWindow()->AddRenderer(vtkSmartPointer<vtkRenderer>::New());
    volYWidget->GetRenderWindow()->AddRenderer(vtkSmartPointer<vtkRenderer>::New());
    volZWidget->GetRenderWindow()->AddRenderer(vtkSmartPointer<vtkRenderer>::New());

//   actionFunctional_Network = new QAction(this);
//   actionFunctional_Network->setObjectName(QString::fromUtf8("actionFunctional_Network"));
//   actionFunctional_Network->setCheckable(true);
//   actionFunctional_Network->setChecked(false);
//   actionFunctional_Network->setEnabled(true);
//   actionEffective_Network = new QAction(this);
//   actionEffective_Network->setObjectName(QString::fromUtf8("actionEffective_Network"));
//   actionEffective_Network->setCheckable(true);
//   actionEffective_Network->setChecked(true);
//
//   connect(actionFunctional_Network, SIGNAL(triggered(bool)),this, SLOT(Change2FunctionalConnectivity()));
//   connect(actionEffective_Network,SIGNAL(triggered(bool)),this,SLOT(Change2EffectiveConnectivity()));

//   setWindowTitle(tr("Dynamic Brain Network Visualization Toolkit"));

    roiDefDial= new myRoiDefDdial();
//   connect(actionDefine_Now, SIGNAL(triggered(bool)), this, SLOT(ShowROIDefDial()));
    connect(roiDefDial,SIGNAL(accepted()),this,SLOT(CreatandVisualizeSingleROI()));
//   connect(actionLoad_Profile, SIGNAL(triggered(bool)), this, SLOT(LoadProfiles()));
//   connect(actionSave_Network, SIGNAL(triggered(bool)), this, SLOT(SaveNetwork()));
//   connect(actionGenerate_Network, SIGNAL(triggered(bool)),this, SLOT(GenerateNetworkViaDefinedROIS()));


    QObject::connect(lineEditX, SIGNAL(textChanged(QString)), this, SLOT(UpdateSliceX(const QString& )));
    QObject::connect(lineEditY, SIGNAL(textChanged(QString)), this, SLOT(UpdateSliceY(const QString& )));
    QObject::connect(lineEditZ, SIGNAL(textChanged(QString)), this, SLOT(UpdateSliceZ(const QString& )));

    connect(lineEditX,SIGNAL(textChanged(QString)),this,SLOT(SetSliderValueX(QString)));
    connect(lineEditY,SIGNAL(textChanged(QString)),this,SLOT(SetSliderValueY(QString)));
    connect(lineEditZ,SIGNAL(textChanged(QString)),this,SLOT(SetSliderValueZ(QString)));

    connect(horizontalSliderX,SIGNAL(valueChanged(int)),this, SLOT(SetLineEditX(int)));
    connect(horizontalSliderY,SIGNAL(valueChanged(int)),this, SLOT(SetLineEditY(int)));
    connect(horizontalSliderZ,SIGNAL(valueChanged(int)),this, SLOT(SetLineEditZ(int)));
}
void CUnit::MyMouseMove(void )
{

    QVTKInteractor* iren= surfWidget->GetInteractor();
    vtkPropPicker* pker=  static_cast< vtkPropPicker*> ( iren->GetPicker());

    if (pker && pker->GetActor()&& iren->GetControlKey()) {
        int cx,cy; //mouse display position;
        pker->GetActor()->SetPickable(0);

        vtkCellPicker* cPicker= vtkCellPicker::New();
        iren->SetPicker(cPicker);
        cPicker->SetTolerance(0.02);
        iren->GetEventPosition(cx,cy);

        roiPicked= this->GetROIID(pker->GetActor());
        if (cPicker->Pick(cx,cy,0,surfRenderer)) {
            if (!isSingleROIMode) {
                double* worldCoordinate=cPicker->GetPickPosition();
                // 	    cout<< worldCoordinate[0]<<" "<<worldCoordinate[1]<<" "<<worldCoordinate[2]<<" "<<cPicker->GetPointId()<<endl;
                allSphereSources[roiPicked]->SetCenter(worldCoordinate);
                allNetworkNodes[roiPicked].x= worldCoordinate[0];
                allNetworkNodes[roiPicked].y= worldCoordinate[1];
                allNetworkNodes[roiPicked].z= worldCoordinate[2];

                allSphereSources[roiPicked]->Update();

                //	    cout<<roiPicked<<endl;

                surfWidget->update();
                this->GetPCA1stCurrentROI(roiPicked,(*allroiAdjs[roiPicked]).GetRadius(),allNodesPCA1st[roiPicked]);

                if (modeType==FunctionalConnectivity)
                    this->UpdateFunctionalNetworkWeights(roiPicked);
                else
                    this->UpdateEffectiveNetworkWeights(roiPicked);


                allroiAdjs[roiPicked]->checkBox_BOLD->setChecked(true);
                allroiAdjs[roiPicked]->checkBox_Fiber->setChecked(true);

                //	    cout<<"UpdateNetworkVisualization"<<endl;
                this->UpdateNetworkVisualization(*allroiAdjs[roiPicked]);

                allROIAttributeStatus[roiPicked] = NeedUpdate;

                this->ChangeROIFIberVisiability(*allroiAdjs[roiPicked]);

                this->UpdateBOLDSignal4Node(roiPicked, allroiAdjs[roiPicked]->GetRadius());

                //update the slices;
                this->UpdateVolViewSlices(roiPicked);

                //update status;
                this->UpdateAtlasInfoAtStatusBar(roiPicked);

                if ( NULL!= pTemplateWnd )
                {
		  cout<<"model: "<<pTemplateWnd->allNetworkNodes.size()<<endl;

                    if (textActor)
                    {
                        fiberRenderer->RemoveActor(textActor);
                        surfRenderer->RemoveActor(textActor);
                        if (NULL==textActor)
                            textActor->Delete();
                        fiberWidget->update();
                        surfWidget->update();
                    }

                    this->UpdateText(roiPicked,Vector3D<double>(worldCoordinate[0],worldCoordinate[1],worldCoordinate[2]));
                }
            } // for network movement;
            else {
                cout<<"Node picked: "<<roiPicked<<endl;
                double* worldCoordinate=cPicker->GetPickPosition();
                allSphereSources[roiPicked]->SetCenter(worldCoordinate);
                allNetworkNodes[roiPicked].x= worldCoordinate[0];
                allNetworkNodes[roiPicked].y= worldCoordinate[1];
                allNetworkNodes[roiPicked].z= worldCoordinate[2];
                allSphereSources[roiPicked]->Update();

                this->GetPCA1stCurrentROI(roiPicked,(*allroiAdjs[roiPicked]).GetRadius(),allNodesPCA1st[roiPicked]);

                this->UpdateBoldsSingleROI(allNetworkNodes[roiPicked],(*allroiAdjs[roiPicked]).GetRadius());
                allroiAdjs[roiPicked]->checkBox_BOLD->setChecked(true);
                allROIAttributeStatus[roiPicked] = NeedUpdate;
                allroiAdjs[roiPicked]->checkBox_Fiber->setChecked(true);
                this->VisualizeFibersCurrentROI(roiPicked);
//     	        this->UpdateBOLDNetROI(roiPicked, allroiAdjs[roiPicked]->GetRadius());
                //update the slices;
                this->UpdateVolViewSlices(roiPicked);

                surfWidget->update();
                fiberWidget->update();
                volXWidget->update();
                volYWidget->update();
                volZWidget->update();
            } // for single roii;
        }
        cPicker->Delete();
        iren->SetPicker(pker);
    }
}
int CUnit::CreatandVisualizeSingleROI(void )
{

//   isSingleROIMode=true;
    if (b0Spacing[0]==0)
    {
        QMessageBox::critical(this,tr("Be patient"),tr("Load dti b0 volume first!"), QMessageBox::Ok);
        return 1;
    }
    if (0==myBolds.xsize())
    {
        QMessageBox::critical(this,tr("Be patient"),tr("Load fmri bolds first!"), QMessageBox::Ok);
        return 1;
    }
    if (NULL==surfRenderer)
    {
        QMessageBox::critical(this,tr("Be patient"),tr("Load surface first!"), QMessageBox::Ok);
        return 1;
    }
    if (NULL==fiberRenderer)
    {
        QMessageBox::critical(this,tr("Be patient"),tr("Load fibers first!"), QMessageBox::Ok);
        return 1;
    }

    if (NULL==myModel)
    {
        myModel= new CJointModel( *mySurface,*myFiber, myBolds, myDtiGM);
        myModel->SetOffSet(b0Offset[0],b0Offset[1],b0Offset[2]);
        myModel->FindFiberSurfIntersections();
        myModel->FindFibersOnCells();
        myModel->WarpingBoldsOnSurface();
    }



    if (0==roiDefDial->tabWidget->currentIndex())
    {

        if (true==roiDefDial->radioButtonSurfID->isChecked())
        {
            QString surfID= roiDefDial->lineEditSurfID->text();
            int id= surfID.toInt();

            if (0> id || id >= mySurface->GetNumOfPoints() )
            {
                QMessageBox::critical(this, tr("Invalid Surface ID"),
                                      tr("The surface ID is invalid.\n"
                                         "Please enter a valid number.") ,QMessageBox::Cancel);
                return 1;

            }

            const Vector3D<float>& surfCoords= mySurface->GetPointCoords(id);
            float radius=roiDefDial->doubleSpinBoxRadius->value();

            if (radius<=0 )
            {
                QMessageBox::critical(this, tr("Invalid radius"),
                                      tr("The radius is invalid.\n"
                                         "Please enter a valid number.") ,QMessageBox::Cancel);
                return 1;

            }


            //remove the previous one:
            if (NULL!=preROIActor)
            {
                surfRenderer->RemoveActor(preROIActor);
                fiberRenderer->RemoveActor(preROIActor);
                volViewX->GetRenderer()->RemoveActor(preROIActor);
                volViewY->GetRenderer()->RemoveActor(preROIActor);
                volViewZ->GetRenderer()->RemoveActor(preROIActor);
                fiberWidget->update();
                surfWidget->update();
            }
            else
                preROIActor=vtkSmartPointer<vtkActor>::New();


            QString roiName= roiDefDial->lineEdit->text();

            VTK_CREATE(vtkSphereSource,sphere);
            sphere->SetCenter(surfCoords.x,surfCoords.y,surfCoords.z);
            sphere->SetRadius(radius);
            sphere->SetPhiResolution(100);
            sphere->SetThetaResolution(100);


            VTK_CREATE(vtkPolyDataMapper, sphereMapper);
            sphereMapper->SetInput(sphere->GetOutput());
            sphereMapper->SetColorModeToDefault();
            currROIActor= vtkSmartPointer<vtkActor>::New();
            currROIActor->SetMapper(sphereMapper);
            currROIActor->GetProperty()->SetColor(1.0000, 0, 0);

            allSphereSources.push_back(sphere);
            allNetworkNodesActors.push_back(currROIActor);
            allNetworkNodes.push_back(surfCoords);
            allROINames.push_back(roiName.toStdString());

            surfRenderer->AddActor(currROIActor);
            surfWidget->update();



            int index= allROINames.size()-1;
            myROIAdj* currentAdj= new myROIAdj;
            currentAdj->SetROIName(allROINames[index]);
            currentAdj->SetROIID(index);
            currentAdj->SetCoordinates(surfCoords.x, surfCoords.y, surfCoords.z);
            currentAdj->checkBox_BOLD->setChecked(false);
            currentAdj->checkBox_Fiber->setChecked(false);
            allroiAdjs.push_back(currentAdj);


            allROIAttributeStatus.push_back(NotDone);
            allROIFIbersVTK.push_back(vector<vtkSmartPointer<vtkIdList> >());
            allROIFiberIDs.push_back(vector<int>());
            this->TrackFibersCurrentROI(*currentAdj);


            vector<float> firstEignVector;
            this->GetPCA1stCurrentROI(index,currentAdj->GetRadius(),firstEignVector);
            allNodesPCA1st.push_back(firstEignVector);
            this->GenerateBOLDCurveofROI(index);

            this->VisualizeCurrentROIOnVolume(surfCoords,radius);
//      this->UpdateBoldsSingleROI(surfCoords,radius);
            this->VisualizeFibersCurrentROI(index);

            preROIActor=currROIActor;
        } // define using surface ID;


    }

    return 0;

}
void CUnit::VisualizeCurrentROIOnVolume(const Vector3D<float>& centerCoords, float radius )
{
    int lx=(int)( (centerCoords.x-b0Offset[0])/b0Spacing[0] +0.5);
    int ly=(int)( (centerCoords.y-b0Offset[1])/b0Spacing[1] +0.5);
    int lz=(int)( (centerCoords.z-b0Offset[2])/b0Spacing[2] +0.5);


    this->SetLineEditX(lx);
    this->SetLineEditY(ly);
    this->SetLineEditZ(lz);

    volViewX->SetSlice(lx);
    volViewY->SetSlice(ly);
    volViewZ->SetSlice(lz);

    if (NULL!= preROIActor)
    {
        volViewX->GetRenderer()->RemoveActor(preROIActor);
        volViewY->GetRenderer()->RemoveActor(preROIActor);
        volViewZ->GetRenderer()->RemoveActor(preROIActor);
    }


//   currROIActor->upda();

    volViewX->GetRenderer()->AddActor(currROIActor);
    volViewY->GetRenderer()->AddActor(currROIActor);
    volViewZ->GetRenderer()->AddActor(currROIActor);


    volXWidget->update();
    volYWidget->update();
    volZWidget->update();

}
void CUnit::VisualizeFibersCurrentROI(int roiID)
{

    if (allROIAttributeStatus[roiID]!= Updated)
    {
        cout<<"tracking roi: "<<roiID<<endl;
        TrackFibersCurrentROI(*allroiAdjs[roiID]);
    }

    VTK_CREATE(vtkCellArray,myArray);
    myArray->Allocate(5000);
    for (int findex= 0; findex < allROIFIbersVTK[roiID].size(); ++ findex)
    {
        myArray->InsertNextCell(allROIFIbersVTK[roiID][findex]);
    }


    fiberData->SetLines(myArray);

    if (NULL!=preROIActor)
        fiberRenderer->RemoveActor(preROIActor);
    fiberRenderer->AddActor(currROIActor);
    fiberWidget->update();

}
void CUnit::UpdateBoldsSingleROI(const KML::Vector3D< float >& centerCoords, float radius)
{
    vector<float> firstEigenVector;
    this->GetPCA1stCurrentROI(centerCoords,radius,firstEigenVector);
    //visualize the first principal component ;
    this->VisualizeBoldSingalROI(firstEigenVector);
}
void CUnit::UpdateBOLDSignal4Node(int roiID, float radius)
{
    vector<float> firstEigenVector;
    this->GetPCA1stCurrentROI(roiID,radius,firstEigenVector);
//   cout<<"UpdateBOLDNetROI"<<roiID<<endl;
    allNodesPCA1st[roiID] = firstEigenVector;
//   for(int index=0; index< firstEigenVector.size(); ++index)
//     cout<<firstEigenVector[index]<<" ";
//   cout<<endl;
//  cout<<firstEigenVector<<endl;
    this->GenerateBOLDCurveofROI(roiID);
    this->ChangeROIBOLDSVisiability(*allroiAdjs[roiID]);
}
void CUnit::GetPCA1stCurrentROI(const KML::Vector3D< float >& centerCoords, float radius, std::vector< float >& firstEigenVector)
{

    cout<<"GetPCA1stCurrentROI"<<endl;
    vector<Vector3D<int> >  allVoxels;


    int lx=(int)( (centerCoords.x-radius-b0Offset[0])/myBolds.xdim() +0.5);
    int ly=(int)( (centerCoords.y-radius-b0Offset[1])/myBolds.ydim() +0.5);
    int lz=(int)( (centerCoords.z-radius-b0Offset[2])/myBolds.zdim() +0.5);

    int ux=(int)( (centerCoords.x+radius-b0Offset[0])/myBolds.xdim() +0.5);
    int uy=(int)( (centerCoords.y+radius-b0Offset[1])/myBolds.ydim() +0.5);
    int uz=(int)( (centerCoords.z+radius-b0Offset[2])/myBolds.zdim() +0.5);

    //defensive
    lx=lx>0 ? lx : 0;
    ly=ly>0 ? ly : 0;
    lz=lz>0 ? lz : 0;

    ux=ux>0 ? ux : 0;
    uy=uy>0 ? uy : 0;
    uz=uz>0 ? uz : 0;


    lx=lx>myBolds.xsize() ? myBolds.xsize()-1 : lx;
    ly=ly>myBolds.ysize() ? myBolds.ysize()-1 : ly;
    lz=lz>myBolds.zsize() ? myBolds.zsize()-1 : lz;

    ux=ux>myBolds.xsize() ? myBolds.xsize()-1 : ux;
    uy=uy>myBolds.ysize() ? myBolds.ysize()-1 : uy;
    uz=uz>myBolds.zsize() ? myBolds.zsize()-1 : uz;


    for (int xid= lx; xid<= ux; ++xid)
        for (int yid = ly; yid<= uy; ++yid)
            for (int zid= lz; zid<= uz; ++zid)
            {
                if ( myBolds(xid,yid,zid,0) )
                {
                    Vector3D<int> roiVoxel(xid,yid,zid);
                    allVoxels.push_back(roiVoxel);
                }
            }

    //now we have all the roi voxels, calculate the pcas of these rois;

    vector<float> evalues(allVoxels.size());
    vector<vector<float> > evectors(allVoxels.size(), vector<float> ( myBolds.tsize(),0));
    this->DoPCA(allVoxels,evalues,evectors);

    firstEigenVector.clear();
    firstEigenVector.resize(evectors.size());
    for (int idx=0; idx< firstEigenVector.size(); idx++)
        firstEigenVector[idx]= evectors[idx][0];
    KML::StatsNormalizeFisher(firstEigenVector);

}
int CUnit::GetPCA1stCurrentROI ( int roiID, float radius, std::vector< float >& firstEigenVector )
{
    extractorSphere->SetCenter( allNetworkNodes[roiID].x,  allNetworkNodes[roiID].y,  allNetworkNodes[roiID].z);
    extractorSphere->SetRadius(radius);
    vtkPolyData* outData= extractor->GetOutput();
    outData->Update();

    set<int> allIDs;
    for (int cellID = 0; cellID < outData->GetNumberOfCells(); ++cellID) {
        allIDs.insert(outData->GetCell(cellID)->GetPointId(0));
        allIDs.insert(outData->GetCell(cellID)->GetPointId(1));
        allIDs.insert(outData->GetCell(cellID)->GetPointId(2));
    } //end of for loop:: pointID

    //cout<<"number of voxels for PCA: "<<allIDs.size()<<endl;
    if (1==allIDs.size())
    {
        firstEigenVector= myModel->GetVertexTimeSeris(*allIDs.begin());
        return 0;
    }
    if (0==allIDs.size())
    {
        this->GetPCA1stCurrentROI(allNetworkNodes[roiID], radius, firstEigenVector);
        return 0;
    }

    vector<vector<float> > theEVectors;
    vector<float> evalues;
    this->DoPCA(allIDs,evalues, theEVectors);
    firstEigenVector.clear();
    firstEigenVector.resize(theEVectors.size());
    for (int idx=0; idx< firstEigenVector.size(); idx++)
        firstEigenVector[idx]= theEVectors[idx][0];

    KML::StatsNormalizeFisher(firstEigenVector);
    return 0;
}
void CUnit::DoPCA(const std::vector< KML::Vector3D<int> >& roiVoxels, std::vector< float >& evalue, vector< std::vector< float > >& evector)
{


    vector<vector<float> > rawData;
    for (int sampleIndex=0; sampleIndex< roiVoxels.size(); ++sampleIndex)
    {

        vector<float> currentSignal;
        for (int t=0; t<myBolds.tsize(); ++t)
        {
            currentSignal.push_back(myBolds( roiVoxels[sampleIndex].x,roiVoxels[sampleIndex].y, roiVoxels[sampleIndex].z, t));
        }

        rawData.push_back(currentSignal);
    }
    pcaMdl.SetData(rawData);
    pcaMdl.DoPCA(evector);

}
int CUnit::DoPCA(const std::set<int>& surfID, std::vector< float >& evalue, vector< std::vector< float > >& evector)
{

    vector<vector<float> > rawData;
    for (set<int>::iterator it= surfID.begin(); it!= surfID.end(); ++it)
    {

        vector<float>& currentBOLD= myModel->GetVertexTimeSeris(*it);
        rawData.push_back(currentBOLD);
    }
    pcaMdl.SetData(rawData);
    pcaMdl.DoPCA(evector);

    return 0;

}
void CUnit::ShowROIDefDial(void) {
    roiDefDial->show();
}
void CUnit::SetSliderValueX(const QString& newText)
{
    horizontalSliderX->setValue(newText.toInt());
}
void CUnit::SetSliderValueY(const QString& newText)
{
    horizontalSliderY->setValue(newText.toInt());
}
void CUnit::SetSliderValueZ(const QString& newText)
{
    horizontalSliderZ->setValue(newText.toInt());
}
void CUnit::SetLineEditX(int val)
{
    QString stringValue;
    lineEditX->setText(stringValue.setNum(val));
}
void CUnit::SetLineEditY(int val)
{
    QString stringValue;
    lineEditY->setText(stringValue.setNum(val));

}
void CUnit::SetLineEditZ(int val)
{
    QString stringValue;
    lineEditZ->setText(stringValue.setNum(val));

}
CUnit::~CUnit()
{
    /* if(NULL!=surfRenderer)
       surfRenderer->Delete();
      if(NULL!=fiberRenderer)
       fiberRenderer->Delete();
     if(NULL!=linkCameraConnection)
       linkCameraConnection->Delete();
     if(NULL!=volViewX)
       volViewX->Delete();
     if(NULL!=volViewY)
       volViewY->Delete();
     if(NULL!=volViewZ)
       volViewZ->Delete();
     if(NULL!= boldView)
       boldView->Delete();
      if(NULL!= myFiber)
        delete myFiber;
      if(NULL!=mySurface)
        delete mySurface;
    //   if(NULL!= fiberData)
    //     fiberData->Delete();
    //   if(NULL!= linePlot)
    //     linePlot->Delete();

     if(NULL!= boldSignalChart)
       boldSignalChart->Delete();
     if(NULL!= boldSignalTable)
       boldSignalTable->Delete();*/

// MyDestroy(roiDefDial);
// MyDestroy(myModel);


}
void CUnit::VisualizeSurface(const std::string& myFileName )
{

    // create a window to make it stereo capable and give it to QVTKWidget
    if (NULL== surfRenderer)
    {

        vtkSmartPointer<vtkRenderWindow> surfWin = vtkSmartPointer<vtkRenderWindow>::New();
        surfWin->StereoCapableWindowOn();
        surfWidget->SetUseTDx(true);
        surfWidget->SetRenderWindow( surfWin);
        surfRenderer=vtkRenderer::New();
        surfWidget->GetRenderWindow()->AddRenderer(surfRenderer);


        //show popup
        QMenu* popupLinkCamera = new QMenu(surfWidget);
        popupLinkCamera->addAction("&LinkCamera");
        popupLinkCamera->addAction("&DisconnectCamera");
        popupLinkCamera->addAction("&AddFiberOverlay");
        popupLinkCamera->addAction("&RemoveFiberOverlay");
        connect(popupLinkCamera, SIGNAL(triggered(QAction*)), this, SLOT(PopupSurfaceWind(QAction*)));
        linkCameraConnection->Connect(surfWidget->GetRenderWindow()->GetInteractor(), vtkCommand::RightButtonPressEvent, this,
                                      SLOT(Popup( vtkObject*, unsigned long, void*, void*, vtkCommand*)), popupLinkCamera, 1.0);


        pPicker= vtkSmartPointer<vtkPropPicker>::New();
        surfWidget->GetInteractor()->SetPicker(pPicker);
        surfWidget->setMouseTracking(true);
        mouseMoveConnection->Connect(surfWidget->GetRenderWindow()->GetInteractor(), vtkCommand::MouseMoveEvent, this, SLOT(MyMouseMove()));
        releaseMouseMoveConnection->Connect(surfWidget->GetRenderWindow()->GetInteractor(), vtkCommand::KeyReleaseEvent, this, SLOT(ReleaseMouseMove()));


    }

    surfRenderer->RemoveAllViewProps();


    vtkSmartPointer<vtkPolyDataReader> surfReader=vtkSmartPointer<vtkPolyDataReader>::New();
    surfReader->SetFileName(myFileName.c_str());
    vtkSmartPointer<vtkPolyDataMapper> surfMapper=vtkSmartPointer<vtkPolyDataMapper>::New();
    surfMapper->SetInput(surfReader->GetOutput());
    surfMapper->SetScalarModeToUseCellData();

    surfData= surfReader->GetOutput();
    extractor->SetInputConnection(surfReader->GetOutputPort());

    surfActor=vtkSmartPointer<vtkActor>::New();
    surfActor->SetMapper(surfMapper);
    surfRenderer->AddActor(surfActor);


    if (NULL!=mySurface && NULL==myModel)
    {
        delete mySurface;
        mySurface=NULL;
    }
    mySurface= new CTriSurface(myFileName.c_str());
}

void CUnit::PopupSurfaceWind(QAction* obj)
{
    if ( obj->iconText()==QString("LinkCamera"))
    {
        cout<<"LinkCamera"<<endl;
        this->LinkCameraofSurfAndFiber(obj);
    } else if ( obj->iconText()==QString("DisconnectCamera"))
    {
        cout<<"DisconnectCamera"<<endl;
        this->DisconnectCamera();
    } else if (obj->iconText()==QString("AddFiberOverlay"))
    {
        cout<<"AddFiberOverlay"<<endl;
        surfRenderer->AddActor(fiberActor);
        surfWidget->update();

    } else if (obj->iconText()==QString("RemoveFiberOverlay"))
    {
        cout<<"RemoveFiberOverlay"<<endl;
        surfRenderer->RemoveActor(fiberActor);
        surfWidget->update();
    }

}

void CUnit::PopupFiberWind(QAction* obj)
{
    if ( obj->iconText()==QString("LinkCamera"))
    {
        cout<<"LinkCamera"<<endl;
        this->LinkCameraofSurfAndFiber(obj);
    } else if ( obj->iconText()==QString("DisconnectCamera"))
    {
        cout<<"DisconnectCamera"<<endl;
        this->DisconnectCamera();
    } else if (obj->iconText()==QString("AddSurfaceOverlay"))
    {
        cout<<"AddSurfaceOverlay"<<endl;
        fiberRenderer->AddActor(surfActor);
        fiberWidget->update();

    } else if (obj->iconText()==QString("RemoveSurfaceOverlay"))
    {
        cout<<"RemoveSurfaceOverlay"<<endl;
        fiberRenderer->RemoveActor(surfActor);
        fiberWidget->update();
    }
}
void CUnit::LinkCameraofSurfAndFiber(QAction* obj)
{


    vtkSmartPointer<vtkCamera> commonCamera= vtkSmartPointer<vtkCamera>::New();

    surfRenderer->SetActiveCamera(commonCamera);
    fiberRenderer->SetActiveCamera(commonCamera);

    surfRenderer->ResetCamera();
    fiberRenderer->ResetCamera();
    linkCameraConnection->Connect(commonCamera, vtkCommand::AnyEvent, fiberWidget, SLOT(update()));
    linkCameraConnection->Connect(commonCamera, vtkCommand::AnyEvent, surfWidget, SLOT(update()));
}
void CUnit::ShowFullScreen(QAction* obj)
{

}
void CUnit::VisualizeFibers(const std::string& myFileName)
{

    if (NULL==fiberRenderer)
    {
        vtkSmartPointer<vtkRenderWindow>  fiberWin = vtkSmartPointer<vtkRenderWindow>::New();
        fiberWin->StereoCapableWindowOn();
        fiberWidget->SetUseTDx(true);
        fiberWidget->SetRenderWindow( fiberWin);

        // add a renderer
        fiberRenderer = vtkRenderer::New();
        fiberWidget->GetRenderWindow()->AddRenderer(fiberRenderer);

        //show popup
        QMenu* popupLinkCameraFiber = new QMenu(fiberWidget);
        popupLinkCameraFiber->addAction("&LinkCamera");
        popupLinkCameraFiber->addAction("&DisconnectCamera");
        popupLinkCameraFiber->addAction("&AddSurfaceOverlay");
        popupLinkCameraFiber->addAction("&RemoveSurfaceOverlay");

        connect(popupLinkCameraFiber, SIGNAL(triggered(QAction*)), this, SLOT(PopupFiberWind(QAction*)));

        linkCameraConnection->Connect(fiberRenderer->GetRenderWindow()->GetInteractor(), vtkCommand::RightButtonPressEvent,   this,
                                      SLOT(Popup( vtkObject*, unsigned long, void*, void*, vtkCommand*)),
                                      popupLinkCameraFiber, 1.0);
    }

    fiberRenderer->RemoveAllViewProps();

    vtkSmartPointer<vtkPolyDataReader> fiberReader=vtkSmartPointer<vtkPolyDataReader>::New();
    fiberReader->SetFileName(myFileName.c_str());
    fiberReader->Update();
    fiberData=fiberReader->GetOutput();
    vtkSmartPointer<vtkPolyDataMapper>  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInput(fiberData);
    fiberActor = vtkSmartPointer<vtkActor>::New();
    fiberActor->SetMapper(mapper);
    fiberRenderer->AddViewProp(fiberActor);


    if (NULL!=myFiber && NULL==myModel)
    {
        delete myFiber;
        myFiber=NULL;
    }
    myFiber= new CFibers(myFileName.c_str());

}
void CUnit::VisualizeB0Volume(const std::string& myFileName )
{
    isVolT1=false;

//   QObject::connect(lineEditX, SIGNAL(textChanged(QString)), this, SLOT(UpdateSliceX(const QString& )));
//   QObject::connect(lineEditY, SIGNAL(textChanged(QString)), this, SLOT(UpdateSliceY(const QString& )));
//   QObject::connect(lineEditZ, SIGNAL(textChanged(QString)), this, SLOT(UpdateSliceZ(const QString& )));
//
//   connect(lineEditX,SIGNAL(textChanged(QString)),this,SLOT(SetSliderValueX(QString)));
//   connect(lineEditY,SIGNAL(textChanged(QString)),this,SLOT(SetSliderValueY(QString)));
//   connect(lineEditZ,SIGNAL(textChanged(QString)),this,SLOT(SetSliderValueZ(QString)));
//
//   connect(horizontalSliderX,SIGNAL(valueChanged(int)),this, SLOT(SetLineEditX(int)));
//   connect(horizontalSliderY,SIGNAL(valueChanged(int)),this, SLOT(SetLineEditY(int)));
//   connect(horizontalSliderZ,SIGNAL(valueChanged(int)),this, SLOT(SetLineEditZ(int)));



    vtkSmartPointer< vtkMetaImageReader > volReader=vtkSmartPointer<vtkMetaImageReader>::New();
    volReader->SetFileName(myFileName.c_str());
    volReader->Update();

    vtkSmartPointer< vtkImageData > img=vtkSmartPointer<vtkImageData>::New();
    img= volReader->GetOutput();

    img->GetDimensions(b0Dimensions);
    img->GetSpacing(b0Spacing);
    img->GetOrigin(b0Offset);


//TODO: set intensity mapper;

//x view;
    //volXWidget->SetUseTDx(true);
    volViewX=vtkImageViewer2::New();
    volViewX->SetInput(volReader->GetOutput());
    volViewX->SetupInteractor(volXWidget->GetInteractor());
//   volViewX->GetRenderWindow()->Finalize();
//   volViewX->GetRenderWindow()->StereoCapableWindowOn();
    volXWidget->SetRenderWindow(volViewX->GetRenderWindow());

    volViewX->SetSliceOrientationToYZ();
    volViewX->SetSlice(b0Dimensions[0]/2);

//   volViewX->GetRenderer()->ResetCamera();
//   volViewX->GetRenderWindow()->StereoCapableWindowOn();
    volViewX->SetColorWindow(6000);
    volViewX->SetColorLevel(3000);



    // y view;
    //volYWidget->SetUseTDx(true);
    volViewY=vtkImageViewer2::New();
    volViewY->SetInput(volReader->GetOutput());
    volViewY->SetupInteractor(volYWidget->GetInteractor());
//   volViewY->GetRenderWindow()->Finalize();
//   volViewY->GetRenderWindow()->StereoCapableWindowOn();
    volYWidget->SetRenderWindow(volViewY->GetRenderWindow());

    volViewY->SetSliceOrientationToXZ();
    volViewY->SetSlice(b0Dimensions[1]/2);
    volViewY->SetColorWindow(6000);
    volViewY->SetColorLevel(3000);
//   volViewY->GetRenderer()->ResetCamera();
//   volViewY->GetRenderWindow()->StereoCapableWindowOn();

    //z view;
    //volZWidget->SetUseTDx(true);
    volViewZ=vtkImageViewer2::New();
    volViewZ->SetInput(volReader->GetOutput());
    volViewZ->SetupInteractor(volZWidget->GetInteractor());
//   volViewZ->GetRenderWindow()->Finalize();
//   volViewZ->GetRenderWindow()->StereoCapableWindowOn();
    volZWidget->SetRenderWindow(volViewZ->GetRenderWindow());

    volViewZ->SetSliceOrientationToXY();
    volViewZ->SetSlice(b0Dimensions[2]/2);
    volViewZ->SetColorWindow(6000);
    volViewZ->SetColorLevel(3000);


    horizontalSliderX->setMaximum(b0Dimensions[0]-1);
    horizontalSliderX->setMinimum(0);
    horizontalSliderY->setMaximum(b0Dimensions[1]-1);
    horizontalSliderY->setMinimum(0);
    horizontalSliderZ->setMaximum(b0Dimensions[2]-1);
    horizontalSliderZ->setMinimum(0);


//   QString sliceID;
//
//   sliceID.setNum(0);
//   lineEditX->setText(sliceID);
//   sliceID.setNum(0);
//   lineEditY->setText(sliceID);
//   sliceID.setNum(0);
//   lineEditZ->setText(sliceID);
//
//   sliceID.setNum(b0Dimensions[0]/2);
//   lineEditX->setText(sliceID);
//   sliceID.setNum(b0Dimensions[1]/2);
//   lineEditY->setText(sliceID);
//   sliceID.setNum(b0Dimensions[2]/2);
//   lineEditZ->setText(sliceID);

    this->SetLineEditX(0);
    this->SetLineEditX(b0Dimensions[0]/2);
    this->SetLineEditY(0);
    this->SetLineEditY(b0Dimensions[1]/2);
    this->SetLineEditZ(0);
    this->SetLineEditZ(b0Dimensions[2]/2);

//   QObject::connect(lineEditX, SIGNAL(textChanged(QString)), this, SLOT(UpdateSliceX(const QString& )));
//   QObject::connect(lineEditY, SIGNAL(textChanged(QString)), this, SLOT(UpdateSliceY(const QString& )));
//   QObject::connect(lineEditZ, SIGNAL(textChanged(QString)), this, SLOT(UpdateSliceZ(const QString& )));
//
//   connect(lineEditX,SIGNAL(textEdited(QString)),this,SLOT(SetSliderValueX(QString)));
//   connect(lineEditY,SIGNAL(textEdited(QString)),this,SLOT(SetSliderValueY(QString)));
//   connect(lineEditZ,SIGNAL(textEdited(QString)),this,SLOT(SetSliderValueZ(QString)));
//
//   connect(horizontalSliderX,SIGNAL(valueChanged(int)),this, SLOT(SetLineEditX(int)));
//   connect(horizontalSliderY,SIGNAL(valueChanged(int)),this, SLOT(SetLineEditY(int)));
//   connect(horizontalSliderZ,SIGNAL(valueChanged(int)),this, SLOT(SetLineEditZ(int)));

}

void CUnit::VisualizeT1Volume(const std::string& myFileName)
{

    isVolT1=true;

    vtkSmartPointer< vtkMetaImageReader > volReader=vtkSmartPointer<vtkMetaImageReader>::New();
    volReader->SetFileName(myFileName.c_str());
    volReader->Update();

    vtkSmartPointer< vtkImageData > img=vtkSmartPointer<vtkImageData>::New();
    img= volReader->GetOutput();

    img->GetDimensions(t1Dimensions);
    img->GetSpacing(t1Spacing);
    img->GetOrigin(t1Offset);

    volViewX=vtkImageViewer2::New();
    volViewX->SetInput(volReader->GetOutput());
    volViewX->SetupInteractor(volXWidget->GetInteractor());
    volXWidget->SetRenderWindow(volViewX->GetRenderWindow());

    volViewX->SetSliceOrientationToYZ();
    volViewX->SetSlice(t1Dimensions[0]/2);
    volViewX->SetColorWindow(6000);
    volViewX->SetColorLevel(3000);

    volViewY=vtkImageViewer2::New();
    volViewY->SetInput(volReader->GetOutput());
    volViewY->SetupInteractor(volYWidget->GetInteractor());
    volYWidget->SetRenderWindow(volViewY->GetRenderWindow());

    volViewY->SetSliceOrientationToXZ();
    volViewY->SetSlice(t1Dimensions[1]/2);
    volViewY->SetColorWindow(6000);
    volViewY->SetColorLevel(3000);

    volViewZ=vtkImageViewer2::New();
    volViewZ->SetInput(volReader->GetOutput());
    volViewZ->SetupInteractor(volZWidget->GetInteractor());
    volZWidget->SetRenderWindow(volViewZ->GetRenderWindow());

    volViewZ->SetSliceOrientationToXY();
    volViewZ->SetSlice(t1Dimensions[2]/2);
    volViewZ->SetColorWindow(6000);
    volViewZ->SetColorLevel(3000);

    horizontalSliderX->setMaximum(t1Dimensions[0]-1);
    horizontalSliderX->setMinimum(0);
    horizontalSliderY->setMaximum(t1Dimensions[1]-1);
    horizontalSliderY->setMinimum(0);
    horizontalSliderZ->setMaximum(t1Dimensions[2]-1);
    horizontalSliderZ->setMinimum(0);

    this->SetLineEditX(0);
    this->SetLineEditX(t1Dimensions[0]/2);
    this->SetLineEditY(0);
    this->SetLineEditY(t1Dimensions[1]/2);
    this->SetLineEditZ(0);
    this->SetLineEditZ(t1Dimensions[2]/2);

}


void CUnit::VisualizeBOLDS(const std::string& myFileName)
{


    int rst=read_volume4D( myBolds, myFileName.c_str());
    int mx=myBolds.xsize()/2;
    int my=myBolds.ysize()/2;
    int mz=myBolds.zsize()/2;


// this->VisualizeBoldSignal(mx,my,mz);

    //BoldWidget->SetUseTDx(true);
    // Set up a 2D scene, add an XY chart to it

    if (NULL==boldView)
    {
        boldView= vtkContextView::New();
        boldView->SetInteractor(  BoldWidget->GetInteractor());
        BoldWidget->SetRenderWindow(boldView->GetRenderWindow());
        BoldWidget->setMouseTracking(false);
        boldSignalChart=  vtkChartXY::New();
        boldView->GetScene()->AddItem(boldSignalChart);
        boldView->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
//       boldView->DisplayHoverTextOff();

        // Create a table with some points in it...
        boldSignalTable= vtkTable::New();
        VTK_CREATE(vtkFloatArray, arrX);
        arrX->SetName("Time Points");
        boldSignalTable->AddColumn(arrX);
        VTK_CREATE(vtkFloatArray, arrC);
        arrC->SetName("BOLDs");
        boldSignalTable->AddColumn(arrC);
        // Test charting with a few more points...
        int numPoints = myBolds.tsize();
        boldSignalTable->SetNumberOfRows(numPoints);
        for (int i = 0; i < numPoints; ++i)
        {
            boldSignalTable->SetValue(i, 0, i);
            boldSignalTable->SetValue(i, 1, /*bolds(mx,my,mz,i)*/ i);
        }
    }


    vector<float> centerBoldsingal( myBolds.tsize());
    for (int t=0; t<myBolds.tsize(); ++t)
    {
        centerBoldsingal[t]=myBolds(mx,my,mz,t);
    }

    KML::StatsNormalizeFisher(centerBoldsingal);
    this->VisualizeBoldSingalROI(centerBoldsingal);


}
void CUnit::VisualizeBoldSingalROI( const std::vector< float >& data)
{

    for (int i = 0; i < myBolds.tsize(); ++i)
    {
        boldSignalTable->SetValue(i, 0, i);
        boldSignalTable->SetValue(i, 1, data[i]);
    }

    if (linePlot!=NULL)
        boldSignalChart->RemovePlot(0);

    linePlot = boldSignalChart->AddPlot(vtkChart::LINE);
    linePlot->SetInput(boldSignalTable, 0, 1);
    linePlot->SetColor(0, 255, 0, 255);
    linePlot->SetWidth(2.0);
    linePlot->GetXAxis()->SetMaximum(myBolds.tsize()+2);
    linePlot->GetXAxis()->SetMinimum(0);

    linePlot->GetYAxis()->SetMaximum(4);
    linePlot->GetYAxis()->SetMinimum(-4);
    linePlot->GetYAxis()->AutoScale();

    boldSignalTable->Update();
    boldSignalChart->Update();
    boldView->Update();

    BoldWidget->update();

}

void CUnit::GenerateBOLDCurveofROI(int roiID)
{
//  cout<<"Number of GenerateBOLDCurveofROI "<< allROINames.size()<<endl;
    if (0!=roiID)
    {
        VTK_CREATE(vtkFloatArray, theData);
        theData->SetName(allROINames[roiID].c_str());
        theData->SetNumberOfValues( myBolds.tsize());
        boldSignalTable->AddColumn(theData);

    }

//  if(0==roiID)
//  {
////     boldSignalChart->RemovePlot(0);
//    boldSignalTable->GetColumn(0)->SetName(allROINames[roiID].c_str());
//  }


// cout<<"data:"<<endl<<allNodesPCA1st[roiID]<<endl;
    for (int i = 0; i < myBolds.tsize(); ++i)
    {
        boldSignalTable->SetValue(i, roiID+1, allNodesPCA1st[roiID][i]);
    }

// /*
//   linePlot = boldSignalChart->AddPlot(vtkChart::LINE);
//   linePlot->SetInput(boldSignalTable, 0, roiID+1);
//   linePlot->SetColor(myColorScheme.GetColorByIndex(roiID).x*255,myColorScheme.GetColorByIndex(roiID).y*255,myColorScheme.GetColorByIndex(roiID).z*255,255);
//   linePlot->SetWidth(2.0);*/
// //
//   boldSignalTable->Update();
//   boldSignalChart->Update();
//   BoldWidget->update();

}
void CUnit::Popup(vtkObject * obj, unsigned long,  void * client_data, void *, vtkCommand * command)
{
    // A note about context menus in Qt and the QVTKWidget
    // You may find it easy to just do context menus on right button up,
    // due to the event proxy mechanism in place.

    // That usually works, except in some cases.
    // One case is where you capture context menu events that
    // child windows don't process.  You could end up with a second
    // context menu after the first one.

    // See QVTKWidget::ContextMenuEvent enum which was added after the
    // writing of this example.

    // get interactor
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkRenderWindowInteractor::SafeDownCast(obj);
    // consume event so the interactor style doesn't get it
    command->AbortFlagOn();
    // get popup menu
    QMenu* popupMenu = static_cast<QMenu*>(client_data);
    // get event location
    int* sz = iren->GetSize();
    int* position = iren->GetEventPosition();
    // remember to flip y
    QPoint pt = QPoint(position[0], sz[1]-position[1]);
    // map to global
    QPoint global_pt = popupMenu->parentWidget()->mapToGlobal(pt);
    // show popup menu at global point
    popupMenu->popup(global_pt);


}
void CUnit::UpdateSliceX(const QString& newSliceNum)
{

    int sliceId= newSliceNum.toInt();
    volViewX->SetSlice(sliceId);
}
void CUnit::UpdateSliceY(const QString& newSliceNum)
{

    volViewY->SetSlice(newSliceNum.toInt());
    volViewY->UpdateDisplayExtent();
}
void CUnit::UpdateSliceZ(const QString& newSliceNum)
{
    volViewZ->SetSlice(newSliceNum.toInt());
    volViewZ->UpdateDisplayExtent();
}
void CUnit::DisconnectCamera(void )
{
    //check whether the camera is the same;
    if ( fiberRenderer->GetActiveCamera() == surfRenderer->GetActiveCamera())
    {
        vtkSmartPointer<vtkCamera> surfCamera= vtkSmartPointer<vtkCamera>::New();
        vtkSmartPointer< vtkCamera > fiberCamera= vtkSmartPointer<vtkCamera>::New();
        fiberRenderer->SetActiveCamera(fiberCamera);
        fiberRenderer->ResetCamera();
        surfRenderer->SetActiveCamera(surfCamera);
        surfRenderer->ResetCamera();
        fiberWidget->update();
        surfWidget->update();
    }
}
int CUnit::VisualizeNetwork(void )
{
    //some initialization;
    if (b0Spacing[0]==0)
    {
        QMessageBox::critical(this,tr("Be patient"),tr("Load dti b0 volume first!"), QMessageBox::Ok);
        return 1;
    }
    if (NULL==surfRenderer)
    {
        QMessageBox::critical(this,tr("Be patient"),tr("Load surface first!"), QMessageBox::Ok);
        return 1;
    }

    if (NULL==fiberRenderer)
    {
        QMessageBox::critical(this,tr("Be patient"),tr("Load fibers first!"), QMessageBox::Ok);
        return 1;
    }


    myColorScheme.SetColorNums(allNetworkNodes.size());;
    if (currROIActor) {

        surfRenderer->RemoveActor(currROIActor);
        surfActor->GetProperty()->SetOpacity(1);
        //surfActor->GetProperty()->SetColor(0.5,0.5,0.5);
        volViewX->GetRenderer()->RemoveActor(currROIActor);
        volViewY->GetRenderer()->RemoveActor(currROIActor);
        volViewZ->GetRenderer()->RemoveActor(currROIActor);

        surfWidget->update();
        fiberWidget->update();
        volXWidget->update();
        volYWidget->update();
        volZWidget->update();

    }

    fiberRenderer->RemoveAllViewProps();

    this->DoNetRendering(allNetworkNodes, edges);

//   actionAdjust_ROIs->setEnabled(true);
//   emit this->NetworkAdded();


    return 0;
}
void CUnit::DoNetRendering(const std::vector< KML::Vector3D< float > >& allROIs, const NEWMAT::Matrix& edges)
{

    this->AddNetworkNodesActors();
    allNodesLinksActors.resize(allROIs.size());
    this->AddNetworkLinksActors();

//   for(int line=1; line<= edges.Nrows(); ++line)
//     for(int col=1; col<= edges.Ncols(); ++col)
//     {
//       if(modeType==EffectiveConnectivity)
// 	this->AddArrow(line-1, col-1, allROIs, edges(line,col));
//       else
// 	this->AddLine(line-1, col-1, allROIs, edges(line,col));
//     }
//
    //now add bolds;
    for (int index=0; index< allNodesPCA1st.size(); ++index )
    {
        this->GenerateBOLDCurveofROI(index);
    }


    BoldWidget->update();
    surfWidget->update();
    fiberWidget->update();
}
int CUnit::AddArrow1(int x, int y, const vector<Vector3D<float> >& roiCenters, float opacity, float spheraRadius, float coneHeight, float cylinderRadius) {


    if (0.05 > opacity)
        return 0;

    if (opacity> 0.15)
        opacity=1;
    else
        opacity*=3 ;

    //this arrow points from x to y;
    Vector3D<float> direction=roiCenters[y]-roiCenters[x];

    float cylinderLength=direction.Norm();
    direction.Normalize();

    Vector3D<float> coneCenter=roiCenters[y]-(spheraRadius+0.5*coneHeight)*direction;
    Vector3D<float> cylinderCenter=0.5*(roiCenters[y]+roiCenters[x])/*-0.5*coneHeight*direction*/;

    vtkSmartPointer< vtkCylinderSource > cylinderSrc=vtkSmartPointer<vtkCylinderSource>::New();
    cylinderSrc->SetHeight(cylinderLength-2*spheraRadius-coneHeight);
    cylinderSrc->SetRadius(cylinderRadius);
    //get the angle and transform the cylinder;



    Vector3D<float> xAxis(1,0,0);
    Vector3D<float> yAxis(0,1,0);
    Vector3D<float> zAxis(0,0,1);

    double a= abs(DotProduct(xAxis,direction));
    double b=abs(DotProduct(zAxis,direction));
    double c=abs(DotProduct(yAxis,direction));

    float cosAngle= CosAngleBetweenVectors(yAxis,direction);
    float angleDegree= acos(cosAngle)/3.14*180;
//	angleDegree=angleDegree>90?  (180-angleDegree) : angleDegree;

    vtkSmartPointer< vtkTransform > transform=vtkSmartPointer<vtkTransform>::New();
    transform->Translate(cylinderCenter.x,cylinderCenter.y,cylinderCenter.z);
    Vector3D<float> axis= CrossProduct(direction,yAxis);

    transform->RotateWXYZ(-angleDegree,axis.x,axis.y,axis.z);

    vtkSmartPointer< vtkTransformPolyDataFilter > polyTransform=vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    polyTransform->SetInput(cylinderSrc->GetOutput());
    polyTransform->SetTransform(transform);
    polyTransform->Update();

    vtkSmartPointer< vtkPolyDataMapper > cylinderMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer< vtkActor > cylinderActor= vtkSmartPointer<vtkActor>::New();
    cylinderMapper->SetInput(polyTransform->GetOutput());
    cylinderActor->SetMapper(cylinderMapper);
    cylinderActor->GetProperty()->SetOpacity(abs(opacity));
    fiberRenderer->AddActor(cylinderActor);




    vtkSmartPointer< vtkConeSource > coneSrc=vtkSmartPointer<vtkConeSource>::New();
    coneSrc->SetDirection(direction.x,direction.y,direction.z);
    coneSrc->SetCenter(coneCenter.x,coneCenter.y,coneCenter.z);
    coneSrc->SetRadius(0.5*coneHeight);
    coneSrc->Update();


    vtkSmartPointer< vtkPolyDataMapper > coneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer< vtkActor > coneActor= vtkSmartPointer<vtkActor>::New();
    coneMapper->SetInput(coneSrc->GetOutput());
    coneActor->SetMapper(coneMapper);
    coneActor->GetProperty()->SetOpacity(opacity);

    allNodesLinksActors[x].push_back(coneActor);
    allNodesLinksActors[x].push_back(cylinderActor);
    allNodesLinksActors[y].push_back(coneActor);
    allNodesLinksActors[y].push_back(cylinderActor);

    fiberRenderer->AddActor(coneActor);
    fiberWidget->update();

    return 0;
}
int CUnit::AddArrow(int x, int y, const vector<Vector3D<float> >& roiCenters, float opacity, float spheraRadius, float coneHeight, float cylinderRadius) {


    if (0.05 > opacity)
        return 0;

    if (opacity> 0.15)
        opacity=1;
    else
        opacity*=3 ;

    float yzScale= opacity/0.05;

    //this arrow points from x to y;
    Vector3D<float> direction=roiCenters[y]-roiCenters[x];

    float arrowLength=direction.Norm()-spheraRadius;
    direction.Normalize();

    Vector3D<float> arrowCenter= roiCenters[x];


    vtkSmartPointer< vtkArrowSource > arrowSrc=vtkSmartPointer<vtkArrowSource>::New();
    arrowSrc->SetTipLength(coneHeight/arrowLength);
    arrowSrc->SetTipResolution(10);
    arrowSrc->SetTipRadius(2/yzScale);

    //get the angle and transform the cylinder;



    Vector3D<float> yAxis(1,0,0);
    float cosAngle= CosAngleBetweenVectors(yAxis,direction);
    float angleDegree= acos(cosAngle)/3.14*180;
//	angleDegree=angleDegree>90?  (180-angleDegree) : angleDegree;

    vtkSmartPointer< vtkTransform > transform=vtkSmartPointer<vtkTransform>::New();
    transform->Translate(arrowCenter.x,arrowCenter.y,arrowCenter.z);
    Vector3D<float> axis= CrossProduct(direction,yAxis);

    transform->RotateWXYZ(-angleDegree,axis.x,axis.y,axis.z);
    transform->Scale(arrowLength,yzScale,yzScale);

    vtkSmartPointer< vtkTransformPolyDataFilter > polyTransform=vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    polyTransform->SetInput(arrowSrc->GetOutput());
    polyTransform->SetTransform(transform);
    polyTransform->Update();

    vtkSmartPointer< vtkPolyDataMapper > cylinderMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer< vtkActor > cylinderActor= vtkSmartPointer<vtkActor>::New();
    cylinderMapper->SetInput(polyTransform->GetOutput());
    cylinderActor->SetMapper(cylinderMapper);
    cylinderActor->GetProperty()->SetOpacity(abs(opacity));
    fiberRenderer->AddActor(cylinderActor);

    allNodesLinksActors[x].push_back(cylinderActor);
    allNodesLinksActors[y].push_back(cylinderActor);

    fiberWidget->update();

    return 0;
}
int CUnit::AddLine(int x, int y, const std::vector< KML::Vector3D< float > >& roiCenters, float opacity, float spheraRadius, float cylinderRadius)
{

    if ( abs(opacity)< 0.3)
        return 0;
    else
        opacity= opacity/abs(opacity) * (abs(opacity)-0.3)/0.7;

    float scaleXZ= abs(opacity)/0.1;
    //this arrow points from x to y;
    Vector3D<float> direction=roiCenters[y]-roiCenters[x];

    float cylinderLength=direction.Norm();
    direction.Normalize();

    Vector3D<float> cylinderCenter=0.5*(roiCenters[y]+roiCenters[x]);

    vtkSmartPointer< vtkCylinderSource > cylinderSrc=vtkSmartPointer<vtkCylinderSource>::New();
    cylinderSrc->SetHeight(cylinderLength-2*spheraRadius);
    cylinderSrc->SetRadius(abs(opacity)/3);
    //get the angle and transform the cylinder;



    Vector3D<float> xAxis(1,0,0);
    Vector3D<float> yAxis(0,1,0);
    Vector3D<float> zAxis(0,0,1);

    double a= abs(DotProduct(xAxis,direction));
    double b=abs(DotProduct(zAxis,direction));
    double c=abs(DotProduct(yAxis,direction));

    float cosAngle= CosAngleBetweenVectors(yAxis,direction);
    float angleDegree= acos(cosAngle)/3.14*180;
//	angleDegree=angleDegree>90?  (180-angleDegree) : angleDegree;

    vtkSmartPointer< vtkTransform > transform=vtkSmartPointer<vtkTransform>::New();
    transform->Translate(cylinderCenter.x,cylinderCenter.y,cylinderCenter.z);
    Vector3D<float> axis= CrossProduct(direction,yAxis);

    transform->RotateWXYZ(-angleDegree,axis.x,axis.y,axis.z);
    transform->Scale(scaleXZ,1,scaleXZ);

    vtkSmartPointer< vtkTransformPolyDataFilter > polyTransform=vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    polyTransform->SetInput(cylinderSrc->GetOutput());
    polyTransform->SetTransform(transform);
    polyTransform->Update();

    vtkSmartPointer< vtkPolyDataMapper > cylinderMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer< vtkActor > cylinderActor= vtkSmartPointer<vtkActor>::New();
    cylinderMapper->SetInput(polyTransform->GetOutput());
    cylinderActor->SetMapper(cylinderMapper);
    cylinderActor->GetProperty()->SetOpacity(abs(opacity));
    if (opacity<0)
        cylinderActor->GetProperty()->SetColor(0,255,0);
    fiberRenderer->AddActor(cylinderActor);

    allNodesLinksActors[x].push_back(cylinderActor);
    allNodesLinksActors[y].push_back(cylinderActor);

    fiberWidget->update();

    return 0;


}
int CUnit::GetRoiCenterDTIVolCoords(const std::string& fileNanme, KML::Vector3D< float >& theCoords)
{

    volume<float> theVolume;
    read_volume(theVolume, fileNanme);
    Vector3D<float> gridCenter;
    int voxCount=0;

    for (int z=0; z< theVolume.zsize(); ++z)
        for (int y=0; y< theVolume.ysize(); ++y)
            for (int x=0; x<theVolume.xsize(); ++x)
            {
                if (theVolume(x,y,z))
                {
                    gridCenter.x+=x;
                    gridCenter.y+=y;
                    gridCenter.z+=z;
                    ++voxCount;
                }
            }

    if (0==voxCount)
    {
        cout<<"No ROI found in volume : "<<fileNanme<<endl;
        exit(1);
    }

    gridCenter/=voxCount;

    int cx= gridCenter.x+0.5;
    int cy= gridCenter.y+0.5;
    int cz= gridCenter.z+0.5;

    theCoords.x= cx* b0Spacing[0]+ b0Offset[0];
    theCoords.y= cy* b0Spacing[1]+ b0Offset[1];
    theCoords.z= cz* b0Spacing[2]+ b0Offset[2];

    return 0;

}

int CUnit::GetRoiCenterSurfCoords(const std::string& fileNanme, Vector3D< float >& theCoords) {

    this->GetRoiCenterDTIVolCoords(fileNanme, theCoords);
    int theId=this->GetRoiCenterSurfID(theCoords);
    theCoords= mySurface->GetPointCoords(theId);

    return 0;

}

int CUnit::GetRoiCenterSurfID(Vector3D< float >& theCoords) {

    float dis=10000000;
    int theId=-1;

    //find the closing surfid;
    for (int sid = 0; sid < mySurface->GetNumOfPoints(); ++sid) {
        Vector3D<float> tmp= theCoords-mySurface->GetPointCoords(sid);
        if (dis> tmp.Norm()) {
            dis= tmp.Norm();
            theId=sid;
        }
    } //end of for loop:: sid
    cout<<"Min distance: "<< dis<<endl;
    return theId;
}

void CUnit::AdjustROIS(void )
{
    this->ShowROIAdjs();
}
void CUnit::LoadBOLDS(void )
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open BOLDS of current subject"),"",tr("volume4D(*.nii.gz);;All Files(*)"));
    std::string myFileName= fileName.toStdString();
    bool hasFile= this->VaryfyInputFile(myFileName);
    if (hasFile)
        this->VisualizeBOLDS(myFileName);
    else
    {
        QMessageBox::critical(this, tr("Input NOT FOUND"), tr("Please input a fmri bolds file"),QMessageBox::Cancel);
    }
}
void CUnit::LoadGM(void )
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open gray matter mask of current subject"),"",tr("volume(*.nii.gz);;All Files(*)"));
    std::string myFileName= fileName.toStdString();
    bool hasFile= this->VaryfyInputFile(myFileName);
    if (hasFile)
        read_volume(myDtiGM,myFileName);
    else
    {
        QMessageBox::critical(this, tr("Input NOT FOUND"), tr("Please input a gm volume file"),QMessageBox::Cancel);
    }
}

void CUnit::LoadFiber(void )
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open surface of current subject"),"",tr("fiber(*.vtk);;All Files(*)"));
    std::string myFileName= fileName.toStdString();
    bool hasFile= this->VaryfyInputFile(myFileName);
    if (hasFile)
    {
        this->VisualizeFibers(myFileName);

    }
    else
    {
        QMessageBox::critical(this, tr("Input NOT FOUND"), tr("Please input a valid fiber file"),QMessageBox::Cancel);
    }
}
void CUnit::LoadSurface(void )
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open surface of current subject"),"",tr("surface (*.vtk);;All Files(*)"));
    std::string myFileName= fileName.toStdString();
    bool hasFile= this->VaryfyInputFile(myFileName);
    if (hasFile) {
        this->VisualizeSurface(myFileName);

    }
    else
    {
        QMessageBox::critical(this, tr("Input NOT FOUND"), tr("Please input a valid surface file"),QMessageBox::Cancel);
    }
}
void CUnit::LoadVolume(void )
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open volume of current subject"),"",tr("volume(*.mhd);;All Files(*)"));
    std::string myFileName= fileName.toStdString();
    bool hasFile= this->VaryfyInputFile(myFileName);
    if (hasFile)
        this->VisualizeB0Volume(myFileName);
    else
    {
        QMessageBox::critical(this, tr("Input NOT FOUND"), tr("Please input a valid volume file"),QMessageBox::Cancel);
    }
}
bool CUnit::VaryfyInputFile(const QString& fileName)
{
    std::string name= fileName.toStdString();
    return this->VaryfyInputFile(name);

}
bool CUnit::VaryfyInputFile(const std::string& fileName)
{
    if (0==fileName.size())
        return false;

    fstream testStrm;
    testStrm.open(fileName.c_str(), ios::in);

    if (testStrm) {
        testStrm.close();
        return true;
    } else
    {
        testStrm.open(fileName.c_str(),ios::in | ios::binary);
        if (testStrm) {
            testStrm.close();
            return true;
        }
    }

    cout<<"Can not find input file :"<<fileName<<endl;
    return false;
}
void CUnit::LoadProfiles(void)
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open profile of current subject"),"",tr("profile(*.prf);;All Files(*)"));
    prfName = fileName.toStdString();
    bool hasFile= this->VaryfyInputFile(prfName);
    if (false==hasFile)
        QMessageBox::critical(this, tr("Input NOT FOUND"), tr("Please input a valid volume file"),QMessageBox::Cancel);
    else
    {
        this->ParserProfile(prfName);
        this->VisualizeSurface( prfMap[surf_KEY]);
        this->VisualizeBOLDS( prfMap[bolds_KEY]);
        this->VisualizeFibers( prfMap[fiber_KEY]);
        this->VisualizeB0Volume( prfMap[dtiB0_KEY]);
        read_volume(myDtiGM, prfMap[dtiseg_KEY]);

        //set up joint model;

        if (NULL!=myModel)
            delete myModel;

        myModel= new CJointModel( *mySurface,*myFiber, myBolds, myDtiGM);
        delete mySurface;
        delete myFiber;
        mySurface= & ( myModel->GetSurface() );
        myFiber= & (myModel->GetFibers());

        myModel->SetOffSet(b0Offset[0],b0Offset[1],b0Offset[2]);
        myModel->FindFiberSurfIntersections();
        myModel->FindFibersOnCells();
//     // myModel->SaveIntersectionPoints("intersections.vtk");
        myModel->WarpingBoldsOnSurface();


        emit this->ProfileLoaded();
        //load network;
        this->LoadNetwork(prfName);

    }
}

void CUnit::SaveProfile(void )
{
    std::string netName=prfName;
    cout<<"saving prfName "<<netName;
    fstream netStrm;
    OpenWriteStrmAscii(netStrm, netName);
    netStrm<<"SURFACE "<<prfMap[surf_KEY]<<endl;
    netStrm<<"FIBERS "<<prfMap[fiber_KEY]<<endl;
    netStrm<<"DTIB0 "<<prfMap[dtiB0_KEY]<<endl;
    netStrm<<"BOLDS "<<prfMap[bolds_KEY]<<endl;
    netStrm<<"DTIGM "<<prfMap[dtiseg_KEY]<<endl;
// 		  netStrm<<"SURFACE "<<prfMap[surf_KEY]<<endl;
// 		  netStrm<<"SURFACE "<<prfMap[surf_KEY]<<endl;
    netStrm<<"NETWORK ";
    if (modeType==FunctionalConnectivity)
        netStrm<<"FunctionalConnectivity"<<endl;
    else
        netStrm<<"EffectiveConnectivity"<<endl;

    netStrm<<"COORDINATES"<<" "<<allNetworkNodes.size()<<endl;
    for (int roiID = 0; roiID < allNetworkNodes.size(); ++roiID) {
        netStrm<<allROINames[roiID]<<" "<<allNetworkNodes[roiID].x<<" "<<allNetworkNodes[roiID].y<<" "<<allNetworkNodes[roiID].z<<endl;
    } //end of for loop:: roiID
    netStrm<<"EDGES  CALCULATED"<<endl;
    if ( FunctionalConnectivity == this->modeType)
        this->CalculateNetEdgePC();
    else
        this->CalculateNetEdgeGC();

    for (int rowID = 0; rowID < allNetworkNodes.size(); ++rowID) {
        for (int colID = 0; colID < allNetworkNodes.size(); ++colID) {
            netStrm<< (edges)(rowID+1,colID+1)<<" ";
        } //end of for loop:: colID
        netStrm<<endl;
    } //end of for loop:: rowID
    netStrm.close();

    //output the rois;
    std::string roivtkName= prfName+".rois.vtk";
    KML::VisualizePointsUsingSphereByVTK(allNetworkNodes,roivtkName);
    cout<<"...done!"<<endl;

}

void CUnit::SaveProfileAs(void )
{

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save current profile"),"",tr("profile(*.prf);;All Files(*)"));
    std::string myFileName= fileName.toStdString();

    if (myFileName.size()) {
        std::string netName=myFileName;
        cout<<"saving..."<<netName;
        fstream netStrm;
        OpenWriteStrmAscii(netStrm, netName);
        netStrm<<"SURFACE "<<prfMap[surf_KEY]<<endl;
        netStrm<<"FIBERS "<<prfMap[fiber_KEY]<<endl;
        netStrm<<"DTIB0 "<<prfMap[dtiB0_KEY]<<endl;
        netStrm<<"BOLDS "<<prfMap[bolds_KEY]<<endl;
        netStrm<<"DTIGM "<<prfMap[dtiseg_KEY]<<endl;
// 		  netStrm<<"SURFACE "<<prfMap[surf_KEY]<<endl;
// 		  netStrm<<"SURFACE "<<prfMap[surf_KEY]<<endl;
        netStrm<<"NETWORK ";
        if (modeType==FunctionalConnectivity)
            netStrm<<"FunctionalConnectivity"<<endl;
        else
            netStrm<<"EffectiveConnectivity"<<endl;

        netStrm<<"COORDINATES"<<" "<<allNetworkNodes.size()<<endl;
        for (int roiID = 0; roiID < allNetworkNodes.size(); ++roiID) {
            netStrm<<allROINames[roiID]<<" "<<allNetworkNodes[roiID].x<<" "<<allNetworkNodes[roiID].y<<" "<<allNetworkNodes[roiID].z<<endl;
        } //end of for loop:: roiID
        netStrm<<"EDGES  CALCULATED"<<endl;
        if ( FunctionalConnectivity == this->modeType)
            this->CalculateNetEdgePC();
        else
            this->CalculateNetEdgeGC();

        for (int rowID = 0; rowID < allNetworkNodes.size(); ++rowID) {
            for (int colID = 0; colID < allNetworkNodes.size(); ++colID) {
                netStrm<< (edges)(rowID+1,colID+1)<<" ";
            } //end of for loop:: colID
            netStrm<<endl;
        } //end of for loop:: rowID
        netStrm.close();

        //output the rois;
        std::string roivtkName= prfName+".rois.vtk";
        KML::VisualizePointsUsingSphereByVTK(allNetworkNodes,roivtkName);

        cout<<"...done!"<<endl;


    }
    else
    {
        QMessageBox::critical(this, tr(""), tr("Please input a valid name for network "),QMessageBox::Cancel);
    }

}

int CUnit::ParserProfile(const std::string& myFileName)
{
    fstream prfStrm;
    OpenReadStrmAscii(prfStrm,myFileName);

    while (! prfStrm.eof() )
    {
        string tmpString;
        prfStrm>>tmpString;
        if (0==tmpString.size())
            continue;

        /*    cout<<tmpString<<endl; */
        if ( tmpString==string("DTIGM") || tmpString==string("SURFACE") || tmpString ==string("DTIB0")  || tmpString == string("BOLDS") || tmpString == string( "FIBERS") || tmpString==t1_KEY || tmpString==t1Label_KEY || tmpString == atlas_KEY)
        {
            string name;
            while (! prfStrm.eof()  &&  prfStrm>> name ) {
                if ( name.size() ) {
                    break;
                }
            };
            if ( name.size())
            {
                if (tmpString == atlas_KEY )
                {
                    prfMap[tmpString] = name;
                }
                else if ( VaryfyInputFile(name) )
                    prfMap[tmpString]=name;
            }
            else
            {
                QMessageBox::critical(this, tr("Profile Error FOUND"), tr("Please correct the profile"),QMessageBox::Cancel);
                return 1;
            }
            continue;
        }
        else if ( tmpString == string("NETWORK") )
        {
            //TODO: a method to process;
        }
        else
        {
            //TODO: a method to process;
        }
    }

    //check whether the profile is complete;


    if ( !prfMap[surf_KEY].size() || !prfMap[dtiB0_KEY].size() ||  !prfMap[bolds_KEY].size() || !prfMap[fiber_KEY].size() )
    {
        QMessageBox::critical(this, tr("Profile Error FOUND"), tr("Missing item in the profile"),QMessageBox::Cancel);
        return 1;
    }

    prfStrm.close();
    return 0;
}
int CUnit::ParserNetworkSettings(const std::string& myFileName)
{
    fstream netDefFile;
    KML::OpenReadStrmAscii(netDefFile, myFileName);
    std::string tmpString;

    //find keywork
    while ( netDefFile >> tmpString )
    {
        if (tmpString == network_KEY ) {
            tmpString.clear();
            netDefFile>>tmpString;
            if (tmpString==string("FunctionalConnectivity")) {
                modeType=FunctionalConnectivity;
//     		this->actionFunctional_Network->setChecked(true);
//     		this->actionEffective_Network->setChecked(false);
            }
            else if (tmpString==string("EffectiveConnectivity")) {
                modeType=EffectiveConnectivity;
//     		this->actionEffective_Network->setChecked(true);
//     		this->actionFunctional_Network->setChecked(false);
            } else if (tmpString==string("StructrualConnectivity"))
            {
                modeType=StructrualConnectivity;
            }
            else
                modeType=NONEType;

            break;
        }
    }

    if ( netDefFile.eof() )
    {
        QMessageBox::critical(this, tr("Invalid network"), tr("Missing network keyword in the profile"),QMessageBox::Cancel);
        return 1;
    }

    while (! netDefFile.eof()) {
        getline(netDefFile,tmpString);
        if (0==tmpString.size())
            continue;
        else
            break;
    } // now have header in eachline;

    stringstream headStrm;
    headStrm<< tmpString;
    int roiNum=0;
    string roiFormat;
    headStrm>> roiFormat>> roiNum;

    //some cleaning jobs if already loaded with one network;
    if (allNetworkNodes.size())
    {
        this->RemoveNetwork();
    }

    edges.ReSize(roiNum,roiNum);

    // now read rois;
    if (string("VOLUMES")== roiFormat) {
        string currentVolName;
        string currentROIName;
        int index=0;
        while (index < roiNum)
        {
            netDefFile>> currentROIName;
            if (0==currentROIName.size() || currentROIName[0]=='#')
            {
                getline( netDefFile, currentROIName);
                currentROIName.clear();
                continue;
            }
            allROINames.push_back(currentROIName);

            netDefFile>> currentVolName;

            Vector3D<float> theCoords;

            this->GetRoiCenterSurfCoords(currentVolName, theCoords);
            //  this->GetRoiCenterDTIVolCoords(currentVolName, theCoords);
            allNetworkNodes.push_back(theCoords);
            vector<float> firstEigenVector;
            this->GetPCA1stCurrentROI(theCoords,5,firstEigenVector);
            allNodesPCA1st.push_back(firstEigenVector);
            ++index;
        }
    }
    else if (string("SURFACEIDS")==roiFormat) {
        int currentROISurfID;
        string currentROIName;
        int index=0;
        while (index < roiNum)
        {
            netDefFile>> currentROIName;
            if (0==currentROIName.size() || currentROIName[0]=='#')
            {
                getline( netDefFile, currentROIName);
                currentROIName.clear();
                continue;
            }
            allROINames.push_back(currentROIName);

            netDefFile>> currentROISurfID;
            Vector3D<float> theCoords;
            //this->GetRoiCenterCoords(currentVolName, theCoords);
            theCoords=mySurface->GetPointCoords(currentROISurfID);
            allNetworkNodes.push_back(theCoords);
            vector<float> firstEigenVector;
            this->GetPCA1stCurrentROI(theCoords,5,firstEigenVector);
            allNodesPCA1st.push_back(firstEigenVector);
            ++index;
        }
    }
    else if (string("COORDINATES")==roiFormat) {
        string currentROIName;
        int index=0;
        while (index < roiNum)
        {
            netDefFile>> currentROIName;
            if (0==currentROIName.size() || currentROIName[0]=='#')
            {
                getline( netDefFile, currentROIName);
                currentROIName.clear();
                continue;
            }
            allROINames.push_back(currentROIName);
            Vector3D<float> theCoords;
            netDefFile>> theCoords.x>>theCoords.y>>theCoords.z;
            int theSurfID= this->GetRoiCenterSurfID(theCoords);
            allNetworkNodes.push_back(mySurface->GetPointCoords(theSurfID));
            vector<float> firstEigenVector;
            this->GetPCA1stCurrentROI(theCoords,5,firstEigenVector);
            allNodesPCA1st.push_back(firstEigenVector);
            ++index;
        }
    }
    else
    {
        QMessageBox::critical(this, tr("Invalid ROI Format"),
                              tr("The ROI Format is invalid.\n"
                                 "Valid formats are VOLUMES SURFACEIDS COORDINATES.") ,QMessageBox::Cancel);
        return 1;
    }



    //read edges;
    while (! netDefFile.eof() && netDefFile>> tmpString )
    {
        if (string("EDGES")==tmpString) {
            tmpString.clear();
            netDefFile>> tmpString;
            if (string("CALCULATED") != tmpString) {
                if ( FunctionalConnectivity == this->modeType)
                    this->CalculateNetEdgePC();
                else
                    this->CalculateNetEdgeGC();
                netDefFile.close();
                return 0;
            }

            break;
        }
    }


    for (int lineid=1; lineid<= roiNum ; ++lineid)
        for ( int colid=1; colid<= roiNum; ++colid)
        {
            netDefFile>> (edges)(lineid, colid);
        }

    netDefFile.close();
    return 0;

}
void CUnit::LoadNetwork(void )
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open network definition"),"",tr("profile(*.net);;All Files(*)"));
    std::string myFileName= fileName.toStdString();
    this->LoadNetwork(myFileName);
}

void CUnit::LoadNetwork(string myFileName) {


    if (false== VaryfyInputFile(myFileName))
    {
        QMessageBox::critical(this,tr("Invalid input file."), tr("network configure file invalid, please check!"), QMessageBox::Cancel);
    }
    else
    {
        if (0!= this->ParserNetworkSettings(myFileName))
        {
            QMessageBox::critical(this,tr("Warning!."), tr("No Network Configure Found!"), QMessageBox::Cancel);
        }
        else
        {
            isSingleROIMode=false;
            this->VisualizeNetwork();
            allROIAttributeStatus= vector<ROIAttrStatus> (allNetworkNodes.size(), NotDone);

            //for sync template;
            allROIFIbersVTK.resize(allNetworkNodes.size());
            allROIFiberIDs.resize(allNetworkNodes.size());
            this->GenerateNetworkControllorAdjs(allNetworkNodesActors);
            this->GenerateNetworkROIPCAs();

            //visulize the bold for the first ROI;
//               int initROI= allNetworkNodes.size()/2;
// 	      allroiAdjs[initROI]->checkBox_BOLD->setChecked(true);
// 	      this->UpdateCheckBoxStatus(allroiAdjs[initROI]);
// 	      this->ChangeROIBOLDSVisiability(initROI);


            //for tracemap ;
            traceMap.SetFiber(myFiber);

            //network added;
            emit this->NetworkAdded();
        }
    }

}
void CUnit::GenerateNetworkControllorAdjs(vector< vtkSmartPointer< vtkActor > >& allROIActors)
{
    if (NULL== roiListArea)
    {

        cout<<"create roilist..."<<endl;
        roiListArea = new QScrollArea;
//	roiListArea->setAttribute(Qt::WA_DeleteOnClose);
        roiListArea->setBackgroundRole(QPalette::Light);
        roiListArea->setWindowTitle(tr("Network Controller"));
        QWidget *window = new QWidget;
        QVBoxLayout *layout = new QVBoxLayout;

        allroiAdjs.resize(allROIActors.size());
        for (int index=0; index < allROIActors.size(); ++index)
        {
            myROIAdj* theROI= new myROIAdj;
            layout->addWidget(theROI);
            allroiAdjs[index]=theROI;
        }

        // must insert first ? ??
        window->setLayout(layout);
        roiListArea->setWidget(window);

        for (int index=0; index< allROIActors.size(); ++index)
        {
            const RGBTYPE& theColor=myColorScheme.GetColorByIndex(index);
            QColor qColorFG;
            qColorFG.setRgbF(theColor.x, theColor.y, theColor.z,1);

            QPalette palette;
            palette.setColor(allroiAdjs[index]->roiIndex->foregroundRole(), qColorFG/*Qt::yellow*/);
            allroiAdjs[index]->roiIndex->setPalette(palette);
            allroiAdjs[index]->roiIndex->setAutoFillBackground(true);

            Vector3D<float>& theCoords= allNetworkNodes[index];
            allroiAdjs[index]->SetROIName(allROINames[index]);
            allroiAdjs[index]->SetROIID(index);
            allroiAdjs[index]->SetCoordinates(theCoords.x, theCoords.y, theCoords.z);

            connect( allroiAdjs[index], SIGNAL(RoiCheckBoxShowChange(int )), this, SLOT(ChangeROISphereVisiability(int)));
            connect( allroiAdjs[index], SIGNAL(RoiCheckBoxFiberChange(int)),this, SLOT(ChangeROIFIberVisiability(int)));
            connect( allroiAdjs[index], SIGNAL(RoiCheckBoxBOLDChange(int)),this, SLOT(ChangeROIBOLDSVisiability(int)));
            connect( allroiAdjs[index], SIGNAL(RoiRadiusChanged(int)), this, SLOT(ChangeROIRadius(int)));

        }//

        this->InitializeCheckBoxesStatusFromNetworkControllor();

    }

}
void CUnit::InitializeCheckBoxesStatusFromNetworkControllor(void )
{
    //update status size;
    if (roiCheckStatus.size()!= allNetworkNodes.size())
    {
        roiCheckStatus.clear();
        fiberCheckStatus.clear();
        boldCheckStatus.clear();
        roiCheckStatus.resize(allNetworkNodes.size(),true);
        fiberCheckStatus.resize(allNetworkNodes.size(),false);
        boldCheckStatus.resize(allNetworkNodes.size(),false);
    }
}

void CUnit::UpdateCheckBoxStatus(myROIAdj* theAdj)
{

    int index=theAdj->GetROIID();
    cout<<"UpdateCheckBoxStatus:"<<index<<endl;
    this->InitializeCheckBoxesStatusFromNetworkControllor();

    for (int index=0; index<  allNetworkNodes.size(); ++index)
    {

//   cout<<roiCheckStatus<<endl;
        roiCheckStatus[index]= allroiAdjs[index]->checkBox->isChecked();
//   cout<<roiCheckStatus<<endl;
//   cout<<fiberCheckStatus<<endl;
        fiberCheckStatus[index]= allroiAdjs[index]->checkBox_Fiber->isChecked();
//   cout<<fiberCheckStatus<<endl;
//   cout<<boldCheckStatus<<endl;
        boldCheckStatus[index]= allroiAdjs[index]->checkBox_BOLD->isChecked();
//   cout<<boldCheckStatus<<endl;
//
    }
    //emit this->ROIAdjCheckBoxChanged(index);
}

void CUnit::GenerateNetworkROIPCAs ( void )
{
    // warp bolds to the cortical surface

    allNodesPCA1st.resize(allNetworkNodes.size());
    float roiSize=5;
    if (allroiAdjs.size())
        roiSize= allroiAdjs[0]->GetRadius();

    for (int roiID=0; roiID< allNetworkNodes.size(); ++roiID)
    {
        //get all the surface IDs in the roi;
        cout<<"generating pca for node: "<< roiID<<endl;;
        this->GetPCA1stCurrentROI(roiID,roiSize,allNodesPCA1st[roiID]);
    }

}
void CUnit::ShowROIAdjs(void )
{
    if (NULL==roiListArea)
        this->GenerateNetworkControllorAdjs( allNetworkNodesActors);
    roiListArea->show();
//   this->ChangeROIBOLDSVisiability(0);
}
void CUnit::ChangeROICoords(myROIAdj& theROI)
{
    int id= theROI.GetROIID();

    //update loations;
    allNetworkNodes[id].x= theROI.Getx();
    allNetworkNodes[id].y= theROI.Gety();
    allNetworkNodes[id].z= theROI.Getz();

    //update the pca1sts;
    //cout<<allNodesPCA1st[id]<<endl;
    this->GetPCA1stCurrentROI(id, theROI.GetRadius(), allNodesPCA1st[id]);
    //cout<<allNodesPCA1st[id]<<endl;

    //update edge and weight;
    if (modeType==FunctionalConnectivity)
        this->UpdateFunctionalNetworkWeights(id);
    else
        this->UpdateEffectiveNetworkWeights(id);

    //update the visualization;
    this->UpdateNetworkVisualization(theROI);


    //change fibers;
    allROIAttributeStatus[id] = NeedUpdate;
    if ( theROI.checkBox_Fiber->isChecked() )
        this->ChangeROIFIberVisiability(theROI);

    //change bols;
    this->UpdateBOLDSignal4Node(id,theROI.GetRadius());

}
void CUnit::ChangeROISphereVisiability(myROIAdj& theROI)
{
    int id= theROI.GetROIID();
    vtkSmartPointer<vtkActor> theActor= allNetworkNodesActors[id];
    if ( /*theROI.checkBox->isChecked()*/ roiCheckStatus[id])
    {
        theActor->GetProperty()->SetOpacity(1);
        //show edges;
        for (int index=0; index < allNodesLinksActors[id].size(); ++index)
        {
            allNodesLinksActors[id][index]->GetProperty()->SetOpacity(1);
        }
    }
    else {
        theActor->GetProperty()->SetOpacity(0);
        //disable edges;
        for (int index=0; index < allNodesLinksActors[id].size(); ++index)
        {
            allNodesLinksActors[id][index]->GetProperty()->SetOpacity(0);
        }
    }
    fiberWidget->update();
    surfWidget->update();
    volXWidget->update();
    volYWidget->update();
    volZWidget->update();

}
void CUnit::UpdateNetworkVisualization(myROIAdj& theROI)
{

    int roiId= theROI.GetROIID();
    //change location;
    allSphereSources[roiId]->SetCenter( allNetworkNodes[roiId].x,allNetworkNodes[roiId].y,allNetworkNodes[roiId].z);


    //remove all realated arrows; and add them again;
    if ( allNodesLinksActors.size())
    {
        for (int index=0; index < allNodesLinksActors[roiId].size(); ++index)
        {
            fiberRenderer->RemoveActor(allNodesLinksActors[roiId][index]);
            allNodesLinksActors[roiId][index]->GetMapper()->RemoveAllInputs();
        }

        allNodesLinksActors[roiId].clear();
    }


    //now add the edges;
    for (int line=1; line<= edges.Nrows(); ++line)
    {
        if (modeType==EffectiveConnectivity) {
            this->AddArrow(line-1, roiId, allNetworkNodes, edges(line,roiId+1),theROI.GetRadius());
            this->AddArrow(roiId,line-1 , allNetworkNodes, edges(roiId+1,line),theROI.GetRadius());
        }
        if (modeType==FunctionalConnectivity)
        {
            this->AddLine(line-1, roiId, allNetworkNodes, edges(line,roiId+1),theROI.GetRadius());
            this->AddLine(roiId,line-1 , allNetworkNodes, edges(roiId+1,line),theROI.GetRadius());
        }
    }

    surfWidget->update();
    fiberWidget->update();
    volXWidget->update();
    volYWidget->update();
    volZWidget->update();
}
void CUnit::ChangeROIBOLDSVisiability(myROIAdj& theROI)
{

//   cout<<theROI.GetROIID()<<endl;
    boldSignalChart->ClearPlots();
    cout<<"Plotting bolds..";
    /*
    if( pTemplateWnd  && allNetworkNodes.size() )
    {
      int index= theROI.GetROIID();

      if( boldCheckStatus[index])
      {
      cout<<index<<" ";

      linePlot = boldSignalChart->AddPlot(vtkChart::LINE);
      linePlot->SetInput(boldSignalTable, 0, index+1);
      linePlot->SetColor(myColorScheme.GetColorByIndex(index).x*255,myColorScheme.GetColorByIndex(index).y*255,myColorScheme.GetColorByIndex(index).z*255,255);
      linePlot->SetWidth(2);
      linePlot->GetXAxis()->SetMaximum(myBolds.tsize()+2);
      linePlot->GetXAxis()->SetMinimum(0);

      linePlot->GetYAxis()->SetMaximum(4);
      linePlot->GetYAxis()->SetMinimum(-4);
      linePlot->GetYAxis()->AutoScale();
      }

    }else*/
    {
        for (int index=0; index< allNetworkNodes.size(); ++index)
        {
            if ( boldCheckStatus[index]  /*allroiAdjs[index]->checkBox_BOLD->isChecked()*/)
            {

                cout<<index<<" ";

                linePlot = boldSignalChart->AddPlot(vtkChart::LINE);
                linePlot->SetInput(boldSignalTable, 0, index+1);
                linePlot->SetColor(myColorScheme.GetColorByIndex(index).x*255,myColorScheme.GetColorByIndex(index).y*255,myColorScheme.GetColorByIndex(index).z*255,255);
                linePlot->SetWidth(2);
                linePlot->GetXAxis()->SetMaximum(myBolds.tsize()+2);
                linePlot->GetXAxis()->SetMinimum(0);

                linePlot->GetYAxis()->SetMaximum(4);
                linePlot->GetYAxis()->SetMinimum(-4);
                linePlot->GetYAxis()->AutoScale();
            }

        }
    }

    cout<<endl;
    boldSignalChart->Update();
    boldSignalTable->Update();
    BoldWidget->update();
}
void CUnit::ChangeROIFIberVisiability(myROIAdj& theROI)
{

    int roiID= theROI.GetROIID();

    //now show all checked fibers;
    VTK_CREATE(vtkCellArray,myArray);
    myArray->Allocate(5000);
    cout<<"ploting fibers: ";

    /*
      if( pTemplateWnd  && allNetworkNodes.size() )
      {
        int index= theROI.GetROIID();
        if( fiberCheckStatus[index])
        {
          cout<<index<<" ";
          for(int findex= 0; findex < allROIFIbersVTK[index].size(); ++ findex)
          {
    	 myArray->InsertNextCell(allROIFIbersVTK[index][findex]);
          }
        }
      }
      else*/
    {
        for (int index=0; index< allNetworkNodes.size(); ++index)
        {

            if (allROIAttributeStatus[index]!= Updated)
            {
                cout<<"tracking..."<<index<<endl;
                TrackFibersCurrentROI(*allroiAdjs[index]);
            }
            if ( fiberCheckStatus[index])
            {
                cout<<index<<" ";
                for (int findex= 0; findex < allROIFIbersVTK[index].size(); ++ findex)
                {
                    myArray->InsertNextCell(allROIFIbersVTK[index][findex]);
                }
            }
        }
    }

    cout<<endl;

    fiberData->SetLines(myArray);
    fiberRenderer->AddActor(fiberActor);
    fiberWidget->update();
    emit this->ROIFiberChanged(roiID);

}
void CUnit::TrackFibersCurrentROI(myROIAdj& theROI)
{

    int roiID= theROI.GetROIID();
    const KML::Vector3D< float >& centerCoords= allNetworkNodes[roiID];
    float radius = theROI.GetRadius();

    vector<vtkSmartPointer<vtkIdList> >  currentroiFibers;
    vector<int> currentROIIds;
    int fiberID=0;

    while (fiberID<myFiber->GetNumFibers())
    {
        bool rightFiber=false;

        vector<int>& currentFiber=myFiber->GetFiber(fiberID);
        for (int elementID=0; elementID< currentFiber.size(); ++elementID)
        {
            int pid= currentFiber[elementID];
            PointCoordType& pcorrd= myFiber->GetPointCoords(pid);
            PointCoordType tmpPoint=pcorrd-centerCoords;
            if (radius> tmpPoint.Norm())
            {
                rightFiber=true;
                break;
            }
            ++elementID;
        }
        if (rightFiber)
        {
            currentROIIds.push_back(fiberID);
            vtkSmartPointer<vtkIdList>  myCellList= vtkSmartPointer<vtkIdList>::New();
            for (int elementID=0; elementID< currentFiber.size(); ++elementID)
            {
                int pid= currentFiber[elementID];
                myCellList->InsertNextId(pid);
            }
            currentroiFibers.push_back(myCellList);
        }


        ++fiberID;
    } //end of while loop;

    allROIFIbersVTK[roiID]= currentroiFibers;
    allROIFiberIDs[roiID]= currentROIIds;
    allROIAttributeStatus[roiID]= Updated;
}
void CUnit::ChangeROIRadius(myROIAdj& theROI)
{

    //update sphere;
    allSphereSources[theROI.GetROIID()]->SetRadius(theROI.GetRadius());
    surfWidget->update();
    UpdateNetworkVisualization(theROI);

    //update fiber ;
    allROIAttributeStatus[theROI.GetROIID()] = NeedUpdate;
    if ( theROI.checkBox_Fiber->isChecked() )
        this->ChangeROIFIberVisiability(theROI);

    //update bolds;
    this->UpdateBOLDSignal4Node(theROI.GetROIID(),theROI.GetRadius());
}
int CUnit::GetROIID(vtkActor* theActor)
{
    for (int index=0; index < allSphereSources.size(); ++index) {
        if (theActor==  allNetworkNodesActors[index]) {
            return index;
        }
    }
    cout<<"ROI Not found";
    return -1; //error information;
}
void CUnit::ReleaseMouseMove(void )
{
    //reset pick status for the roi picked;
    vtkRenderWindowInteractor* rwi= surfWidget->GetInteractor();
    if ( rwi->GetKeyCode()=='s')
    {
        for (int roiID = 0; roiID < allNetworkNodesActors.size(); ++roiID) {
            allNetworkNodesActors[roiID]->PickableOn();
            if (! isSingleROIMode)
                this->UpdateNetworkVisualization(*allroiAdjs[roiID]);
        } //end of for loop:: roiID

    }

}
void CUnit::ChangeROIBOLDSVisiability(int roiID)
{
    cout<<"ChangeROIBOLDSVisiability: "<<roiID<<endl;
    myROIAdj* theAdj= static_cast<myROIAdj*> (QObject::sender());
    this->UpdateCheckBoxStatus(theAdj);
    this->ChangeROIBOLDSVisiability(*allroiAdjs[roiID]);
    emit this->ROIAdjCheckBoxChanged(roiID);
}
void CUnit::ChangeROICoords(int roiID)
{
    this->ChangeROICoords(*allroiAdjs[roiID]);
}
void CUnit::ChangeROIFIberVisiability(int roiID)
{
    myROIAdj* theAdj= static_cast<myROIAdj*> (QObject::sender());
    this->UpdateCheckBoxStatus(theAdj);
    this->ChangeROIFIberVisiability(*allroiAdjs[roiID]);
    emit this->ROIAdjCheckBoxChanged(roiID);
}
void CUnit::ChangeROIRadius(int roiID)
{
    this->ChangeROIRadius(*allroiAdjs[roiID]);
}
void CUnit::ChangeROISphereVisiability(int roiID)
{

    myROIAdj* theAdj= static_cast<myROIAdj*> (QObject::sender());
    this->UpdateCheckBoxStatus(theAdj);
    this->ChangeROISphereVisiability(*allroiAdjs[roiID]);
    emit this->ROIAdjCheckBoxChanged(roiID);
}

int CUnit::Change2EffectiveConnectivity(void )
{
//   if(actionEffective_Network->isChecked())
//     return 0;

    modeType=EffectiveConnectivity;

    if (edges.Ncols()) {
        this->CalculateNetEdgeGC();
        this->RemoveNetworkLinksActors();
        this->AddNetworkLinksActors();
        fiberWidget->update();
    }

    if (NULL!=pTemplateWnd)
        pTemplateWnd->Change2EffectiveConnectivity();

    return 0;
}

void CUnit::RemoveNetworkLinksActors(void )
{
    //remove all cones and lines;
    for (int roiId=0; roiId < allNetworkNodes.size(); ++roiId)
    {
        if (allNodesLinksActors.size())
        {
            for (int index=0; index < allNodesLinksActors[roiId].size(); ++index)
            {
                fiberRenderer->RemoveActor(allNodesLinksActors[roiId][index]);
                allNodesLinksActors[roiId][index]->GetMapper()->RemoveAllInputs();
            }

            allNodesLinksActors[roiId].clear();

        }
    }

    allNodesLinksActors.clear();

}


void CUnit::AddNetworkLinksActors(void )
{
    allNodesLinksActors.clear();
    allNodesLinksActors.resize(allNetworkNodes.size());

    switch (modeType) {
    case NONEType :
        break;
    case FunctionalConnectivity :
    {
        // add new lines;
        for (int lineID=1; lineID <=  allNetworkNodes.size(); ++lineID)
            for (int colID=lineID; colID <= allNetworkNodes.size(); ++colID)
            {
                if (lineID!=colID && abs((edges)(lineID,colID)) > 0.3 )	{
                    this->AddLine(lineID-1, colID-1, allNetworkNodes,((edges)(lineID,colID)));
                }
            }
        break;
    }
    case EffectiveConnectivity :
    {
        for (int lineID=1; lineID <=  allNetworkNodes.size(); ++lineID)
            for (int colID=lineID; colID <= allNetworkNodes.size(); ++colID)
            {
                if (lineID!=colID && abs((edges)(lineID,colID)) > 0.03 )
                {
                    this->AddArrow(lineID-1, colID-1, allNetworkNodes,abs((edges)(lineID,colID)));
                }
            }

        break;
    }
    case StructrualConnectivity :
    {

        break;
    }
    default :
    {
        QMessageBox::critical(this,tr("Invalid Network Type."), tr("network type invalid, please check!"), QMessageBox::Cancel);

    }//end of default;

    }

}


int CUnit::Change2FunctionalConnectivity(void )
{
//   if(actionFunctional_Network->isChecked())
//     return 0;

    modeType=FunctionalConnectivity;

    if ( edges.Ncols()) {

        //calculate the correlation matrix;

        this->CalculateNetEdgePC();
        this->RemoveNetworkLinksActors();
        this->AddNetworkLinksActors();


        fiberWidget->update();
    }

    if (NULL!=pTemplateWnd)
        pTemplateWnd->Change2FunctionalConnectivity();

    return 0;
}

int CUnit::Change2StructuralConnectivity(void )
{
    modeType=StructrualConnectivity;
    if (NULL!=pTemplateWnd)
        pTemplateWnd->Change2StructuralConnectivity();
    return 0;
}


int CUnit::Change2NoneConnectivity(void )
{
    modeType=NONEType;
    this->RemoveNetworkLinksActors();
    fiberWidget->update();
    surfWidget->update();

    if (NULL!=pTemplateWnd)
        pTemplateWnd->Change2NoneConnectivity();
    return 0;
}
void CUnit::AddNetworkNodesActors(void )
{
    allNetworkNodesActors.clear();
    float roiSize=5;
    if (allroiAdjs.size())
        roiSize= allroiAdjs[0]->GetRadius();


    for (int index=0; index< allNetworkNodes.size(); ++index)
    {
        const Vector3D<float>& surfCoords= allNetworkNodes[index];
        VTK_CREATE(vtkSphereSource,sphere);
        sphere->SetCenter(surfCoords.x,surfCoords.y,surfCoords.z);
        sphere->SetRadius(roiSize);
        sphere->SetPhiResolution(30);
        sphere->SetThetaResolution(30);
        allSphereSources.push_back(sphere);

        VTK_CREATE(vtkPolyDataMapper,mapper);
        mapper->SetInput(sphere->GetOutput());
        VTK_CREATE(vtkActor,actor);
        actor->SetMapper(mapper);
        const RGBTYPE& theColor=myColorScheme.GetColorByIndex(index);
        actor->GetProperty()->SetColor(theColor.x,theColor.y, theColor.z);
        fiberRenderer->AddActor(actor);
        surfRenderer->AddActor(actor);
        allNetworkNodesActors.push_back(actor);
    }

}

void CUnit::RemoveNetworkNodesActors(void )
{

    for (int roiId=0; roiId < allNetworkNodesActors.size(); ++roiId)
    {
        fiberRenderer->RemoveActor(allNetworkNodesActors[roiId]);
        surfRenderer->RemoveActor(allNetworkNodesActors[roiId]);
        volViewX->GetRenderer()->RemoveActor(allNetworkNodesActors[roiId]);
        volViewY->GetRenderer()->RemoveActor(allNetworkNodesActors[roiId]);
        volViewZ->GetRenderer()->RemoveActor(allNetworkNodesActors[roiId]);

        //TODO::BODES view change;
    }

    allNetworkNodesActors.clear();
}

void CUnit::RemoveNetwork(void )
{

    cout<<"clean previous rois"<<endl;
    //remove picker if necessary;

    this->RemoveNetworkLinksActors();
    this->RemoveNetworkNodesActors();
    fiberRenderer->AddActor(fiberActor);

    surfWidget->update();
    fiberWidget->update();
    volViewZ->UpdateDisplayExtent();;
    volViewY->UpdateDisplayExtent();
    volViewX->UpdateDisplayExtent();

    allROINames.clear();
    allNetworkNodes.clear();
    allNodesPCA1st.clear();
    allSphereSources.clear();
    allROIAttributeStatus.clear();
    fiberCheckStatus.clear();
    roiCheckStatus.clear();
    boldCheckStatus.clear();
    allROIFIbersVTK.clear();
    allROIFiberIDs.clear();
    allroiAdjs.clear();
    MyDestroy(roiListArea);
    roiListArea=NULL;

    isSingleROIMode=true;
    emit this->NetworkRemoved();
}


void CUnit::UpdateEffectiveNetworkWeights ( int roiID )
{
    //update related weights;
    for (int lineID=1; lineID<=  allNodesPCA1st.size(); ++lineID)
        if (lineID!=roiID+1)
        {
            this->GetGCStrength( allNodesPCA1st[ lineID-1], allNodesPCA1st[roiID],  (edges)(lineID,roiID+1), (edges)(roiID+1, lineID));
        }
        else
        {
            (edges)(roiID+1, lineID)=0;
        }
}
void CUnit::UpdateFunctionalNetworkWeights ( int roiID )
{
    //update related weights;
    for (int lineID=1; lineID<=  allNodesPCA1st.size(); ++lineID)
        if (lineID!=roiID+1)
        {
            this->GetPearsonCorrelation( allNodesPCA1st[ lineID-1], allNodesPCA1st[roiID],(edges)(lineID,roiID+1));
            (edges)(roiID+1, lineID) =  (edges)(lineID,roiID+1) ;
        }
        else
        {
            (edges)(roiID+1, lineID)=0;
        }
}
void CUnit::GetGCStrength(std::vector< float >& x, std::vector< float >& y, double& x2y, double& y2x)
{
    vector<float> gcRst=KML::GCA(x, y, 1);
    x2y=gcRst[2];
    y2x=gcRst[3];

}
void CUnit::GetPearsonCorrelation(std::vector< float >& x, std::vector< float >& y, double& coorel)
{
    coorel= CoorelOfTwoSeries(x,y,0,x.size());
}
void CUnit::SaveNetwork() {
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save current network"),"",tr("network(*.net);;All Files(*)"));
    std::string myFileName= fileName.toStdString();

    if (myFileName.size()) {
        std::string netName=myFileName+".net";
        cout<<"saving..."<<netName;
        fstream netStrm;
        OpenWriteStrmAscii(netStrm, netName);
        netStrm<<"NETWORK ";
        if (modeType==FunctionalConnectivity)
            netStrm<<"FunctionalConnectivity"<<endl;
        else
            netStrm<<"EffectiveConnectivity"<<endl;

        netStrm<<"COORDINATES"<<" "<<allNetworkNodes.size()<<endl;
        for (int roiID = 0; roiID < allNetworkNodes.size(); ++roiID) {
            netStrm<<allROINames[roiID]<<" "<<allNetworkNodes[roiID].x<<" "<<allNetworkNodes[roiID].y<<" "<<allNetworkNodes[roiID].z<<endl;
        } //end of for loop:: roiID
        netStrm<<"EDGES  CALCULATED"<<endl;
        if ( FunctionalConnectivity == this->modeType)
            this->CalculateNetEdgePC();
        else
            this->CalculateNetEdgeGC();

        for (int rowID = 0; rowID < allNetworkNodes.size(); ++rowID) {
            for (int colID = 0; colID < allNetworkNodes.size(); ++colID) {
                netStrm<< (edges)(rowID+1,colID+1)<<" ";
            } //end of for loop:: colID
            netStrm<<endl;
        } //end of for loop:: rowID
        netStrm.close();

        //output the rois;
        std::string roivtkName= myFileName+".vtk";
        KML::VisualizePointsUsingSphereByVTK(allNetworkNodes,roivtkName);

        cout<<"...done!"<<endl;


    }
    else
    {
        QMessageBox::critical(this, tr(""), tr("Please input a valid name for network "),QMessageBox::Cancel);
    }
}
void CUnit::UpdateVolViewSlices(int roiID) {

    if (allNetworkNodesActors[roiID] != roiPickedActor )
    {
        volViewX->GetRenderer()->RemoveActor(roiPickedActor);
        volViewY->GetRenderer()->RemoveActor(roiPickedActor);
        volViewZ->GetRenderer()->RemoveActor(roiPickedActor);
        roiPickedActor=allNetworkNodesActors[roiID] ;
    }

    volViewX->GetRenderer()->AddActor(allNetworkNodesActors[roiID]);
    volViewY->GetRenderer()->AddActor(allNetworkNodesActors[roiID]);
    volViewZ->GetRenderer()->AddActor(allNetworkNodesActors[roiID]);

    Vector3D<float>& coord= allNetworkNodes[roiID];

    int xSlice,ySlice,zSlice;

    if (isVolT1)
    {
        xSlice= (coord.x- t1Offset[0])/t1Spacing[0]+0.5;
        ySlice= (coord.y- t1Offset[1])/t1Spacing[1]+0.5;
        zSlice= (coord.z- t1Offset[2])/t1Spacing[2]+0.5;
    }
    else
    {
        xSlice= (coord.x- b0Offset[0])/b0Spacing[0]+0.5;
        ySlice= (coord.y- b0Offset[1])/b0Spacing[1]+0.5;
        zSlice= (coord.z- b0Offset[2])/b0Spacing[2]+0.5;
    }


    this->horizontalSliderX->setValue(xSlice);
    this->horizontalSliderY->setValue(ySlice);
    this->horizontalSliderZ->setValue(zSlice);

    volViewX->SetSlice(xSlice);
    volViewX->UpdateDisplayExtent();
    volViewY->SetSlice(ySlice);
    volViewY->UpdateDisplayExtent();
    volViewZ->SetSlice(zSlice);
    volViewZ->UpdateDisplayExtent();
}
void CUnit::GenerateNetworkViaDefinedROIS() {

    string myFileName("_definedNetwork.net");
    fstream netStrm;
    OpenWriteStrmAscii(netStrm, myFileName);
    netStrm<<"NETWORK ";
    if (modeType==FunctionalConnectivity) {
        netStrm<<"FunctionalConnectivity"<<endl;
        this->CalculateNetEdgePC();
    }
    else {
        netStrm<<"EffectiveConnectivity"<<endl;
        this->CalculateNetEdgeGC();
    }

    netStrm<<"COORDINATES"<<" "<<allNetworkNodes.size()<<endl;
    for (int roiID = 0; roiID < allNetworkNodes.size(); ++roiID) {
        netStrm<<allROINames[roiID]<<" "<<allNetworkNodes[roiID].x<<" "<<allNetworkNodes[roiID].y<<" "<<allNetworkNodes[roiID].z<<endl;
    } //end of for loop:: roiID
    netStrm<<"EDGES"<<endl;
    for (int rowID = 0; rowID < allNetworkNodes.size(); ++rowID) {
        for (int colID = 0; colID < allNetworkNodes.size(); ++colID) {
            netStrm<< (edges)(rowID+1,colID+1)<<" ";
        } //end of for loop:: colID
        netStrm<<endl;
    } //end of for loop:: rowID
    netStrm.close();

    this->LoadNetwork(myFileName);

}

void CUnit::CalculateNetEdgeGC(void) {
    edges.ReSize(allNetworkNodes.size(),allNetworkNodes.size());
    //calculate the weights;
    int netSize=allNetworkNodes.size();
    for (int lineID=1; lineID<=  netSize; ++lineID)
        for (int colID=1; colID<=netSize; ++colID)
        {
            if (lineID!=colID)
                this->GetGCStrength( allNodesPCA1st[ lineID-1], allNodesPCA1st[colID-1],  (edges)(lineID,colID), (edges)(colID, lineID));
            else
                (edges)(lineID,colID)=0;
        }
}

void CUnit::CalculateNetEdgePC(void) {
    edges.ReSize(allNetworkNodes.size(),allNetworkNodes.size());
    int netSize=allNetworkNodes.size();
    for (int lineID=1; lineID<=  netSize; ++lineID)
        for (int colID=lineID; colID<=netSize; ++colID)
        {
            if (lineID!=colID)
            {
                this->GetPearsonCorrelation( allNodesPCA1st[ lineID-1], allNodesPCA1st[colID-1],  (edges)(lineID,colID));
                (edges)(colID, lineID) =  (edges)(lineID,colID) ;
            }
            else
                (edges)(lineID,colID)=0;
        }
}

void CUnit::ChangeOpacityFiber(float opacity)
{
    fiberActor->GetProperty()->SetOpacity(opacity);
    fiberWidget->update();
    surfWidget->update();
}
void CUnit::ChangeOpacitySurface(float opacity)
{
    surfActor->GetProperty()->SetOpacity(opacity);
    surfWidget->update();
    fiberWidget->update();
}

void CUnit::SetTemplate(CUnit* ptrTemplate)
{
    pTemplateWnd= ptrTemplate;
}
void CUnit::UpdateText(int pointPicked, KML::Vector3D< double > worldCoordinate )
{

    this->ChangeROIFIberVisiability(*allroiAdjs[roiPicked]);
    pTemplateWnd->ChangeROIFIberVisiability(*(pTemplateWnd->allroiAdjs[roiPicked]));
    //now get the distance; and project it onto the screen;
    float theDisance = this->CalculateDistanceHTraceMap(pointPicked);
    QString disText = QString::number(theDisance);
    cout<<"Bundle Distance:"<<theDisance<<" "<<disText.toStdString()<<endl;;


    string disTextStd = disText.toStdString();
    string disTextStd4= disTextStd.substr(0,4);
    textActor = vtkSmartPointer<vtkTextActor>::New();
    textActor->SetInput(disTextStd4.c_str());
    vtkTextProperty* tprop = textActor->GetTextProperty();
    tprop->SetFontFamilyToArial ();
    tprop->BoldOn();
    tprop->SetLineSpacing(1.0);
    tprop->SetFontSize(24);
    tprop->SetColor(1.0,0.0,0.0);
    textActor->SetDisplayPosition( 20, 20 );

    fiberRenderer->AddActor2D(textActor);

    surfWidget->update();
    fiberWidget->update();
}

float CUnit::CalculateDistanceHTraceMap(int pointPicked)
{
    //set fiber bundles;
//   this->traceMap.SetSegLength(3);
//   this->traceMap.SetAngleStep(60);
    this->traceMap.SetBundle(allROIFiberIDs[pointPicked]);
    this->traceMap.SetSeedPoint(allNetworkNodes[pointPicked]);
//   cout<<allROIFiberIDs[pointPicked]<<endl;

    pTemplateWnd->traceMap.SetBundle(pTemplateWnd->allROIFiberIDs[pointPicked]);
    pTemplateWnd->traceMap.SetSeedPoint(pTemplateWnd->allNetworkNodes[pointPicked]);
//   cout<<pTemplateWnd->allROIFiberIDs[pointPicked]<<endl;
//   pTemplateWnd->traceMap.SetSegLength(3);
//   pTemplateWnd->traceMap.SetAngleStep(60);
    this->traceMap.DoTraceMap();
    pTemplateWnd->traceMap.DoTraceMap();

    const vector<float>& feature1= traceMap.GetTraceMapFeatures();
    const vector<float>& feature2= pTemplateWnd->traceMap.GetTraceMapFeatures();

    return KML::L1Distance<float>(feature1,feature2);
}
void CUnit::EnableBOLDs(bool option)
{
    isBOLDsEnabled=option;
}

void CUnit::Switch2DTIVol(void )
{
    this->VisualizeT1Volume(prfMap[dtiB0_KEY]);
}
void CUnit::Switch2T1Vol(bool checked )
{
    if (checked)
        this->VisualizeT1Volume(prfMap[t1_KEY]);
    else
        this->Switch2DTIVol();
}

void CUnit::UpdateAtlasInfoAtStatusBar(int roiPicked)
{
    if (prfMap[atlas_KEY].size()  &&  prfMap[t1Label_KEY].size())
    {

        if (0==atlasImg.xsize())
        {
            read_volume(atlasImg,prfMap[t1Label_KEY]);
            this->ReadMetaInfo4T1();
        }
        //determine label;

        int label=0;
        int x,y,z;
        x= (allNetworkNodes[roiPicked].x-t1Offset[0])/t1Spacing[0]+0.5;
        y= (allNetworkNodes[roiPicked].y-t1Offset[1])/t1Spacing[1]+0.5;
        z= (allNetworkNodes[roiPicked].z-t1Offset[2])/t1Spacing[2]+0.5;
        label= atlasImg(x,y,z);

        emit UpdateStatusBar(prfMap[atlas_KEY],label);
    } //doing nothing if atlas not set;

}


void CUnit::ReadMetaInfo4Dti(void )
{
    VectorType tmpVector;
    KML::GetSizesFromMHD(prfMap[dtiB0_KEY].c_str(),tmpVector);
    b0Dimensions[0]=tmpVector.x;
    b0Dimensions[1]=tmpVector.y;
    b0Dimensions[2]=tmpVector.z;


    KML::GetOffsetFromMHD(prfMap[dtiB0_KEY].c_str(),tmpVector);
    b0Offset[0]=tmpVector.x;
    b0Offset[1]=tmpVector.y;
    b0Offset[2]=tmpVector.z;


    KML::GetDimsFromMHD(prfMap[dtiB0_KEY].c_str(),tmpVector);
    b0Spacing[0]=tmpVector.x;
    b0Spacing[1]=tmpVector.y;
    b0Spacing[2]=tmpVector.z;

}
void CUnit::ReadMetaInfo4T1(void )
{
    VectorType tmpVector;
    KML::GetSizesFromMHD(prfMap[t1_KEY].c_str(),tmpVector);
    t1Dimensions[0]=tmpVector.x;
    t1Dimensions[1]=tmpVector.y;
    t1Dimensions[2]=tmpVector.z;


    KML::GetOffsetFromMHD(prfMap[t1_KEY].c_str(),tmpVector);
    t1Offset[0]=tmpVector.x;
    t1Offset[1]=tmpVector.y;
    t1Offset[2]=tmpVector.z;


    KML::GetDimsFromMHD(prfMap[t1_KEY].c_str(),tmpVector);
    t1Spacing[0]=tmpVector.x;
    t1Spacing[1]=tmpVector.y;
    t1Spacing[2]=tmpVector.z;
}

