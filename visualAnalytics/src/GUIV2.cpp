/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/
#include "GUIV2.h"
#include <algorithm>
#include "unit.h"
#include "vtkTextSource.h"
#include "vtkTextActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include <iostream>
using namespace std;
using namespace KML;


GUIV2::GUIV2()
{

    speedAdj= new CQValueChanger(0,1000);
    speedAdj->setToolTip(QString::fromUtf8("Adjust Animation Speed (1ms)"));
    opacityFiberAdj= new CQValueChanger(0,1);
    opacityFiberAdj->setToolTip(QString::fromUtf8("Adjust Fiber Opacity (1/1000)"));
    opacitySurfaceAdj= new CQValueChanger(0,1);
    opacitySurfaceAdj->setToolTip(QString::fromUtf8("Adjust Surface Opacity (1/1000)"));

    templateWindow=NULL;



    this->setupUi(this);
    this->InitializeWindow();



    //value changes;
    connect(speedAdj->horizontalSlider, SIGNAL(valueChanged(int)),stackedWidget, SLOT(setSpeed(int)));

    connect(actionAdjust_Fiber_Opacity,SIGNAL(triggered(bool)), this, SLOT(AdjustOpacityFiber()));
    connect(actionAdjust_Surface_Opacity,SIGNAL(triggered(bool)), this, SLOT(AdjustOpacitySurface()));

    connect(actionWarpSlides,SIGNAL(triggered(bool)),stackedWidget,SLOT(setWrap(bool)));
    connect(actionVertical, SIGNAL(triggered(bool)),stackedWidget,SLOT(setVerticalMode(bool)));


    //connections;
    connect(action_Quit,SIGNAL(triggered()), qApp, SLOT(quit()));
    connect(actionAdd_Sliding_Window,SIGNAL(triggered(bool)), this, SLOT(AddWindow()));
    connect(actionClose_Current_Sliding_Window, SIGNAL(triggered(bool)), this, SLOT(RemoveCurrentWindow()));
    connect(actionAdjust_Speed, SIGNAL(triggered(bool)), this, SLOT(AdjustSpeed()));


    //change type;
    connect(actionEffective_Network,SIGNAL(triggered(bool)),this,SLOT(ChangeConnectivityType(bool)));
    connect(actionFunctional_Network,SIGNAL(triggered(bool)),this,SLOT(ChangeConnectivityType(bool)));
    connect(actionStructrual_Network,SIGNAL(triggered(bool)),this,SLOT(ChangeConnectivityType(bool)));
    connect(actionNone, SIGNAL(triggered(bool)), this, SLOT(ChangeConnectivityType(bool)));

    // for animation;
    connect(actionNextProfile,SIGNAL(triggered(bool)),stackedWidget, SLOT(slideInNext()));
    connect(actionPreviousProfile,SIGNAL(triggered(bool)),stackedWidget, SLOT(slideInPrev()));
    connect(stackedWidget, SIGNAL(animationFinished()), this, SLOT(WindowChanged()));

    // for template;
    connect(actionDock_Template,SIGNAL(triggered(bool)), this, SLOT(DockWindow4Template()));
    //connect(stackedWidget, SIGNAL(animationFinished()), this, SLOT(SyncTemplate2CurrentWnd()));
//     this->statusBar()->showMessage(QString::fromUtf8("Really status"),200000);

    //build atlases;
    string atlasPath = string(getenv("FSL_DIR"))+string("/atlases/");
    BuildAtlases(atlasPath+"HarvardOxford-Cortical.xml",m_atlases["HarvardOxford-Cortical"]);
    BuildAtlases(atlasPath+"Talairach.xml",m_atlases["Talairach"]);
    BuildAtlases(atlasPath+"MNI.xml",m_atlases["MNI"]);

}
GUIV2::~GUIV2()
{

}

void GUIV2::NetworkRemoved(void )
{
    actionAdjust_ROIs->setDisabled(true);
    actionRemove_Network->setDisabled(true);
}

void GUIV2::NetworkAdded(void )
{
    actionAdjust_ROIs->setDisabled(false);
    actionRemove_Network->setDisabled(false);
    actionDock_Template->setDisabled(false);
    this->ChangeWindTitle();
    if (templateWindow)
    {
        this->SyncTemplate2CurrentWnd();
    }
}

int GUIV2::ChangeConnectivityType(bool status)
{
    QAction* networkType=static_cast<QAction*>(QObject::sender()) ;

    if (networkType== actionNone) {
        if (status)
        {
            cout<<"disable all networkTypes"<<endl;
            actionFunctional_Network->setDisabled(true);
            actionEffective_Network->setDisabled(true);
            actionStructrual_Network->setDisabled(true);
            currentWindow->Change2NoneConnectivity();
            modeType=NONEType;
        } else {
            cout<<"enable all networkTypes"<<endl;
            actionFunctional_Network->setDisabled(false);
            actionEffective_Network->setDisabled(false);
            actionStructrual_Network->setDisabled(false);
            if (actionStructrual_Network->isChecked())
                currentWindow->Change2StructuralConnectivity();
            else if (actionFunctional_Network->isChecked())
                currentWindow->Change2FunctionalConnectivity();
            else if (actionEffective_Network->isChecked())
                currentWindow->Change2EffectiveConnectivity();
            else
            {
                cout<<"error, no network selected."<<endl;
            }
        }

    }

    if (networkType==actionStructrual_Network && status) {
        cout<<"StructrualConnectivity"<<endl;
        actionEffective_Network->setChecked(false);
        actionFunctional_Network->setChecked(false);
        actionStructrual_Network->setChecked(true);
        modeType=StructrualConnectivity;
    }

    if (networkType==actionFunctional_Network && status) {
        cout<<"FunctionalConnectivity"<<endl;
        actionEffective_Network->setChecked(false);
        actionFunctional_Network->setChecked(true);
        actionStructrual_Network->setChecked(false);
        modeType=FunctionalConnectivity;

    }
    else if ( networkType == actionEffective_Network && status ) {
        actionFunctional_Network->setChecked(false);
        actionEffective_Network->setChecked(true);
        actionStructrual_Network->setChecked(false);
        modeType=EffectiveConnectivity;
        cout<<"EffectiveConnectivity"<<endl;
    }
}

void GUIV2::AddWindow()
{
//   cout<<"previous count "<<stackedWidget->count()<<endl;

    int newID= stackedWidget->count();
    QString objName("profile_");
    objName+=QString::number(newID);

    cout<<"add new window: "<<newID<<" "<<objName.toStdString()<<endl;
    CUnit* newSlideWidget=new CUnit();
    newSlideWidget->setObjectName(objName);
    stackedWidget->addWidget(newSlideWidget);
    stackedWidget->slideInIdx(stackedWidget->count()-1);
}
void GUIV2::RemoveCurrentWindow()
{
//   cout<<"previous count "<<stackedWidget->count()<<endl;
    cout<<"remove window "<<stackedWidget->currentIndex()<<" "<<stackedWidget->currentWidget()->objectName().toStdString()<<endl;
    this->RemoveConnections4Window(currentWindow);

    if (templateWindow && currentWindow->allNetworkNodes.size())
        RemoveSyncTemplate2CurrentWnd();


    stackedWidget->removeWidget(currentWindow);

    if (0==stackedWidget->count()) //the last one;
    {
        currentWindow= new CUnit();
        this->setWindowTitle(QString::fromStdString(("Visual Analytic Toolkit for Brain Network")));
		stackedWidget->addWidget(currentWindow);
		stackedWidget->slideInIdx(0);
        previousWindow=currentWindow;
        BuildConnections4Window(currentWindow);
    }
    else
    {
	stackedWidget->slideInPrev();
        currentWindow= static_cast<CUnit*> (stackedWidget->currentWidget());
        //build new connections;
        BuildConnections4Window(currentWindow);
        SyncTemplate2CurrentWnd();
    }

}

void GUIV2::WindowChanged()
{
    //remove old connections;
    RemoveConnections4Window(currentWindow);
    previousWindow=currentWindow;
    cout<<"change window to "<<stackedWidget->currentIndex()<<" "<<stackedWidget->currentWidget()->objectName().toStdString()<<endl;
    currentWindow= static_cast<CUnit*> (stackedWidget->currentWidget());

    //build new connections;
    BuildConnections4Window(currentWindow);
    SyncTemplate2CurrentWnd();

    //change window title;

    this->ChangeWindTitle();
    //change network type of current window to the main control;
    if (currentWindow->allNetworkNodes.size())
    {

        switch (modeType) {
        case NONEType :
            currentWindow->Change2NoneConnectivity();
            break;
        case FunctionalConnectivity :
            currentWindow->Change2FunctionalConnectivity() ;
            break;
        case EffectiveConnectivity :
            currentWindow->Change2EffectiveConnectivity();
            break;
        case StructrualConnectivity :
            break;
        default :
            ;
        } //end of switch;
    }//end of if;

// emit currentWindow->ROIAdjCheckBoxChanged(0);
}

void GUIV2::InitializeWindow()
{

    stackedWidget->removeWidget(profile2);
    stackedWidget->setSpeed(200);
    stackedWidget->setCurrentIndex(0);
    this->setWindowTitle(QString::fromStdString(("Visual Analytic Toolkit for Brain Network")));
    currentWindow=profile1;
    previousWindow=profile1;
    BuildConnections4Window(currentWindow);

    //for template window;
    templateWindow= NULL;
}

void GUIV2::AdjustSpeed()
{
    speedAdj->show();
}
void GUIV2::AdjustOpacityFiber()
{
    opacityFiberAdj->show();
}
void GUIV2::AdjustOpacitySurface()
{
    opacitySurfaceAdj->show();
}

void GUIV2::BuildConnections4Window(CUnit* theWnd)
{
    theWnd->pTemplateWnd = templateWindow; 
    connect(actionOpen_Surface,SIGNAL(triggered()),theWnd,SLOT(LoadSurface()));
    connect(actionOpen_Fiber, SIGNAL(triggered()), theWnd,SLOT(LoadFiber()));
    connect(actionOpen_Volume, SIGNAL(triggered()), theWnd, SLOT(LoadVolume()));
    connect(actionOpen_Bolds,SIGNAL(triggered()),theWnd, SLOT(LoadBOLDS()));
    connect(actionOpen_GM_Mask,SIGNAL(triggered()),theWnd,SLOT(LoadGM()));

    connect(actionDisconnect_Camera,SIGNAL(triggered()), theWnd,SLOT(DisconnectCamera()));
    connect(actionLoad_Network, SIGNAL(triggered(bool)), theWnd, SLOT(LoadNetwork()));
    connect(actionAdjust_ROIs, SIGNAL(triggered(bool)), theWnd, SLOT(ShowROIAdjs()));
    connect(actionFunctional_Network, SIGNAL(triggered(bool)),theWnd, SLOT(Change2FunctionalConnectivity()));
    connect(actionEffective_Network,SIGNAL(triggered(bool)),theWnd,SLOT(Change2EffectiveConnectivity()));
    connect(actionDefine_Now, SIGNAL(triggered(bool)), theWnd, SLOT(ShowROIDefDial()));


    connect(actionLoad_Profile, SIGNAL(triggered(bool)), theWnd, SLOT(LoadProfiles()));
    connect(actionSave_Network, SIGNAL(triggered(bool)), theWnd, SLOT(SaveNetwork()));
    connect(actionGenerate_Network, SIGNAL(triggered(bool)),theWnd, SLOT(GenerateNetworkViaDefinedROIS()));
    connect(theWnd, SIGNAL(NetworkAdded()),this,SLOT(NetworkAdded()));
    connect(theWnd,SIGNAL(NetworkRemoved()),this,SLOT(NetworkRemoved()));
    connect(actionNone,SIGNAL(triggered(bool)), theWnd, SLOT(Change2NoneConnectivity()));
    connect(actionStructrual_Network, SIGNAL(triggered(bool)), theWnd, SLOT(Change2StructuralConnectivity()));
    connect(actionRemove_Network,SIGNAL(triggered(bool)),theWnd, SLOT(RemoveNetwork()));
    connect(opacityFiberAdj, SIGNAL(ValueChanged(float)),theWnd, SLOT(ChangeOpacityFiber(float)));
    connect(opacitySurfaceAdj,  SIGNAL(ValueChanged(float)),theWnd, SLOT(ChangeOpacitySurface(float)));
    connect(actionSave_Profile_As, SIGNAL(triggered(bool)), theWnd,SLOT(SaveProfileAs()));
    connect(actionSave_Profile, SIGNAL(triggered(bool)),theWnd, SLOT(SaveProfile()));

//     if(currentWindow->allNetworkNodes.size() && templateWindow )
//       this->SyncTemplate2CurrentWnd();
    connect(theWnd,SIGNAL(ROIAdjCheckBoxChanged(int)),this,SLOT(SyncTemplate2CurrentWnd(int)));


    connect(actionSwitch2T1,SIGNAL(triggered(bool)),theWnd, SLOT(Switch2T1Vol(bool)));
    connect(theWnd,SIGNAL(UpdateStatusBar(string,int)),this,SLOT(UpdateStatusInfo(string,int)));
//     connect(theWnd, SIGNAL(ROIViewChanged(int)), this, SLOT(SyncTemplate2CurrentWnd(int)));
}

void GUIV2::RemoveConnections4Window(CUnit* theWnd)
{

    disconnect(actionOpen_Surface,SIGNAL(triggered()),theWnd,SLOT(LoadSurface()));
    disconnect(actionOpen_Fiber, SIGNAL(triggered()), theWnd,SLOT(LoadFiber()));
    disconnect(actionOpen_Volume, SIGNAL(triggered()), theWnd, SLOT(LoadVolume()));
    disconnect(actionOpen_Bolds,SIGNAL(triggered()),theWnd, SLOT(LoadBOLDS()));
    disconnect(actionDisconnect_Camera,SIGNAL(triggered()), theWnd,SLOT(DisconnectCamera()));
    disconnect(actionLoad_Network, SIGNAL(triggered(bool)), theWnd, SLOT(LoadNetwork()));
    disconnect(actionAdjust_ROIs, SIGNAL(triggered(bool)), theWnd, SLOT(ShowROIAdjs()));
    disconnect(actionFunctional_Network, SIGNAL(triggered(bool)),theWnd, SLOT(Change2FunctionalConnectivity()));
    disconnect(actionEffective_Network,SIGNAL(triggered(bool)),theWnd,SLOT(Change2EffectiveConnectivity()));
    disconnect(actionDefine_Now, SIGNAL(triggered(bool)), theWnd, SLOT(ShowROIDefDial()));
    disconnect(actionLoad_Profile, SIGNAL(triggered(bool)), theWnd, SLOT(LoadProfiles()));
    disconnect(actionSave_Network, SIGNAL(triggered(bool)), theWnd, SLOT(SaveNetwork()));
    disconnect(actionGenerate_Network, SIGNAL(triggered(bool)),theWnd, SLOT(GenerateNetworkViaDefinedROIS()));
    disconnect(theWnd, SIGNAL(NetworkAdded()),this,SLOT(NetworkAdded()));
    disconnect(theWnd,SIGNAL(NetworkRemoved()),this,SLOT(NetworkRemoved()));
    disconnect(opacityFiberAdj, SIGNAL(ValueChanged(float)),theWnd, SLOT(ChangeOpacityFiber(float)));
    disconnect(opacitySurfaceAdj,  SIGNAL(ValueChanged(float)),theWnd, SLOT(ChangeOpacitySurface(float)));
    disconnect(actionNone,SIGNAL(triggered(bool)), theWnd, SLOT(Change2NoneConnectivity()));
    disconnect(actionStructrual_Network, SIGNAL(triggered(bool)), theWnd, SLOT(Change2StructuralConnectivity()));
    disconnect(actionRemove_Network,SIGNAL(triggered(bool)),theWnd, SLOT(RemoveNetwork()));
    disconnect(actionSave_Profile_As, SIGNAL(triggered(bool)), theWnd,SLOT(SaveProfileAs()));
    disconnect(actionSave_Profile, SIGNAL(triggered(bool)),theWnd, SLOT(SaveProfile()));
    disconnect(actionOpen_GM_Mask,SIGNAL(triggered(bool)),theWnd,SLOT(LoadGM()));
    disconnect(actionSwitch2T1,SIGNAL(triggered(bool)),theWnd, SLOT(Switch2T1Vol(bool)));
    disconnect(theWnd,SIGNAL(UpdateStatusBar(string,int)),this,SLOT(UpdateStatusInfo(string,int)));
    disconnect(theWnd,SIGNAL(ROIAdjCheckBoxChanged(int)),this,SLOT(SyncTemplate2CurrentWnd(int)));

}


void GUIV2::DockWindow4Template()
{
    if (NULL==templateWindow)
    {
        templateWindow= new CUnit();
        templateWindow->pTemplateWnd = templateWindow;
        templateWindow->setWindowTitle(QString::fromUtf8("Template"));
        templateWindow->LoadProfiles();
        if (templateWindow->prfName.size())
        {
            templateWindow->GenerateNetworkControllorAdjs( templateWindow->allNetworkNodesActors);

            ///connect some of the main controls with template;
            connect(actionDisconnect_Camera,SIGNAL(triggered()), templateWindow,SLOT(DisconnectCamera()));
            connect(actionFunctional_Network, SIGNAL(triggered(bool)),templateWindow, SLOT(Change2FunctionalConnectivity()));
            connect(actionEffective_Network,SIGNAL(triggered(bool)),templateWindow,SLOT(Change2EffectiveConnectivity()));
            connect(actionStructrual_Network, SIGNAL(triggered(bool)), templateWindow, SLOT(Change2StructuralConnectivity()));
            connect(opacityFiberAdj, SIGNAL(ValueChanged(float)),templateWindow, SLOT(ChangeOpacityFiber(float)));
            connect(opacitySurfaceAdj,  SIGNAL(ValueChanged(float)),templateWindow, SLOT(ChangeOpacitySurface(float)));
            connect(actionNone,SIGNAL(triggered(bool)),templateWindow, SLOT(Change2NoneConnectivity()));
            connect(actionSwitch2T1,SIGNAL(triggered(bool)),templateWindow,SLOT(Switch2T1Vol(bool)));


            ///build roi_adj connections with the template;
            this->SyncTemplate2CurrentWnd();

        }
        else
        {
            delete templateWindow;
            templateWindow=NULL;
            this->DockWindow4Template();
        }
    }

    templateWindow->show();

}
void GUIV2::SyncTemplate2CurrentWnd()
{

    if (templateWindow && currentWindow->allNetworkNodes.size())
    {

        ///sync the inner states between template and current windows;

        templateWindow->roiCheckStatus=currentWindow->roiCheckStatus;
        templateWindow->fiberCheckStatus=currentWindow->fiberCheckStatus;
        templateWindow->boldCheckStatus=currentWindow->boldCheckStatus;

        int index=0;
        templateWindow->ChangeROIFIberVisiability(*templateWindow->allroiAdjs[index]);
        templateWindow->ChangeROIBOLDSVisiability(*templateWindow->allroiAdjs[index]);
        templateWindow->ChangeROISphereVisiability(*templateWindow->allroiAdjs[index]);
	if(currentWindow->roiPicked>=0)
	  templateWindow->UpdateVolViewSlices(currentWindow->roiPicked);


	currentWindow->SetTemplate(templateWindow);



        templateWindow->fiberWidget->update();
        templateWindow->surfWidget->update();
        templateWindow->BoldWidget->update();

    }

}

void GUIV2::SyncTemplate2CurrentWnd(int index)
{
    if (templateWindow)
    {

        ///sync the inner states between template and current windows;

        templateWindow->roiCheckStatus=currentWindow->roiCheckStatus;
        templateWindow->fiberCheckStatus=currentWindow->fiberCheckStatus;
        templateWindow->boldCheckStatus=currentWindow->boldCheckStatus;

//   if(templateWindow->fiberCheckStatus[index])
        templateWindow->ChangeROIFIberVisiability(*templateWindow->allroiAdjs[index]);
//   if(templateWindow->boldCheckStatus[index])
        templateWindow->ChangeROIBOLDSVisiability(*templateWindow->allroiAdjs[index]);
//   if(templateWindow->roiCheckStatus[index])
        templateWindow->ChangeROISphereVisiability(*templateWindow->allroiAdjs[index]);

	templateWindow->UpdateVolViewSlices(index);
		
        templateWindow->fiberWidget->update();
        templateWindow->surfWidget->update();
        templateWindow->BoldWidget->update();

    }

}

void GUIV2::RemoveSyncTemplate2CurrentWnd()
{
//     for(int index=0; index< templateWindow->allNetworkNodes.size(); ++index)
//     {
//       disconnect( currentWindow->allroiAdjs[index], SIGNAL(RoiCheckBoxShowChange(int)), templateWindow, SLOT(ChangeROISphereVisiability(int)));
//       disconnect( currentWindow->allroiAdjs[index], SIGNAL(RoiCheckBoxFiberChange(int)),templateWindow, SLOT(ChangeROIFIberVisiability(int)));
//       disconnect( currentWindow->allroiAdjs[index], SIGNAL(RoiCheckBoxBOLDChange(int)),templateWindow, SLOT(ChangeROIBOLDSVisiability(int)));
//     }
}

void GUIV2::ChangeWindTitle(void )
{
    string currentObjName= currentWindow->prfName;
    if (currentObjName.size())
        this->setWindowTitle(QString::fromStdString(string("Visual Analytic Toolkit for Brain Network : ")+currentObjName));
    else
        this->setWindowTitle(QString::fromStdString(("Visual Analytic Toolkit for Brain Network")));
}

int GUIV2::BuildAtlases(string xmlFile, std::map< int, string >& atlasMap)
{

    fstream fstrm;
    fstrm.open(xmlFile.c_str(),ios::in);
    if (NULL==fstrm)
    {
        cerr<<"can not open file: "<<xmlFile<<endl;
        return 1;
    }
    string tmpString;

    while (!fstrm.eof())
    {
        fstrm>>tmpString;
        if (tmpString=="<data>")
            break;
    }

    while (!fstrm.eof())
    {
        getline(fstrm,tmpString);
        if (tmpString.size()==0)
            continue;

        if (tmpString.size()<20)
            break;

        int index_s= tmpString.find_first_of("\"");
        int index_e= tmpString.find_first_of("\"",index_s+1);
        index_s+=1;
        string index_str= tmpString.substr(index_s,index_e-index_s);
        int index_int= atoi(index_str.c_str());

        int label_s= tmpString.find_first_of(">");
        label_s+=1;
        int label_e= tmpString.size()-8;
        string theLabel = tmpString.substr(label_s,label_e-label_s);

        atlasMap[index_int]=theLabel;
    }

    return 0;
}


void GUIV2::UpdateStatusInfo(string atlasName, int label)
{
    // check whether or not label =0;
    string contents=atlasName+" :: ";

    if (label)
    {
        contents+= m_atlases[atlasName][label-1];
    }
    else
    {
        contents+="N/A";
    }

    this->statusbar->showMessage(QString::fromStdString(contents));
    cout<<contents<<endl;
}
