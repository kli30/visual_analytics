/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/
#ifndef __UNIT__H
#define __UNIT__H

#include "ui_unit.h"
#include <QMainWindow>
#include <QWidget>
#include <string>
#include "newimageall.h"
#include "jointModel.h"
#include "vtkImplicitFunction.h"
#include "vtkSphere.h"
#include <vtkExtractPolyDataGeometry.h>
#include "newimageall.h"
#include "triSurface.h"
#include "fibers.h"
#include "vtkSmartPointer.h"
#include <vtkActor.h>
#include "PCACCAAnalysis.h"
#include <vtkTable.h>
#include "ColorSchemeRGB.h"
#include <map>
#include "roiAdj.h"
#include "roiDefDial.h"
#include "vtkTextSource.h"
#include "traceMap.h"
#include "kmPCA.h"


using std::string;
using namespace NEWIMAGE;
using namespace KML;
using namespace Ui;

class vtkSphereSource;
class vtkRenderer;
class vtkEventQtSlotConnect;
class vtkObject;
class vtkCommand;
class QAction;
class QString;
class vtkImageViewer2;
class QLineEdit;
class vtkContextView;
class vtkChartXY;
class vtkPlot;
class vtkPolyData;
class vtkConeSource;
class vtkCylinderSource;
class QScrollArea;
class vtkIdList;
class vtkPropPicker;
class vtkPolyDataMapper;
class vtkTextActor;


class CUnit : public QWidget, public Ui::Form
{
    Q_OBJECT
    friend class GUIV2;
public:

    void SetTemplate(CUnit* ptrTemplate);
    enum ROIAttrStatus { NotDone, Updated, NeedUpdate};
    enum ConnectivityType {NONEType, FunctionalConnectivity, EffectiveConnectivity , StructrualConnectivity};
    CUnit();
    ~CUnit();

public slots:

    void SaveProfileAs(void);
    void SaveProfile(void);
    void AdjustROIS(void);
    void LoadSurface(void);
    void LoadFiber(void);
    void LoadVolume(void);
    void LoadGM(void);
    void LoadBOLDS(void);
    void LoadProfiles(void);
    void LoadNetwork(void);
    void Popup(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command);
    void ShowFullScreen(QAction* obj);

    void UpdateSliceX(const QString& newSliceNum);
    void UpdateSliceY(const QString& newSliceNum);
    void UpdateSliceZ(const QString& newSliceNum);

    void SetSliderValueX(const QString& newText);
    void SetSliderValueY(const QString& newText);
    void SetSliderValueZ(const QString& newText);

    void SetLineEditX(int val);
    void SetLineEditY(int val);
    void SetLineEditZ(int val);

    void VisualizeBoldSingalROI( const std::vector< float >& data);
    void ShowROIDefDial(void);
    int CreatandVisualizeSingleROI(void);
    void LinkCameraofSurfAndFiber(QAction*);
    void DisconnectCamera(void);
    void GenerateNetworkControllorAdjs(vector<vtkSmartPointer<vtkActor> >& allROIActors);
    void ChangeROICoords(myROIAdj& theROI);
    void ChangeROICoords(int roiID);
    void ChangeROISphereVisiability(myROIAdj& theROI);
    void ChangeROIFIberVisiability(myROIAdj& theROI);
    void ChangeROIBOLDSVisiability(myROIAdj& theROI);

    //exclusive for GUI interaction;
    void ChangeROISphereVisiability(int roiID);
    void ChangeROIFIberVisiability(int roiID);
    void ChangeROIBOLDSVisiability(int roiID);
    //exclusive for GUI interaction;

    void ChangeROIRadius( myROIAdj& theROI);
    void ChangeROIRadius( int roiID);
    void MyMouseMove(void);
    void ReleaseMouseMove(void);
    void ShowROIAdjs(void);
    void PopupSurfaceWind(QAction* obj);
    void PopupFiberWind(QAction* obj);
    void GenerateBOLDCurveofROI(int roiID);

    //about network;

    int  Change2FunctionalConnectivity(void);
    int  Change2EffectiveConnectivity(void);
    int  Change2StructuralConnectivity(void);
    int  Change2NoneConnectivity(void);
    void SaveNetwork(void);
    void RemoveNetwork(void);
    void GenerateNetworkViaDefinedROIS(void);

    //about the visualization;
    void ChangeOpacitySurface(float opacity);
    void ChangeOpacityFiber(float opacity);

    void UpdateCheckBoxStatus(myROIAdj* sender);

    //switch to t1 volumes;
    void Switch2T1Vol(bool checked);
    void Switch2DTIVol(void);


signals:
    void NetworkAdded(void);
    void NetworkRemoved(void);
    void ProfileLoaded(void);
    void ROIAdjCheckBoxChanged(int index);
    void ROIFiberChanged(int index);
    void UpdateStatusBar(string atlasName,int label);

private:

    //for tracemap Distance;
    CUnit* pTemplateWnd;
    CTraceMap traceMap;
    CPCA pcaMdl;


    QWidget* m_parient;
    vtkRenderer* surfRenderer;
    vtkRenderer* fiberRenderer;
    vtkEventQtSlotConnect* linkCameraConnection;
    vtkEventQtSlotConnect* mouseMoveConnection;
    vtkEventQtSlotConnect* releaseMouseMoveConnection;
    QAction * actionAdjust_ROIs;

    int roiPicked;
    bool isSingleROIMode;
    bool isBOLDsEnabled;
    bool isVolT1;
    vtkSmartPointer<vtkActor> roiPickedActor;

    vtkImageViewer2* volViewX;
    vtkImageViewer2* volViewY;
    vtkImageViewer2* volViewZ;
    vtkContextView* boldView;
    vtkPlot* linePlot;
    vtkSmartPointer<vtkPropPicker> pPicker;


    ConnectivityType modeType;

    volume4D<float> myBolds;
    volume<float> myDtiGM;
    volume<short> atlasImg;
    myRoiDefDdial* roiDefDial;
    double b0Offset[3];
    double b0Spacing[3];
    int b0Dimensions[3];
    double t1Offset[3];
    double t1Spacing[3];
    int t1Dimensions[3];

    std::string prfName;
    CFibers* myFiber;
    CTriSurface* mySurface;
    CJointModel* myModel;
    CColorSchemeRGB myColorScheme;

    vector<bool> fiberCheckStatus;
    vector<bool> roiCheckStatus;
    vector<bool> boldCheckStatus;

    vtkSmartPointer<vtkTextActor> textActor;
    vtkSmartPointer<vtkTextSource> textSource;
    vtkSmartPointer<vtkPolyDataMapper> textMapper;

    vtkSmartPointer<vtkActor> preROIActor;
    vtkSmartPointer<vtkActor> currROIActor;
    vtkSmartPointer<vtkActor> surfActor;
    vtkSmartPointer<vtkPolyData> surfData;
    vtkSmartPointer<vtkActor> fiberActor;
    vtkSmartPointer<vtkPolyData> fiberData;
    vtkSmartPointer<vtkExtractPolyDataGeometry> extractor;
    vtkSmartPointer<vtkSphere> extractorSphere;

    vtkTable* boldSignalTable;
    vtkChartXY* boldSignalChart;

    //network
    std::vector<Vector3D<float> > allNetworkNodes;
    Matrix edges;
    Matrix functionalEdges;
    Matrix effectiveEdges;
    vector<vtkSmartPointer<vtkActor> > allNetworkNodesActors;
    vector<vtkSmartPointer<vtkSphereSource> > allSphereSources;
    vector< vector<vtkSmartPointer<vtkActor> >  >allNodesLinksActors;
    vector< vector<float> > allNodesPCA1st;
    vector<string> allROINames;

    //need to keep record the cone and rod;
    vector<ROIAttrStatus> allROIAttributeStatus;
    vector<vector<vtkSmartPointer<vtkIdList> > > allROIFIbersVTK;
    vector<vector<int> > allROIFiberIDs;
    vector<myROIAdj*> allroiAdjs;

    QScrollArea* roiListArea;

    //profile
    map<std::string, std::string> prfMap;
    string surf_KEY;
    string dtiB0_KEY;
    string bolds_KEY;
    string fiber_KEY;
    string network_KEY;
    string dtiseg_KEY;
    string t1_KEY;
    string t1Label_KEY;
    string atlas_KEY;

    //private methods;

    void VisualizeCurrentROIOnVolume(const KML::Vector3D< float >& centerCoords, float radius);
    void VisualizeFibersCurrentROI(int roiID);
    void UpdateBOLDSignal4Node(int roiID, float radius);
    void UpdateBoldsSingleROI(const KML::Vector3D< float >& centerCoords, float radius);

    void DoPCA(const std::vector< KML::Vector3D< int > >& theROI, std::vector< float >& evalue, vector< std::vector< float > >& evector);
    int DoPCA( const std::set< int >& surfID, std::vector< float >& evalue, vector< std::vector< float > >& evector);
    int GetPCA1stCurrentROI(int roiID, float radius, std::vector< float >& firstEigenVector);
    void GetPCA1stCurrentROI(const KML::Vector3D< float >& centerCoords, float radius,vector<float>& firstEig);


    void DoNetRendering(const std::vector<Vector3D<float> >& allROIs, const Matrix& edges );
    int  GetRoiCenterDTIVolCoords(const std::string& fileNanme, Vector3D<float>& theCoords);
    int GetRoiCenterSurfCoords(const std::string& fileNanme, Vector3D< float >& theCoords);
    int GetRoiCenterSurfID(Vector3D< float >& volCoords );

    bool VaryfyInputFile(const QString& fileName);
    bool VaryfyInputFile(const std::string& fileName);

    void VisualizeSurface(const std::string& myFileName);
    void VisualizeFibers(const std::string& myFileName);
    void VisualizeBOLDS(const std::string& myFileName);
    void VisualizeB0Volume(const std::string& myFileName);
    void VisualizeT1Volume(const std::string& myFileName);


    int VisualizeNetwork(void);
    int AddArrow1(int x, int y, const vector< Vector3D< float > >& roiCenters, float opacity , float spheraRadius = 4, float coneHeight = 4, float cylinderRadius = 0.3);
    int AddArrow(int x, int y, const vector< Vector3D< float > >& roiCenters, float opacity , float spheraRadius = 4, float coneHeight = 4, float cylinderRadius = 0.3);
    int AddLine(int x, int y,const vector< Vector3D< float > >& roiCenters, float opacity , float spheraRadius = 4,   float cylinderRadius = 0.15 );
    void RemoveNetworkLinksActors(void);
    void AddNetworkLinksActors(void);
    void AddNetworkNodesActors(void);
    void RemoveNetworkNodesActors(void);

    int ParserProfile(const std::string& myFileName);
    int ParserNetworkSettings(const std::string& myFileName);
    void UpdateNetworkVisualization(myROIAdj& theROI);
    void TrackFibersCurrentROI(myROIAdj& theROI);
    int GetROIID(vtkActor* theActor);
    void UpdateFunctionalNetworkWeights(int roiID);
    void UpdateEffectiveNetworkWeights(int roiID);
    void GenerateNetworkROIPCAs(void);
//   void ChangeConnectivityType(QAction* networkType);
    void GetGCStrength(std::vector< float >& x, std::vector< float >& y, double& x2y, double& y2x);
    void GetPearsonCorrelation(std::vector< float >& x, std::vector< float >& y, double& coorel);
    void UpdateVolViewSlices(int roiID);
    void CalculateNetEdgeGC(void);
    void CalculateNetEdgePC(void);
    void LoadNetwork(string netConfigName);
    void UpdateText(int pointPicked,Vector3D<double> coord);
    float CalculateDistanceHTraceMap(int pointPicked);
    void EnableBOLDs(bool option);

    //initialize the box status ;
    void InitializeCheckBoxesStatusFromNetworkControllor(void);

    //update the atlas information on the status bar;
    void UpdateAtlasInfoAtStatusBar(int roiPicked);

    //read meta information
    void ReadMetaInfo4Dti(void);
    void ReadMetaInfo4T1(void);

};

#endif
