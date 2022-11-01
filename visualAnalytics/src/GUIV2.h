/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/
#ifndef _GUI_h
#define _GUI_h

#include <QMainWindow>
#include "ui_GUIV2.h"
#include "valueChanger.h"
#include "traceMap.h"


class CUnit;
class QAction;
class QString;
class myRoiDefDdial;
class QScrollArea;

using namespace std;
using namespace KML;

class GUIV2 : public QMainWindow, public Ui::MainWindow
{
    Q_OBJECT

public:
    GUIV2();
    ~GUIV2();

public slots:

    int ChangeConnectivityType(bool status);
    void NetworkRemoved(void);
    void NetworkAdded(void);
    void ChangeWindTitle(void);

    void InitializeWindow();

    // for animation;
    void AddWindow();
    void RemoveCurrentWindow();
    void WindowChanged();
    void AdjustSpeed();

    //for visualization;
    void AdjustOpacitySurface();
    void AdjustOpacityFiber();

    // dock a certain windows for comparion;
    void DockWindow4Template();
    void SyncTemplate2CurrentWnd();
    void SyncTemplate2CurrentWnd(int index);

    //update status information;
    void UpdateStatusInfo(string atlasName, int label);


private:

    enum ROIAttrStatus { NotDone, Updated, NeedUpdate};
    enum ConnectivityType {NONEType, FunctionalConnectivity, EffectiveConnectivity , StructrualConnectivity};
    ConnectivityType modeType;

    CUnit* currentWindow;
    CUnit* previousWindow;
    CUnit* templateWindow;

    CQValueChanger* speedAdj;
    CQValueChanger* opacityFiberAdj;
    CQValueChanger* opacitySurfaceAdj;


    void BuildConnections4Window(CUnit* theWnd);
    void RemoveConnections4Window(CUnit* theWnd);
    void RemoveSyncTemplate2CurrentWnd();

    map<string,map<int,string> > m_atlases;
    int BuildAtlases(string xmlFile, std::map< int, string >& atlasMap);

};

#endif // _GUI_h

