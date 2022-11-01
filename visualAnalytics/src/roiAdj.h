/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#ifndef __ROI_ADJ__H
#define __ROI_ADJ__H

#include "ui_roiAdj.h"
#include <QMainWindow>
#include <QWidget>
#include <string>
using std::string;



class myROIAdj : public QWidget, public Ui::adjROIW
{
    Q_OBJECT
public:
    myROIAdj();
    ~myROIAdj();
    void SetROIName(string name);
    const string& GetROIName(void) const;
    int GetROIID(void);
    void SetROIID(int);
    const float& Getx(void);
    const float& Gety(void);
    const float& Getz(void);
    const float& GetRadius(void);
    void SetCoordinates(float x, float y, float z);
public slots:
    void ChangeDispFromCoords( );
    void UpdateCoordX(void);
    void UpdateCoordY(void);
    void UpdateCoordZ(void);
    void ChangeCoordsFromDisp(void);
    void CheckBoxShowChanged(void);
    void CheckBoxFiberChanged(void);
    void CheckBoxBOLDChanged(void);
    void ChangeRadius(void);

signals:
    void RoiCoordChanged(int );
    void RoiRadiusChanged(int );
    void RoiCheckBoxShowChange(int );
    void RoiCheckBoxFiberChange(int);
    void RoiCheckBoxBOLDChange(int);

private:
    float cx;
    float cy;
    float cz;
    string roiName;
    int roiID;
    float radius;

};



#endif