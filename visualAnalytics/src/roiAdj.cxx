/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#include "roiAdj.h"
#include "ui_roiAdj.h"
#include "kaimingCommon.h"
#include "sstream"

using namespace std;
using namespace KML;


myROIAdj::myROIAdj()
{
    setupUi(this);
    cx=0;
    cy=0;
    cz=0;
    radius=5;

    ChangeDispFromCoords();
    connect(xchange, SIGNAL(editingFinished()),this,SLOT(UpdateCoordX()));
    connect(ychange, SIGNAL(editingFinished()),this,SLOT(UpdateCoordY()));
    connect(zchange, SIGNAL(editingFinished()),this,SLOT(UpdateCoordZ()));
    connect(coordDisp, SIGNAL(textEdited(QString)), this,SLOT(ChangeCoordsFromDisp()));
    connect(checkBox, SIGNAL(stateChanged(int)), this, SLOT(CheckBoxShowChanged()));
    connect(checkBox_BOLD, SIGNAL(stateChanged(int)),this, SLOT(CheckBoxBOLDChanged()));
    connect(checkBox_Fiber, SIGNAL(stateChanged(int)),this, SLOT(CheckBoxFiberChanged()));
    connect(radiusBox, SIGNAL(editingFinished()), this, SLOT(ChangeRadius()));


}
myROIAdj::~myROIAdj()
{

}

void myROIAdj::ChangeDispFromCoords()
{
    //display the coords in lineEdit;
    stringstream theStrStrm;
    theStrStrm.precision(5);
    theStrStrm<<cx<<"   "<<cy<<"   "<<cz;
    string dispContent;
    getline(theStrStrm, dispContent);
    coordDisp->setText(QString::fromStdString(dispContent));

}

void myROIAdj::UpdateCoordX(void )
{
    if ( xchange->value() != 0 )
    {

        cx+= xchange->value();
        this->ChangeDispFromCoords();
        xchange->setValue(0);
        emit RoiCoordChanged(roiID);

    }
}
void myROIAdj::UpdateCoordY(void )
{
    if ( ychange->value() !=0)
    {
        cy+= ychange->value();
        this->ChangeDispFromCoords();
        ychange->setValue(0);
        emit RoiCoordChanged(roiID);
    }

}
void myROIAdj::UpdateCoordZ(void )
{
    if (0 != zchange->value() )
    {
        cz+= zchange->value();
        this->ChangeDispFromCoords();
        zchange->setValue(0);
        emit RoiCoordChanged(roiID);
    }
}

void myROIAdj::ChangeCoordsFromDisp(void )
{

    string theCoord= coordDisp->text().toStdString();
    stringstream theStrStrm;
    theStrStrm<<theCoord;
    theStrStrm>>cx>>cy>>cz;
    emit RoiCoordChanged(roiID);

}
void myROIAdj::SetROIName(string name)
{
    this->roiIndex->setText(name.c_str());
    this->roiName=name;
}
const std::string& myROIAdj::GetROIName(void ) const
{
    return this->roiName;

}

int myROIAdj::GetROIID(void )
{
    return roiID;
}

void myROIAdj::SetROIID(int id)
{
    this->roiID= id;
}
const float& myROIAdj::Getx(void )
{
    return cx;
}

const float& myROIAdj::Gety(void )
{
    return cy;
}
const float& myROIAdj::Getz(void )
{
    return cz;
}


void myROIAdj::SetCoordinates(float x, float y, float z)
{
    cx=x;
    cy=y;
    cz=z;
    this->ChangeDispFromCoords();
}

void myROIAdj::CheckBoxShowChanged(void )
{
    emit RoiCheckBoxShowChange(roiID);
}

void myROIAdj::CheckBoxBOLDChanged(void )
{
    emit RoiCheckBoxBOLDChange(roiID);
}
void myROIAdj::CheckBoxFiberChanged(void )
{
    emit RoiCheckBoxFiberChange(roiID);
}

const float& myROIAdj::GetRadius(void )
{
    return radius;
}

void myROIAdj::ChangeRadius(void )
{
    if ( this->radiusBox->value() != this->radius )
    {
        this->radius= radiusBox->value();
        emit RoiRadiusChanged(roiID);
    }

}//end of namespace;