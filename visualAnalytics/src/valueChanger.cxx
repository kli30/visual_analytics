/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#include "valueChanger.h"
#include <QtGui>
#include <QMenu>

CQValueChanger::CQValueChanger(float min, float max)
{
    this->setupUi(this);
    this->setWindowTitle(QString::fromUtf8("Value Changer"));
    minV=min;
    maxV=max;
    connect(horizontalSlider, SIGNAL(valueChanged(int)), this,SLOT(ChangeValue()));
}
CQValueChanger::~CQValueChanger()
{

}
float CQValueChanger::GetValue(void )
{
    return horizontalSlider->value()*(maxV-minV)*0.001;
}
void CQValueChanger::SetMax(float maxValue)
{
    maxV=maxValue;

}
void CQValueChanger::SetMin(float minValue)
{
    minV=minValue;
}

void CQValueChanger::ChangeValue(void)
{
    emit this->ValueChanged(this->GetValue());
}



