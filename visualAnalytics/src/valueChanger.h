/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#ifndef __SPEED_ADJ__H
#define __SPEED_ADJ__H
#include "ui_valueChanger.h"
#include <QWidget>

class CQValueChanger: public QWidget, public Ui::QValueChanger
{
  Q_OBJECT
  public:
    CQValueChanger(float min,float max);
    virtual ~CQValueChanger();
    float GetValue(void);
    void SetMax(float maxValue);
    void SetMin(float minValue);
    
private:
  float maxV;
  float minV; 
  
public slots:
  void ChangeValue(void); 
signals:
  void ValueChanged(float);
}; 

#endif
