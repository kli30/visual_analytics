/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#ifndef __ROIDEFDIAL__H
#define __ROIDEFDIAL__H
#include "ui_defineROI.h"
#include <QMainWindow>

class myRoiDefDdial: public QDialog, public Ui::Dialog
{
  Q_OBJECT
  public:
    myRoiDefDdial(void);
    virtual ~myRoiDefDdial();
}; 

#endif