/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#include "qapplication.h"
#include "GUIV2.h"
#include "unit.h"
 
int main(int argc, char** argv)
{
  QApplication app(argc, argv);
  GUIV2 widget;
  widget.show();
  return app.exec();
}

