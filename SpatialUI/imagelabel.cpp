#include "imagelabel.h"

#include <QMouseEvent>
#include <QPixmap>

#include <iostream>

using namespace std;

   ImageLabel::ImageLabel ( QWidget * parent, Qt::WindowFlags f  )
     : QLabel(parent, f)
   {
   }

   ImageLabel::ImageLabel ( const QString & text, QWidget * parent , Qt::WindowFlags f  )
     : QLabel(text, parent, f)
   {
   }


 ImageLabel::~ImageLabel()
 {
 }

void	
ImageLabel::mouseMoveEvent ( QMouseEvent * ev )
{
  QLabel::mouseMoveEvent(ev);
  const QPixmap* pixmap = QLabel::pixmap();
  if (pixmap == NULL) return;

  int imageX = (int)((double)ev->x()/(double)width()*(double)pixmap->width());
  int imageY = (int)((double)ev->y()/(double)height()*(double)pixmap->height());
/*
  cout 
    //<< "x: " << (int)((double)ev->x()) << " "
    //<< "y: " << (int)((double)ev->y()) << " "
    << "x: " << imageX  << " "
    << "y: " << imageY  << endl;
    */
  emit positionChanged(imageX, imageY);

}