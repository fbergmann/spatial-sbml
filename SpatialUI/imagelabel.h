#ifndef IMAGE_LABEL_H
#define IMAGE_LABEL_H

#include <QLabel>
#include <QString>

class QMouseEvent;
class QWidget;

class ImageLabel : public QLabel
{
  Q_OBJECT

public: 
  
  ImageLabel ( QWidget * parent = 0, Qt::WindowFlags f = 0 );
  ImageLabel ( const QString & text, QWidget * parent = 0, Qt::WindowFlags f = 0 );    
  virtual ~ImageLabel();

signals:
  void positionChanged(int x, int y);

protected:

  virtual void	mouseMoveEvent ( QMouseEvent * ev );


};


#endif // IMAGE_LABEL_H