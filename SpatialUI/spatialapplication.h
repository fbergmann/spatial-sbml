#ifndef SPATIAL_APPLICATION_H
#define SPATIAL_APPLICATION_H

#include <QApplication>
#include <QString>

class SpatialMainWindow;

class SpatialApplication: public QApplication
{
  Q_OBJECT

public:
  SpatialApplication(int & argc, char ** argv);

  virtual ~SpatialApplication();

  virtual bool event(QEvent * event);

  void setMainWindow(SpatialMainWindow * mainWindow);

private:
  SpatialMainWindow * mainWindow;

  QString file;

  bool starting;
};

#endif // SPATIAL_APPLICATION_H