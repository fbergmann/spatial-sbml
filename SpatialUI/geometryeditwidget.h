#ifndef GEOMETRY_EDIT_WIDGET_H
#define GEOMETRY_EDIT_WIDGET_H


#include <QDialog>
#include <QString>

class QAbstractButton;

namespace Ui {
    class GeometryEditWidget;
}

class SpatialSimulator;

class GeometryEditWidget : public QDialog 
{
  Q_OBJECT

public:

  GeometryEditWidget( SpatialSimulator* simulator = NULL,  QWidget * parent = NULL, Qt::WindowFlags f = 0 );

  SpatialSimulator* getSimulator() { return mpSimulator; }
  
  void setSimulator(SpatialSimulator* simulator);
  void setLastDir(const QString& lastDir) { msLastDir = lastDir;}
  public slots:

    void handleButtons(QAbstractButton * button);
    void	compartmentChanged ( int index);

private:

  void updateUI();
  
  SpatialSimulator* mpSimulator;
  QString msLastDir;
    Ui::GeometryEditWidget *ui;

};


#endif // GEOMETRY_EDIT_WIDGET_H