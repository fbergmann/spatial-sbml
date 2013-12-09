#include "geometryeditwidget.h"
#include "ui_geometryeditwidget.h"

#include "spatialsimulator.h"

#include <sbml/Model.h>
#include <sbml/Compartment.h>
#include <sbml/packages/spatial/common/SpatialExtensionTypes.h>
#include <sbml/packages/spatial/extension/SpatialModelPlugin.h>

#include "spatialstructs.h"

#include <QDialog>
#include <QImage>
#include <QPicture>
#include <QLabel>
#include <QAbstractButton>
#include <QPushButton>
#include <QDialogButtonBox>
#include <QFileDialog>

GeometryEditWidget::GeometryEditWidget( SpatialSimulator *simulator, QWidget * parent , Qt::WindowFlags f )
  : QDialog(parent, f), mpSimulator(simulator), msLastDir(""), mNeedReload(false)
{
  ui = new Ui::GeometryEditWidget();
  ui->setupUi(this);

  connect(ui->lstCompartments, SIGNAL(currentRowChanged (int)), this, SLOT(compartmentChanged(int)));
  //connect(ui->boxOkCancel, SIGNAL(clicked(QAbstractButton*)), this, SLOT(handleButtons(QAbstractButton*)));
  connect(ui->boxOpenSave, SIGNAL(clicked(QAbstractButton*)), this, SLOT(handleButtons(QAbstractButton*)));
  connect(ui->cmdFlipReorder, SIGNAL(clicked()), this, SLOT(flipVolumeOrder()));
  updateUI();
}

bool GeometryEditWidget::needReload() const
{
  return mNeedReload;
}

void GeometryEditWidget::setSimulator(SpatialSimulator *simulator)
{
  mpSimulator = simulator;
  updateUI();
}

void	GeometryEditWidget::compartmentChanged ( int compIndex)
{
  if (mpSimulator == NULL) 
    return;

  int maxX = mpSimulator->getXDim();
  int maxY = mpSimulator->getYDim();
  QImage image(maxX, maxY, QImage::Format_RGB32);

  int length;
  int* geometry = mpSimulator->getGeometry(mpSimulator->getModel()->getCompartment(compIndex)->getId(), length);
  int* values = mpSimulator->getBoundary(mpSimulator->getModel()->getCompartment(compIndex)->getId(), length);
//  boundaryType* bounds = mpSimulator->getBoundaryType(mpSimulator->getModel()->getCompartment(compIndex)->getId(), length);
  double *xVal = mpSimulator->getX(length);
  double *yVal = mpSimulator->getY(length);

  double flipHorizontally = ui->chkFlipHorizontally->isChecked() ? maxX - 1 : 0;

  double flipVertically = ui->chkFlipVertically->isChecked() ? maxY - 1 : 0;
  double signH = flipHorizontally  != 0 ? -1 : 1;
  double signV = flipVertically  != 0 ? -1 : 1;


  int Xdiv = maxX; int Ydiv = maxY;
  int Xindex = 2 * Xdiv - 1, Yindex = 2 * Ydiv - 1;
  int X, Y;
  int index;
  int stride = 2;
  for (Y = 0; Y < Yindex; Y += stride ) 
  {
    for (X = 0; X < Xindex; X += stride ) 
    {
      index = Y * Xindex + X;
      int x =  xVal[index];
      int y =  yVal[index];
      if (geometry[index] != 0)
      {
        image.setPixel(flipHorizontally + signH*x, 
          flipVertically + signV*y ,
          qRgb(255,255,255));
      }
      else
      {
        image.setPixel(flipHorizontally + signH*x, 
          flipVertically + signV*y ,
          qRgb(0,0,0));

      }
      if (values[index] != 0)
      {
        image.setPixel(flipHorizontally + signH*x, 
          flipVertically + signV*y ,
          qRgb(255,0,0));
      }

      //switch (bounds[index])
      //{
      //case Xp:
      //   /*image.setPixel(flipHorizontally + signH*x, 
      //    flipVertically + signV*y ,
      //    QColor(Qt::yellow).rgb());*/
      //  break;
      //case Xm:
      //   image.setPixel(flipHorizontally + signH*x, 
      //    flipVertically + signV*y ,
      //    QColor(Qt::green).rgb());
      //  break;
      //case Yp:
      //   image.setPixel(flipHorizontally + signH*x, 
      //    flipVertically + signV*y ,
      //    QColor(Qt::magenta).rgb());
      //  break;
      //case Ym:
      //   image.setPixel(flipHorizontally + signH*x, 
      //    flipVertically + signV*y ,
      //    QColor(Qt::blue).rgb());
      //  break;

      //default:
      //  /*image.setPixel(flipHorizontally + signH*x, 
      //    flipVertically + signV*y ,
      //    QColor(Qt::yellow).rgb());*/
      //  break;
      //}

    }
  }

  ui->imageLabel->setPixmap(QPixmap::fromImage(image));
  ui->imageLabel->repaint();
}

bool isRedish (const QColor& color)
{

  int red = color.red();
  int green = color.green();
  int blue = color.blue();
  return  red > 200 && green < 50 && blue < 50;
}

void GeometryEditWidget::handleButtons(QAbstractButton * button)
{
  if (mpSimulator == NULL) 
    return;

  QDialogButtonBox::StandardButton standard = ui->boxOpenSave->standardButton ( button );

  //if (standard == QDialogButtonBox::NoButton)
  //  standard = ui->boxOpenSave->standardButton ( button );

  if (standard == QDialogButtonBox::Open)  
  {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Geometry file"), msLastDir, tr("Image files (*.png *.jpg *.jpeg)"));
    if (fileName.isEmpty())
      return;
    QFileInfo info(fileName);
    msLastDir = info.absoluteDir().absolutePath();
    QImage temp; temp.load(fileName);


    int maxX = mpSimulator->getXDim();
    int maxY = mpSimulator->getYDim();

    int length;
    int* values = mpSimulator->getBoundary(mpSimulator->getModel()->getCompartment(ui->lstCompartments->currentRow())->getId(), length);
    int* geometry = mpSimulator->getGeometry(mpSimulator->getModel()->getCompartment(ui->lstCompartments->currentRow())->getId(), length);
    boundaryType* bounds = mpSimulator->getBoundaryType(mpSimulator->getModel()->getCompartment(ui->lstCompartments->currentRow())->getId(), length);
    double *xVal = mpSimulator->getX(length);
    double *yVal = mpSimulator->getY(length);

    double flipHorizontally = ui->chkFlipHorizontally->isChecked() ? maxX - 1 : 0;

    double flipVertically = ui->chkFlipVertically->isChecked() ? maxY - 1 : 0;
    double signH = flipHorizontally  != 0 ? -1 : 1;
    double signV = flipVertically  != 0 ? -1 : 1;


    int Xdiv = maxX; int Ydiv = maxY;
    int Xindex = 2 * Xdiv - 1, Yindex = 2 * Ydiv - 1;
    int X, Y;
    int index;
    int stride = 2;
    for (Y = 0; Y < Yindex; Y += stride ) 
    {
      for (X = 0; X < Xindex; X += stride ) 
      {
        index = Y * Xindex + X;
        int x =  xVal[index];
        int y =  yVal[index];

        QColor pixel = QColor(temp.pixel(flipHorizontally + signH*x, 
          flipVertically + signV*y));


        if (isRedish(pixel))
        {
          values[index] = 1;          
          if (x - 1 >= 0  && 
            QColor(temp.pixel(flipHorizontally + signH*(x-1), 
          flipVertically + signV*y)) == QColor(Qt::black))
          {
            bounds[index].isBofXm = true;
            bounds[index].isBofXp = false;
            bounds[index].isBofYm = false;
            bounds[index].isBofYp = false;
            bounds[index].isBofZm = false;
            bounds[index].isBofZp = false;
          }          
          else if  (y - 1 >= 0  && 
            QColor(temp.pixel(flipHorizontally + signH*(x), 
          flipVertically + signV*(y-1))) == QColor(Qt::black))
          {
            bounds[index].isBofXm = false;
            bounds[index].isBofXp = false;
            bounds[index].isBofYm = true;
            bounds[index].isBofYp = false;
            bounds[index].isBofZm = false;
            bounds[index].isBofZp = false;
          }
          else if  (y + 1 < maxY  && 
            QColor(temp.pixel(flipHorizontally + signH*(x), 
            flipVertically + signV*(y+1)))== QColor(Qt::black))
          {
            bounds[index].isBofXm = false;
            bounds[index].isBofXp = false;
            bounds[index].isBofYm = false;
            bounds[index].isBofYp = true;
            bounds[index].isBofZm = false;
            bounds[index].isBofZp = false;
          }
          else if  (x + 1 < maxX  && 
            QColor(temp.pixel(flipHorizontally + signH*(x+1), 
          flipVertically + signV*(y))) == QColor(Qt::black))
          {
            bounds[index].isBofXm = false;
            bounds[index].isBofXp = true;
            bounds[index].isBofYm = false;
            bounds[index].isBofYp = false;
            bounds[index].isBofZm = false;
            bounds[index].isBofZp = false;
          }
        }
        else
        {
          values[index] = 0;
          bounds[index].isBofXm = false;
          bounds[index].isBofXp = true;
          bounds[index].isBofYm = false;
          bounds[index].isBofYp = false;
          bounds[index].isBofZm = false;
          bounds[index].isBofZp = false;
          if (pixel != Qt::black)
          {
            geometry[index] = 1;
          }
          else
          {
            geometry[index] = 0;
          }

        }
      }
    }

    compartmentChanged(ui->lstCompartments->currentRow());

    
  }
  else if (standard == QDialogButtonBox::Save)
  {
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Geometry file"), msLastDir, tr("Image files (*.png *.jpg *.jpeg)"));
    
    if (fileName.isEmpty())
      return;

    QFileInfo info(fileName);
    if (info.completeSuffix().isEmpty())
      fileName += ".png";
    msLastDir = info.absoluteDir().absolutePath();
    QPixmap temp(*ui->imageLabel->pixmap());
    if (!temp.save(fileName))
    {
      cerr << "saving: " << fileName.ascii() << " failed!!!!!" << endl;
    }
  }


}


void GeometryEditWidget::flipVolumeOrder()
{
  if (mpSimulator == NULL) 
    return;
  mpSimulator->flipVolumeOrder();
  mNeedReload = true;
}

void GeometryEditWidget::updateUI()
{
  ui->lstCompartments->clear();
  ui->imageLabel->setPixmap(QPixmap(0,0));

  if (mpSimulator == NULL) 
    return;

  for (unsigned int i = 0; i < mpSimulator->getModel()->getNumCompartments(); ++i)
  {
    const Compartment* comp = mpSimulator->getModel()->getCompartment(i);

    if (comp->isSetName())
      ui->lstCompartments->addItem(comp->getName().c_str());
    else 
      ui->lstCompartments->addItem(comp->getId().c_str());

  }

}
