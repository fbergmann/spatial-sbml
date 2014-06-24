
#include "spatialmainwidget.h"
#include "spatialsimulator.h"
#include "concentrationpalette.h"
#include "displayitem.h"
#include "ui_spatialmainwidget.h"
#include "simulationthread.h"
#include "geometryeditwidget.h"

#include <SpatialSBML/datahelper.hh>

#include <QApplication>
#include <QtGui>
#include <QObject>
#include <QDir>
#include <QFile>
#include <QListWidget>
#include <QFileDialog>
#include <QMessageBox>
#include <QMenu>
#include <QMenuBar>
#include <QToolBar>
#include <QStatusBar>


#include <sbml/SBMLTypes.h>
#include <sbml/packages/spatial/extension/SpatialSpeciesRxnPlugin.h>
#include <sbml/packages/spatial/extension/SpatialModelPlugin.h>
#include <sbml/packages/spatial/extension/SpatialParameterPlugin.h>

#include <string>
#include <string.h>

#ifdef USE_SBW_INTEGRATION
#include <stdlib.h>
#endif // USE_SBW_INTEGRATION

#define SPATIAL_ANNOTATION_URL "http://fbergmann.github.io/spatial-sbml/settings"


using namespace std;
LIBSBML_CPP_NAMESPACE_USE;

SpatialMainWindow::SpatialMainWindow() : thread(NULL), updating(false), pickerX(0), pickerY(0), maxX(501), maxY(501)
#ifdef USE_SBW_INTEGRATION
  , mpSBWModule(NULL)
  , mSBWAnalyzerModules()
  , mSBWAnalyzerServices()
  , mSBWActionMap()
  , mpSBWActionGroup(NULL)
  , mpSBWMenu(NULL)
  , mpSBWAction(NULL)
  , mSBWIgnoreShutdownEvent(true)
#endif // USE_SBW_INTEGRATION

{
  QWidget *widget = new QWidget();
  ui = new Ui::SpatialMainWindow();
  ui->setupUi(widget);
  setCentralWidget(widget);

  thread = new SimulationThread();
  thread->mpDisplayItems = &displayItems;
  connect(thread, SIGNAL(visualize(double, const QImage&)), this, SLOT(visualize(double, const QImage&)));

  lstSpecies = ui->lstSpecies;

  connect(ui->lstPalettes, SIGNAL(activated(int)), this, SLOT(palatteIndexChanged(int)));
  connect(ui->cmdRemove, SIGNAL(clicked ()), this, SLOT(removeSpecies()));
  connect(ui->cmdAdd, SIGNAL(clicked ()), this, SLOT(addSpecies()));
  connect(ui->cmdExportConc, SIGNAL(clicked ()), this, SLOT(exportConcentration()));
  connect(ui->cmdImportConc, SIGNAL(clicked ()), this, SLOT(importConcentration()));
  connect(ui->lstAssignments, SIGNAL(currentRowChanged(int)), this, SLOT(selectedSpeciesChanged(int)));
  connect(ui->lblImage, SIGNAL(positionChanged(int, int)), this, SLOT(updatePosition(int,int)));
  connect(ui->tblParameters, SIGNAL(cellChanged(int, int)), this, SLOT(parameterChanged(int,int)));
  connect(ui->chkHideBC, SIGNAL(toggled(bool)), this, SLOT(toggledHideBC(bool)));
  connect(ui->chkHideDiff, SIGNAL(toggled(bool)), this, SLOT(toggledHideBC(bool)));

  createActions();
  createMenus();
  createToolBars();
  createStatusBar();

  readSettings();

#ifdef USE_SBW_INTEGRATION

  sbwConnect();
  sbwRegister();
    
#endif // USE_SBW_INTEGRATION


  QDir dir (qApp->applicationDirPath());
  QFileInfoList foo = dir.entryInfoList(QStringList() << "*.txt");
  for (int i=0;i<foo.size();i++) {
    palettes.push_back(new ConcentrationPalette(qApp->applicationDirPath() +"/"+ foo[i].fileName()));
    ui->lstPalettes->addItem(foo[i].fileName());
  }

  setCurrentFile("");
  setUnifiedTitleAndToolBarOnMac(true);
  setWindowIcon(QIcon(":/images/ICON_Spatial_128x128.png"));

  if (palettes.size() == 0)
  {
    ui->lblImage->setText("No Pallette found, this should not be happening, please reinstall.");
    this->setEnabled(false);
  }

}

void SpatialMainWindow::updatePosition(int x, int y)
{
  pickerX = x;
  pickerY = y;
  ui->txtXPos->setText(QString::number(x));
  ui->txtYPos->setText(QString::number(y));

  int transformedX = x;
  int transformedY = maxY-1 - y;

  if (thread == NULL) return;

  if (  (QApplication::mouseButtons().operator&((int)Qt::LeftButton) == (int) Qt::LeftButton)  && 
    ui->cmdEnableDose->isChecked()
    )
  {          
    thread->applyDose(transformedX, transformedY,  ui->lstDose->currentText(), ui->txtDose->text().toDouble());
  }

  // update info 
  vector<pair<string, double> > currentValues = thread->getConcentrationsAt(transformedX, transformedY);
  if (ui->lstConcentrations->count() != currentValues.size())
  {
    ui->lstConcentrations->clear();
    for (unsigned int i = 0; i < currentValues.size(); ++i)
      ui->lstConcentrations->addItem("");
  }

  for (unsigned int i = 0; i < currentValues.size(); ++i)
  {
    pair<string, double> current = currentValues[i];
    ui->lstConcentrations->item(i)->setText(tr("%1 = %2")
      .arg(current.first.c_str())
      .arg(current.second));
  }

}

void SpatialMainWindow::selectedSpeciesChanged(int row)
{
}

void SpatialMainWindow::importConcentration()
{
  QString fileName = QFileDialog::getOpenFileName(this, tr("Import Concentration"), lastDir, tr("Dump files (*.dmp)"));
  if (fileName.isEmpty())
    return;

  lastDir = QFileInfo(fileName).absoluteDir().absolutePath();

  const auto& currentSpecies = ui->lstDose->currentText();
  //const auto& index = ui->lstDose->currentIndex();

  DataHelper helper(fileName.toStdString());
  if (!helper.isValid()) return;

  for (int i = 0; i < maxX; ++i)
    for(int j = 0; j < maxY; ++j)
    {

      int transformedX = i;
      int transformedY = maxY-1 - j;

      thread->applyDose(transformedX, transformedY, currentSpecies, helper(transformedX,transformedY));
    }

}

void SpatialMainWindow::exportConcentration()
{

  QString fileName = QFileDialog::getSaveFileName(this, tr("Export  Concentration"), lastDir, tr("Dump files (*.dmp)"));
  if (fileName.isEmpty())
    return;

  lastDir = QFileInfo(fileName).absoluteDir().absolutePath();

  const QString& currentSpecies = ui->lstDose->currentText();
  //int index = ui->lstDose->currentIndex();

  DataHelper helper(thread->mpSimulator->getXDim(), thread->mpSimulator->getYDim());
  
  for (int x = 0; x < thread->mpSimulator->getXDim(); ++x)
  {
    for (int y = 0; y < thread->mpSimulator->getYDim(); ++y)
    {
      const std::vector<std::pair<std::string, double> >& all = thread->getConcentrationsAt(x, y);
      std::vector<std::pair<std::string, double> >::const_iterator it = all.cbegin();
      bool wrote = false;
      while( it != all.end())
      {
        if ((*it).first == currentSpecies.toStdString())
        {
          helper(x, y) = it->second;
          wrote = true;
          break;
        }
        ++it;
      }
      if (!wrote)
          helper(x, y) = 0;
    }
  }

  helper.writeToFile(fileName.toStdString());
  
}

void SpatialMainWindow::addSpecies()
{
  QString currentSpecies = ui->lstSpecies->currentText();
  QString currentPalette = ui->lstPalettes->currentText();
  int index = ui->lstPalettes->currentIndex();
  DisplayItem* item = NULL;

  vector<DisplayItem*>::iterator iter;
  for(iter = displayItems.begin(); iter != displayItems.end(); ++iter)
  {
    if ((*iter)->getId() == currentSpecies)
    {
      item = *iter;
      break;
    }  
  }

  if (item != NULL)
  {
    item->setPalette(palettes[index]);
  }
  else
  {
    item = new DisplayItem(
      currentSpecies.toStdString(),
      palettes[index]
    );
    displayItems.push_back(item);

    ui->lstAssignments->addItem(currentSpecies);
  }

  item->getPalette()->setMaxValue(ui->txtMaxValue->text().toDouble());

}

void SpatialMainWindow::removeSpecies()
{
  int current = ui->lstAssignments->currentRow();
  QListWidgetItem* item = ui->lstAssignments->takeItem (current);
  if (item == NULL) return;
  QString selected = item->text();

  bool wasRunning = !thread->mPaused;
  pause();

  vector<DisplayItem*>::reverse_iterator iter;

  for(iter = displayItems.rbegin(); iter != displayItems.rend(); ++iter)
  {
    if ((*iter)->getId() == selected)
    {
      displayItems.erase((iter+1).base());
      break;
    }
  }



  if (wasRunning) play();

}

void SpatialMainWindow::palatteIndexChanged(int index)
{
}

void SpatialMainWindow::closeEvent(QCloseEvent *event)
{
  writeSettings();

  stop();

  if (thread != NULL)
  {
    if (thread->isRunning())
    {
      thread->stop();
    }
    delete thread;
  }

  event->accept();
}

void SpatialMainWindow::toggledHideBC(bool)
{
  fillParameters();
}

void SpatialMainWindow::newFile()
{
  setCurrentFile("");

  if (thread != NULL)
  {
    if (thread->isRunning())
    {
      thread->terminate();
      thread->wait();
    }

    delete thread;
    thread = NULL;
  }

}

void SpatialMainWindow::open()
{
  QString fileName = QFileDialog::getOpenFileName(this, tr("Open SBML file"), lastDir, tr("SBML files (*.xml *.sbml)"));
  if (!fileName.isEmpty())
    loadFile(fileName);
}

bool SpatialMainWindow::save()
{
  //if (curFile.isEmpty()) 
  //{
  return saveAs();
  //} 
  //else 
  //{
  //  return saveFile(curFile);
  //}
}

bool SpatialMainWindow::saveImage()
{
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save Image"), lastDir, tr("Image files (*.png *.jpeg *.jpg)"));
  if (fileName.isEmpty())
    return false;

  return saveImageFile(fileName);
}

bool SpatialMainWindow::saveAs()
{
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save SBML file"), lastDir,  tr("SBML files (*.xml *.sbml)"));
  if (fileName.isEmpty())
    return false;

  return saveFile(fileName);
}

void SpatialMainWindow::about()
{
  QMessageBox::about(this, tr("About Spatial UI"),
    tr(
    "The <b>Spatial UI</b> application opens a file with the SBML spatial extension and allows to "
    "simulate it using Akira's simulator.<br/><br/>"
    "It uses <a href='http://sbml.org/Software/libSBML'>libSBML</a> Version: %1")
    .arg(getLibSBMLDottedVersion())
    );
}

void SpatialMainWindow::createActions()
{
  newAct = new QAction(QIcon(":/images/file-128x128.png"), tr("&New"), this);
  newAct->setShortcuts(QKeySequence::New);
  newAct->setStatusTip(tr("Create a new file"));
  connect(newAct, SIGNAL(triggered()), this, SLOT(newFile()));

  openAct = new QAction(QIcon(":/images/folder-open-128x128.png"), tr("&Open..."), this);
  openAct->setShortcuts(QKeySequence::Open);
  openAct->setStatusTip(tr("Open an existing file"));
  connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

  saveAct = new QAction(QIcon(":/images/floppy-128x128.png"), tr("&Save"), this);
  saveAct->setShortcuts(QKeySequence::Save);
  saveAct->setStatusTip(tr("Save the document to disk"));
  connect(saveAct, SIGNAL(triggered()), this, SLOT(save()));

  saveImageAct = new QAction(QIcon(":/images/floppy-128x128.png"), tr("&Save Image"), this);
  saveImageAct->setStatusTip(tr("Save the current Image to disk"));
  connect(saveImageAct, SIGNAL(triggered()), this, SLOT(saveImage()));


  editGeometryAct = new QAction(tr("Edit &Geometry"), this);
  editGeometryAct ->setStatusTip(tr("View / Edit geometries"));
  connect(editGeometryAct, SIGNAL(triggered()), this, SLOT(editGeometry()));

  saveAsAct = new QAction(tr("Save &As..."), this);
  saveAsAct->setShortcuts(QKeySequence::SaveAs);
  saveAsAct->setStatusTip(tr("Save the document under a new name"));
  connect(saveAsAct, SIGNAL(triggered()), this, SLOT(saveAs()));

  exitAct = new QAction(tr("E&xit"), this);
  exitAct->setShortcuts(QKeySequence::Quit);
  exitAct->setStatusTip(tr("Exit the application"));
  connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));


  aboutAct = new QAction(tr("&About"), this);
  aboutAct->setStatusTip(tr("Show the application's About box"));
  connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));

  aboutQtAct = new QAction(tr("About &Qt"), this);
  aboutQtAct->setStatusTip(tr("Show the Qt library's About box"));
  connect(aboutQtAct, SIGNAL(triggered()), qApp, SLOT(aboutQt()));

  togglePlayAct = new QAction(QIcon(":/images/media-playback-start.png"), tr("&Play"), this);
  togglePlayAct->setStatusTip(tr("Start / Pause simulation"));
  togglePlayAct->setShortcut (Qt::Key_Space);
  connect(togglePlayAct, SIGNAL(triggered()), this, SLOT(togglePlay()));

  stopAct = new QAction(QIcon(":/images/media-playback-stop.png"), tr("&Stop"), this);
  stopAct->setStatusTip(tr("Stop simulation"));
  stopAct->setShortcut (Qt::Key_Backspace);
  connect(stopAct, SIGNAL(triggered()), this, SLOT(stop()));

  restartAct = new QAction(QIcon(":/images/media-skip-backward.png"), tr("&Restart"), this);
  restartAct->setStatusTip(tr("Restart simulation"));
  connect(restartAct, SIGNAL(triggered()), this, SLOT(restart()));

}

void SpatialMainWindow::editGeometry()
{
  if (thread == NULL || thread->mpSimulator == NULL)
    return;
  bool wasRunning = thread->isRunning();
  if (wasRunning)
    pause();

  GeometryEditWidget widget(thread->mpSimulator);
  widget.setLastDir(lastDir);
  widget.exec();

  if (widget.needReload())
  {
    restart();
    return;
  }

  if (wasRunning)
    play();

}

void SpatialMainWindow::createMenus()
{
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(newAct);
  fileMenu->addAction(openAct);
  fileMenu->addAction(saveAct);
  fileMenu->addAction(saveAsAct);
  fileMenu->addSeparator();
  fileMenu->addAction(saveImageAct);
  fileMenu->addSeparator();
  fileMenu->addAction(exitAct);

  editMenu = menuBar()->addMenu(tr("&Edit"));
  editMenu->addAction(restartAct);
  editMenu->addAction(togglePlayAct);
  editMenu->addAction(stopAct);
  editMenu->addSeparator();
  editMenu->addAction(editGeometryAct);

  QMenu *windowsMenu = menuBar()->addMenu(tr("&Window"));
  windowsMenu->addAction(ui->dockAssignment->toggleViewAction());
  windowsMenu->addAction(ui->dockInfo->toggleViewAction());
  windowsMenu->addAction(ui->dockDose->toggleViewAction());
  windowsMenu->addAction(ui->dockParameters->toggleViewAction());


#ifdef USE_SBW_INTEGRATION
  // create and populate SBW menu
  mpSBWMenu = menuBar()->addMenu(tr("S&BW"));
  mpSBWAction = mpSBWMenu->menuAction();
#endif // USE_SBW_INTEGRATION

  
  menuBar()->addSeparator();

  helpMenu = menuBar()->addMenu(tr("&Help"));
  helpMenu->addAction(aboutAct);
  helpMenu->addAction(aboutQtAct);
}

void SpatialMainWindow::createToolBars()
{
  fileToolBar = addToolBar(tr("File"));
  fileToolBar->addAction(newAct);
  fileToolBar->addAction(openAct);
  fileToolBar->addAction(saveAct);

  editToolBar = addToolBar(tr("Edit"));
  editToolBar->addAction(restartAct);
  editToolBar->addAction(togglePlayAct);
  editToolBar->addAction(stopAct);

  editToolBar->addSeparator();
  txtStepsize = new QLineEdit("0.01", this);  
  txtStepsize->setMaximumWidth(50);

  connect(txtStepsize, SIGNAL(textEdited(const QString&)), this, SLOT(stepChanged(const QString&)));

  editToolBar->addWidget(new QLabel("Step Size: ", this));

  QAction *stepSizeAction = editToolBar->addWidget(txtStepsize);
  stepSizeAction->setStatusTip(tr("Set the timestep used by the spatial simulator"));

  editToolBar->addSeparator();

  txtUpdate = new QLineEdit("100", this);  
  txtUpdate->setMaximumWidth(50);

  connect(txtUpdate, SIGNAL(textEdited(const QString&)), thread, SLOT(updateFrequencyChanged(const QString&)));

  editToolBar->addWidget(new QLabel("Update Interval: ", this));
  QAction *updateAction = editToolBar->addWidget(txtUpdate);
  updateAction ->setStatusTip(tr("Set the frequency for updating the UI"));



}

void SpatialMainWindow::createStatusBar()
{
  statusBar()->showMessage(tr("Ready"));
}

void SpatialMainWindow::readSettings()
{
  QSettings settings("COPASI", "Spatial UI");
  QPoint pos = settings.value("pos", QPoint(200, 200)).toPoint();
  QSize size = settings.value("size", QSize(400, 400)).toSize();
  lastDir = settings.value("lastDir", QString()).toString();

  settings.beginGroup("dockAssignment");
  ui->dockAssignment->setFloating(settings.value("docked").toBool());
  ui->dockAssignment->resize(settings.value("size", QSize(1, 1)).toSize());
  ui->dockAssignment->move(settings.value("pos", QPoint(200, 200)).toPoint());
  addDockWidget((Qt::DockWidgetArea)settings.value("dockarea", Qt::RightDockWidgetArea).toInt(), ui->dockAssignment);
  ui->dockAssignment->setVisible(settings.value("visible", true).toBool());
  settings.endGroup();

  settings.beginGroup("dockInfo");
  ui->dockInfo->setFloating(settings.value("docked").toBool());
  ui->dockInfo->resize(settings.value("size", QSize(1, 1)).toSize());
  ui->dockInfo->move(settings.value("pos", QPoint(200, 200)).toPoint());
  addDockWidget((Qt::DockWidgetArea)settings.value("dockarea", Qt::RightDockWidgetArea).toInt(), ui->dockInfo);
  ui->dockInfo->setVisible(settings.value("visible", true).toBool());
  settings.endGroup();

  settings.beginGroup("dockDose");
  ui->dockDose->setFloating(settings.value("docked").toBool());
  ui->dockDose->resize(settings.value("size", QSize(1, 1)).toSize());
  ui->dockDose->move(settings.value("pos", QPoint(200, 200)).toPoint());
  addDockWidget((Qt::DockWidgetArea)settings.value("dockarea", Qt::RightDockWidgetArea).toInt(), ui->dockDose);
  ui->dockDose->setVisible(settings.value("visible", true).toBool());
  settings.endGroup();

  settings.beginGroup("dockParameters");
  ui->dockParameters->setFloating(settings.value("docked").toBool());
  ui->dockParameters->resize(settings.value("size", QSize(1, 1)).toSize());
  ui->dockParameters->move(settings.value("pos", QPoint(200, 200)).toPoint());
  addDockWidget((Qt::DockWidgetArea)settings.value("dockarea", Qt::RightDockWidgetArea).toInt(), ui->dockParameters);
  ui->dockParameters->setVisible(settings.value("visible", true).toBool());
  settings.endGroup();

  resize(size);
  move(pos);
}

void SpatialMainWindow::writeSettings()
{
  QSettings settings("COPASI", "Spatial UI");
  settings.setValue("pos", pos());
  settings.setValue("size", size());
  settings.setValue("lastDir", lastDir);

  settings.beginGroup("dockAssignment");
  settings.setValue("dockarea", dockWidgetArea(ui->dockAssignment));
  settings.setValue("docked", ui->dockAssignment->isFloating());
  settings.setValue("size", ui->dockAssignment->size());
  settings.setValue("pos", ui->dockAssignment->pos());
  settings.setValue("visible", ui->dockAssignment->isVisible());
  settings.endGroup();

  settings.beginGroup("dockInfo");
  settings.setValue("dockarea", dockWidgetArea(ui->dockInfo));
  settings.setValue("docked", ui->dockInfo->isFloating());
  settings.setValue("size", ui->dockInfo->size());
  settings.setValue("pos", ui->dockInfo->pos());
  settings.setValue("visible", ui->dockInfo->isVisible());
  settings.endGroup();

  settings.beginGroup("dockDose");
  settings.setValue("dockarea", dockWidgetArea(ui->dockDose));
  settings.setValue("docked", ui->dockDose->isFloating());
  settings.setValue("size", ui->dockDose->size());
  settings.setValue("pos", ui->dockDose->pos());
  settings.setValue("visible", ui->dockDose->isVisible());
  settings.endGroup();

  settings.beginGroup("dockParameters");
  settings.setValue("dockarea", dockWidgetArea(ui->dockParameters));
  settings.setValue("docked", ui->dockParameters->isFloating());
  settings.setValue("size", ui->dockParameters->size());
  settings.setValue("pos", ui->dockParameters->pos());
  settings.setValue("visible", ui->dockParameters->isVisible());
  settings.endGroup();

}

XMLNode* getAnnotationNode(Model* model, const std::string& ns)
{
  if (model == NULL || !model->isSetAnnotation()) return NULL;

  XMLNode* parent = model->getAnnotation();
  bool again = true;
  while(again)
  {
    again = false;
    for (size_t i = 0;i < parent->getNumChildren(); ++i)
    {
      XMLNode* current = &parent->getChild(i);
      if (current->getName() == "annotation")
      {
        again = true;
        parent = current;
        break;
      }

      if (current->hasNamespaceURI(ns))
        return current;
    }
  }

  return NULL;
}

void SpatialMainWindow::initializeDisplay(Model* model)
{

  displayItems.clear();
  ui->lstAssignments->clear();

  XMLNode* node = getAnnotationNode(model, SPATIAL_ANNOTATION_URL);

  if (node != NULL)
  {
    const XMLNode& update = node->getChild("update");
    if (update.getName() == "update")
    {
      txtStepsize->setText(update.getAttrValue("step").c_str());
      txtUpdate->setText(update.getAttrValue("freq").c_str());
    }

    const XMLNode& items = node->getChild("items");
    if (items.getName() == "items")
    {
      for (size_t i = 0;i < items.getNumChildren(); ++i)
      {
        const XMLNode& item = items.getChild(i);

        const string&id = item.getAttrValue("sbmlId");
        int index = ui->lstPalettes->findText(item.getAttrValue("name").c_str());

        DisplayItem* dItem = new DisplayItem( id, palettes[index]  );
        dItem->getPalette()->setMaxValue(QString(item.getAttrValue("max").c_str()).toDouble());
        displayItems.push_back(dItem);
        ui->lstAssignments->addItem(id.c_str());

      }
    }
  }

  #pragma region // fallback
  if (node == NULL)
  {
    if (lstSpecies->count() > 0 && palettes
      .size() > 0)
    {
      displayItems.push_back(new DisplayItem(
        lstSpecies->itemText(0).toStdString(),
        palettes[0]
      ));
      ui->lstAssignments->addItem(lstSpecies->itemText(0));
    }

    if (lstSpecies->count() > 1 && palettes
      .size() > 1)
    {
      displayItems.push_back(new DisplayItem(
        lstSpecies->itemText(1).toStdString(),
        palettes[1]
      ));
      ui->lstAssignments->addItem(lstSpecies->itemText(1));
    }

    if (lstSpecies->count() > 2 && palettes
      .size() > 2)
    {
      displayItems.push_back(new DisplayItem(
        lstSpecies->itemText(2).toStdString(),
        palettes[2]
      ));
      ui->lstAssignments->addItem(lstSpecies->itemText(2));
    }
  }
#pragma endregion

}

void SpatialMainWindow::saveDisplayToModel(Model* model)
{
  if (model == NULL) return;

  XMLNode node(XMLTriple("spatialInfo", SPATIAL_ANNOTATION_URL, ""), XMLAttributes());
  node.addAttr("xmlns", SPATIAL_ANNOTATION_URL);

  // save stepsize
  XMLNode update(XMLTriple("update", SPATIAL_ANNOTATION_URL, ""), XMLAttributes());
  update.addAttr("step", txtStepsize->text().toStdString());
  update.addAttr("freq", txtUpdate->text().toStdString());

  node.addChild(update);

  // save assignments
  if (displayItems.size() > 0)
  {
    XMLNode items(XMLTriple("items", SPATIAL_ANNOTATION_URL, ""), XMLAttributes());
    std::vector<DisplayItem*>::const_iterator it = displayItems.cbegin();

    while(it != displayItems.cend())
    {
      DisplayItem* current = *it;
      XMLNode item(XMLTriple("item", SPATIAL_ANNOTATION_URL, ""), XMLAttributes());
      item.addAttr("sbmlId", current->getId());
      QString file = current->getPalette()->getFilename();
      file = file.replace(qApp->applicationDirPath() +"/", QString(""));
      item.addAttr("name", file.toStdString());
      item.addAttr("max", QString::number(current->getPalette()->getMaxValue()).toStdString());
      items.addChild(item);
      ++it;
    }
    node.addChild(items);
  }

  model->removeTopLevelAnnotationElement("spatialInfo", SPATIAL_ANNOTATION_URL, false);
  model->appendAnnotation(&node);
}


void SpatialMainWindow::loadFromDocument(SBMLDocument* toLoad)
{

  Model* model = doc->getModel();

  if (model == NULL)
  {
    QMessageBox::critical(this, tr("Spatial UI"),
      tr("fatal errors while reading file %1\n%2.")
      .arg(curFile)
      .arg(doc->getErrorLog()->toString().c_str()));    
    return;
  }


  SpatialModelPlugin *modelPlugin = (SpatialModelPlugin*)model ->getPlugin("spatial");
  bool haveSpatial = modelPlugin != NULL && 
    doc->getPkgRequired("spatial");

  if (!haveSpatial)
  {
    QMessageBox::critical(this, tr("Spatial UI"),
      tr("It would seem that the given model '%1' does not use the spatial package. Please load one that does. \n")
      .arg(curFile)
      );
    return;
  }

  lstSpecies->clear();
  ui->lstDose->clear();

  maxX = modelPlugin->getGeometry()->getCoordinateComponent("x")->getBoundaryMax()->getValue() + 1;
  maxY = modelPlugin->getGeometry()->getCoordinateComponent("y")->getBoundaryMax()->getValue() + 1;


  for(unsigned int i = 0; i < model->getNumSpecies(); i++)
  {
    Species* species = model->getSpecies(i);
    SpatialSpeciesRxnPlugin* plugin = (SpatialSpeciesRxnPlugin*)species->getPlugin("spatial");
    if (plugin != NULL)
    {
      if (plugin->getIsSpatial())
      {
        lstSpecies->addItem(species->getId().c_str());
        ui->lstDose->addItem(species->getId().c_str());
      }
    }
  }

  initializeDisplay(model);

  restart();
  //thread->mpSimulator->flipVolumeOrder();
  //restart();

  fillParameters();

}

void SpatialMainWindow::fillParameters()
{
  if (thread == NULL || thread->mpSimulator == NULL) return;
  const Model* model = thread->mpSimulator->getModel();
  if (model == NULL) return;
  bool hideBC = ui->chkHideBC->isChecked();
  bool hideDiff = ui->chkHideDiff->isChecked();
  ui->tblParameters->clear();
  ui->tblParameters->setRowCount(model->getNumParameters());
  ui->tblParameters->setColumnCount(1);
  ui->tblParameters->setHorizontalHeaderLabels(QStringList() << "Value");
#if (QT_VERSION < QT_VERSION_CHECK(5, 0, 0))
  ui->tblParameters->horizontalHeader()->setResizeMode(0, QHeaderView::Stretch);
  ui->tblParameters->verticalHeader()->setResizeMode(QHeaderView::ResizeToContents);
#endif
  int count = 0;
  for (unsigned int i =0; i < model->getNumParameters(); ++i)
  {
    const Parameter* param = model->getParameter(i);
    const SpatialParameterPlugin* plug = dynamic_cast<const SpatialParameterPlugin*>(param->getPlugin("spatial"));
    if (hideBC && plug != NULL && plug->getType() == SBML_SPATIAL_BOUNDARYCONDITION)
      continue;
    if (hideDiff && plug != NULL && plug->getType() == SBML_SPATIAL_DIFFUSIONCOEFFICIENT)
      continue;
    if (plug != NULL && plug->getType() == SBML_SPATIAL_SPATIALSYMBOLREFERENCE)
      continue;
    QTableWidgetItem* item = new QTableWidgetItem(param->getId().c_str()); 
    ui->tblParameters->setVerticalHeaderItem(count, item);
    ui->tblParameters->setItem(count, 0, new QTableWidgetItem( QString::number(param->getValue())));
    ++count;
  }
  ui->tblParameters->setRowCount(count);

  ui->tblParameters->repaint();
}

void SpatialMainWindow::loadFromString(const std::string& sbml)
{
  updating = true;
  stop();

  doc = readSBMLFromString(sbml.c_str());

  setCurrentFile("fromSBW.xml");
  loadFromDocument(doc);

  updating = false;

  statusBar()->showMessage(tr("File loaded"), 2000);
}

void SpatialMainWindow::loadFile(const QString &fileName)
{
  updating = true;
  stop();

  QFile file(fileName);
  if (!file.open(QFile::ReadOnly | QFile::Text)) {
    QMessageBox::warning(this, tr("Application"),
      tr("Cannot read file %1:\n%2.")
      .arg(fileName)
      .arg(file.errorString()));
    return;
  }

  lastDir = QFileInfo(fileName).absoluteDir().absolutePath();
  doc = readSBMLFromFile(fileName.toStdString().c_str());

  setCurrentFile(fileName);
  loadFromDocument(doc);

  updating = false;

  statusBar()->showMessage(tr("File loaded"), 2000);
}

void SpatialMainWindow::parameterChanged(int row, int column)
{
  if (updating) return;
  bool wasRunning = !thread->mPaused;

  if (row < 0 || row > (int) doc->getModel()->getNumParameters())
    return;

  double newValue = ui->tblParameters->item(row, column)->text().toDouble();

  pause();

  const std::string& id = ui->tblParameters->
    verticalHeaderItem(row)->text().toStdString();

  thread->mpSimulator->setParameter(id, newValue);

  if (wasRunning) play();

}

void SpatialMainWindow::restart()
{
  stop();

  if (thread != NULL)
  {
    if (thread->isRunning())
    {
      thread->stop();
    }

    thread->mpSimulator = new SpatialSimulator(doc, maxX, maxY);
    thread->mStep = txtStepsize->text().toDouble();
    thread->mTime = 0;
    thread->start();

  } 
  else
  {
    thread = new SimulationThread();
    thread->mpDisplayItems = &displayItems;
    connect(thread, SIGNAL(visualize(double, const QImage&)), this, SLOT(visualize(double, const QImage&)));

  }

  togglePlay();

}

void SpatialMainWindow::stepChanged(const QString &newStep)
{
  if (thread != NULL)
    thread->mStep = newStep.toDouble();
}

bool SpatialMainWindow::saveImageFile(const QString &fileName)
{
  if (thread != NULL)
  {
    thread->getImage().save(fileName);
  }
  return true;
}

bool SpatialMainWindow::saveFile(const QString &fileName)
{
  if (thread != NULL)
  {
    saveDisplayToModel(doc->getModel());
    writeSBMLToFile(doc, fileName.toStdString().c_str());    
    setCurrentFile(fileName);
  }
  return true;
}

void SpatialMainWindow::setCurrentFile(const QString &fileName)
{
  curFile = fileName;
  setWindowModified(false);

  QString shownName = curFile;
  if (curFile.isEmpty())
    shownName = "untitled.xml";
  setWindowFilePath(shownName);
}

QString SpatialMainWindow::strippedName(const QString &fullFileName)
{
  return QFileInfo(fullFileName).fileName();
}

void SpatialMainWindow::pause()
{
  if (thread == NULL ) return;
  if (!thread->mPaused )
  {
    thread->mPaused = true;
    togglePlayAct->setText(tr("Start"));
    togglePlayAct->setIcon(QIcon(":/images/media-playback-start.png"));
  }
}

void SpatialMainWindow::play()
{
  if (thread == NULL ) return;
  if (thread->mPaused)
  {
    thread->mPaused = false;
    togglePlayAct->setText(tr("Pause"));
    togglePlayAct->setIcon(QIcon(":/images/media-playback-pause.png"));
  }
}

void SpatialMainWindow::togglePlay()
{
  if (thread == NULL ) return;

  if (!thread->mPaused)
  {
    pause();
  }
  else 
  {
    play();
  }
}

void SpatialMainWindow::stop()
{
  if (thread == NULL) return;

  thread->mPaused = true;

  togglePlayAct->setText(tr("Start"));
  togglePlayAct->setIcon(QIcon(":/images/media-playback-start.png"));
  if (thread != NULL)
  {
    thread->stop();
  }
}

void SpatialMainWindow::visualize(double time, const QImage& image)
{
  statusBar()->showMessage(tr("Time: %1").arg(time ));
  ui->lblImage->setPixmap(QPixmap::fromImage(image));
}

void SpatialMainWindow::resizeEvent( QResizeEvent *e )
{
  QWidget::resizeEvent( e );
}

void SpatialMainWindow::showEvent( QShowEvent *e )
{
  QWidget::showEvent( e );
}


#ifdef USE_SBW_INTEGRATION
// Create 2 custom events, one containing the filename to an SBML document to be loaded
// into COPASI
SpatialMainWindow::QSBWSBMLEvent::QSBWSBMLEvent(const std::string & SBMLModel):
  QEvent((QEvent::Type)65433),
  mSBML(SBMLModel)
{}

const std::string & SpatialMainWindow::QSBWSBMLEvent::getSBMLModel() const
{return mSBML;}

SpatialMainWindow::QSBWShutdownEvent::QSBWShutdownEvent():
  QEvent((QEvent::Type)65434)
{}

void SpatialMainWindow::registerMethods(SystemsBiologyWorkbench::MethodTable< SpatialMainWindow > & table)
{
  table.addMethod(&SpatialMainWindow::sbwAnalysis,
                  "void doAnalysis(string)",
                  false,
                  "Imports a given SBML Model into CopasiUI.");

  table.addMethod(&SpatialMainWindow::sbwGetSBML,
                  "string getSBML()",
                  false,
                  "Retrieves the currently in CopasiUI loaded model in SBML format.");
}

//virtual
void SpatialMainWindow::onShutdown()
{
  QApplication::postEvent(this, new QSBWShutdownEvent());
}

void SpatialMainWindow::customEvent(QEvent * event)
{
  // handle the file event, that is import the SBML file
  switch ((int) event->type())
    {
      case 65433:

        try
          {
            QSBWSBMLEvent *sbwEvent = static_cast< QSBWSBMLEvent *>(event);
            loadFromString(sbwEvent->getSBMLModel());
          }
        catch (...)
          {}

        break;

      case 65434:

        if (!mSBWIgnoreShutdownEvent)
          exit(0);

        mSBWIgnoreShutdownEvent = true;
        break;
    }
}

void SpatialMainWindow::sbwConnect()
{
  delete mpSBWModule;
  mpSBWModule = NULL;

  try
    {
      // Let us define how COPASI will look to the rest of SBW
      std::string FullName("Spatial SBML");

      // By belonging to the Analysis category, we tell all other modules that
      // COPASI can take SBML files and do *something* with them
      std::string Category("Analysis");
      std::string Description("Spatial Analyzer - Imports an SBML model into SpatialSBML");

      mpSBWModule =
        new SystemsBiologyWorkbench::ModuleImpl(FullName, FullName,
            SystemsBiologyWorkbench::UniqueModule,
            Description);

      mpSBWModule->addServiceObject(FullName, FullName, Category, this, Description);

      // this lets SBW ask COPASI to shut down
      SBW::addListener(this);

      // here we start the SBW services and give over to QT's main loop
      mpSBWModule->enableModuleServices();
    }

  catch (...)
    {      
      delete mpSBWModule;
      mpSBWModule = NULL;
    }

  // Update the SBW Menu
  sbwRefreshMenu();
}

void SpatialMainWindow::sbwDisconnect()
{
  if (mpSBWModule != NULL)
    {
      try
        {
          SBWLowLevel::disconnect();
          delete mpSBWModule;
        }

      catch (...)
        {}

      mpSBWModule = NULL;
    }
}

void SpatialMainWindow::sbwRegister()
{  
  if (mpSBWModule != NULL)
    {
      try
        {
          // Set the commandline so that SBW knows how to call us.
          const std::string& Self = QCoreApplication::arguments().at(0).toStdString();

          mpSBWModule->setCommandLine(Self.c_str());
          mpSBWModule->registerModule();

          sbwRefreshMenu();
        }

      catch (...)
        {}
    }
}

void SpatialMainWindow::sbwUnregister(const std::string & moduleName) const
{
  try
    {
      SystemsBiologyWorkbench::Module Module = SBW::getModuleInstance("BROKER");
      SystemsBiologyWorkbench::Service Service = Module.findServiceByName("BROKER");

      DataBlockWriter Arguments;
      Arguments.add(moduleName);
      Service.getMethod("void unregisterModule(string)").call(Arguments);
    }

  catch (...)
    {}
}

// get a list of all SBW analyzers and stick them into a menu
void SpatialMainWindow::sbwRefreshMenu()
{
  if (mpSBWMenu == NULL) 
    return;

  bool Visible = true;
  bool IsSBWRegistered = false;

  mSBWAnalyzerModules.clear();
  mSBWAnalyzerServices.clear();
  mpSBWMenu->clear();
  mSBWActionMap.clear();

  if (mpSBWActionGroup != NULL)
    {
      disconnect(mpSBWActionGroup, SIGNAL(triggered(QAction *)), this, SLOT(sbwSlotMenuTriggered(QAction *)));
      mpSBWActionGroup->deleteLater();
      mpSBWActionGroup = NULL;
    }

  mpSBWActionGroup = new QActionGroup(this);
  connect(mpSBWActionGroup, SIGNAL(triggered(QAction *)), this, SLOT(sbwSlotMenuTriggered(QAction *)));

  QAction * pAction;

  try
    {
      std::vector<DataBlockReader> Services = sbwFindServices("Analysis", true);
      std::vector<DataBlockReader>::const_iterator it = Services.begin();
      std::vector<DataBlockReader>::const_iterator end = Services.end();

      QMap< QString, int > SortedNames;
      QStringList ModuleList;
      QStringList ServiceList;

      std::string Self = QCoreApplication::arguments().at(0).toStdString();

      int i = 0;

      for (; it != end; ++it)
        {
          SystemsBiologyWorkbench::ServiceDescriptor Service(*it);
          SystemsBiologyWorkbench::ModuleDescriptor Module = Service.getModuleDescriptor();

          std::string ModuleName = Module.getName();
          std::string ServiceName = Service.getName();
          std::string MenuName = Service.getDisplayName();

          // Check whether the registered service is provided by Spatial SBML
          if (ServiceName.compare(0, 12, "Spatial SBML") == 0)
            {
              std::string CommandLine = Module.getCommandLine();

              // Check whether the registered module points to current CopasiUI
              if (CommandLine.compare(0, Self.length(), Self) == 0)
                {
                  // Check whether the versions match
                  if (ServiceName != "Spatial SBML")
                    {
                      // We transparently update the registration information
                      sbwUnregister(ModuleName);
                      return sbwRegister();
                    }

                  IsSBWRegistered = true;
                  continue;
                }
              else
                {
                  // Check whether the CommandLine is still valid
                  std::string::size_type Length = CommandLine.find(" -sbwmodule");

                  if (Length != std::string::npos)
                    {
                      CommandLine = CommandLine.substr(0, Length);
                    }

                  if (!QFile(CommandLine.c_str()).exists())
                    {
                      sbwUnregister(ModuleName);
                      continue;
                    }
                }
            }

          SortedNames[MenuName.c_str()] = i++;
          ModuleList.append(ModuleName.c_str());
          ServiceList.append(ServiceName.c_str());
        }

      // Add the option to register in SBW
      if (!IsSBWRegistered)
        {
          pAction = new QAction("Register", mpSBWActionGroup);
          mpSBWMenu->addAction(pAction);
          mSBWActionMap[pAction] = SortedNames.size();

          mpSBWMenu->addSeparator();
        }

      QMap< QString, int >::const_iterator itMap = SortedNames.begin();
      QMap< QString, int >::const_iterator endMap = SortedNames.end();

      for (i = 0; itMap != endMap; ++itMap, i++)
        {
          mSBWAnalyzerModules.append(ModuleList[itMap.value()]);
          mSBWAnalyzerServices.append(ServiceList[itMap.value()]);

          pAction = new QAction(itMap.key(), mpSBWActionGroup);
          mpSBWMenu->addAction(pAction);
          mSBWActionMap[pAction] = i;
        }

      if (mSBWAnalyzerModules.empty())
        Visible = false;
    }

  catch (...)
    {
      Visible = false;
    }

  if (!Visible)
    menuBar()->removeAction(mpSBWAction);

  return;
}

void SpatialMainWindow::sbwSlotMenuTriggeredFinished(bool success)
{
}

void SpatialMainWindow::sbwSlotMenuTriggered(QAction * pAction)
{

  mSBWActionId = mSBWActionMap[pAction];

  if (mSBWActionId == mSBWAnalyzerModules.size())
    {
      sbwRegister();
    }
  else
    {
      try
        {
          int nModule = SBWLowLevel::getModuleInstance(mSBWAnalyzerModules[mSBWActionId].toStdString().c_str());
          int nService = SBWLowLevel::moduleFindServiceByName(nModule, mSBWAnalyzerServices[mSBWActionId].toStdString().c_str());
          int nMethod = SBWLowLevel::serviceGetMethod(nModule, nService, "void doAnalysis(string)");

          DataBlockWriter args;
          args << writeSBMLToString(doc);
          SBWLowLevel::methodSend(nModule, nService, nMethod, args);
        }

      catch (SBWException * pE)
        {
          QMessageBox::critical(this, "SBW Error",
                                 pE->getMessage().c_str(),
                                 QMessageBox::Ok | QMessageBox::Default,
                                 QMessageBox::NoButton);
        }
      
    }
}

std::vector< DataBlockReader > SpatialMainWindow::sbwFindServices(const std::string & category,
    const bool & recursive)
{
  std::vector< DataBlockReader > result;

  try
    {
      DataBlockWriter oArguments;
      oArguments.add(category);
      oArguments.add(recursive);

      Module oModule = SBW::getModuleInstance("BROKER");
      Service oService = oModule.findServiceByName("BROKER");
      oService.getMethod("{}[] findServices(string, boolean)").call(oArguments) >> result;
    }

  catch (...)
    {
      result.clear();
    }

  return result;
}

// Here we get an SBML document
SystemsBiologyWorkbench::DataBlockWriter SpatialMainWindow::sbwAnalysis(SystemsBiologyWorkbench::Module /*from*/, SystemsBiologyWorkbench::DataBlockReader reader)
{
  try
    {
      std::string sSBMLModel;
      reader >> sSBMLModel;

      QSBWSBMLEvent *event = new QSBWSBMLEvent(sSBMLModel);
      QApplication::postEvent(this, event);
    }

  catch (...)
    {
      throw new SystemsBiologyWorkbench::SBWApplicationException("Error in doAnalysis");
    }

  // and yes ... every SBW method has to return something
  return SystemsBiologyWorkbench::DataBlockWriter();
}

SystemsBiologyWorkbench::DataBlockWriter SpatialMainWindow::sbwGetSBML(SystemsBiologyWorkbench::Module /*from*/, SystemsBiologyWorkbench::DataBlockReader /*reader*/)
{
  QMutexLocker Locker(&mSBWMutex);
  mSBWCallFinished = false;


  if (!mSBWCallFinished)
    {
      mSBWWaitSlot.wait(&mSBWMutex);
    }

  SystemsBiologyWorkbench::DataBlockWriter result;

  if (mSBWSuccess)
    {
      try
        {
          // write the current model as SBML and return it
          result << writeSBMLToString(doc);

          return result;
        }

      catch (...)
        {
          throw new SystemsBiologyWorkbench::SBWApplicationException("Error getting the SBML.");
        }
    }

  throw new SystemsBiologyWorkbench::SBWApplicationException("Error getting the SBML.");

  // This will never be reached.
  return result;
}

void SpatialMainWindow::sbwSlotGetSBMLFinished(bool success)
{
  QMutexLocker Locker(&mSBWMutex);

  mSBWCallFinished = true;
  mSBWSuccess = success;

  mSBWWaitSlot.wakeAll();
}
#else
//void SpatialMainWindow::slotSBWMenuTriggered(QAction * /* pAction */) {}
void SpatialMainWindow::customEvent(QEvent * /* event */) {}

#endif // USE_SBW_INTEGRATION