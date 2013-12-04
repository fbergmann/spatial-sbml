#include <QApplication>
#include <QtGui>
#include <QObject>
#include <QDir>
#include <QListWidget>

#include "spatialmainwidget.h"
#include "spatialsimulator.h"
#include "concentrationpalette.h"
#include "displayitem.h"
#include "ui_spatialmainwidget.h"
#include "simulationthread.h"
#include "geometryeditwidget.h"

#include <sbml/SBMLTypes.h>
#include <sbml/packages/spatial/extension/SpatialSpeciesRxnPlugin.h>
#include <sbml/packages/spatial/extension/SpatialModelPlugin.h>

#include <string>
#include <string.h>

using namespace std;
LIBSBML_CPP_NAMESPACE_USE;

SpatialMainWindow::SpatialMainWindow() : thread(NULL), updating(false), pickerX(0), pickerY(0), maxX(501), maxY(501)
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
  connect(ui->lstAssignments, SIGNAL(currentRowChanged(int)), this, SLOT(selectedSpeciesChanged(int)));
  connect(ui->lblImage, SIGNAL(positionChanged(int, int)), this, SLOT(updatePosition(int,int)));
  connect(ui->tblParameters, SIGNAL(cellChanged(int, int)), this, SLOT(parameterChanged(int,int)));

  createActions();
  createMenus();
  createToolBars();
  createStatusBar();

  readSettings();

  QDir dir (qApp->applicationDirPath());
  QFileInfoList foo = dir.entryInfoList(QStringList() << "*.txt");
  for (int i=0;i<foo.size();i++) {
    palettes.push_back(new ConcentrationPalette(qApp->applicationDirPath() +"/"+ foo[i].fileName()));
    ui->lstPalettes->addItem(foo[i].fileName());
  }

  setCurrentFile("");
  setUnifiedTitleAndToolBarOnMac(true);
  setWindowIcon(QIcon(":/images/ICON_Spatial_128x128.png"));

}

void SpatialMainWindow::updatePosition(int x, int y)
{
  pickerX = x;
  pickerY = y;
  ui->txtXPos->setText(QString::number(x));
  ui->txtYPos->setText(QString::number(y));

  int transformedX = x;
  int transformedY = maxY-1 - y;

  if (  (QApplication::mouseButtons().operator&((int)Qt::LeftButton) == (int) Qt::LeftButton)  && 
    ui->cmdEnableDose->isChecked() &&  thread != NULL
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
      currentSpecies.ascii(),
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

  bool wasRunning = thread->mPaused;
  if (wasRunning) pause();

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
  //if (curFile.isEmpty()) {
  return saveAs();
  //} else {
  //  return saveFile(curFile);
  //}
}

bool SpatialMainWindow::saveAs()
{
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save Image"), lastDir, tr("Image files (*.png *.jpeg *.jpg)"));
  if (fileName.isEmpty())
    return false;

  return saveFile(fileName);
}

void SpatialMainWindow::about()
{
  QMessageBox::about(this, tr("About Spatial UI"),
    tr("The <b>Spatial UI</b> application opens a file with the SBML spatial extension and allows to simulate it using Akira's simulator."));
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
  doc = readSBML(fileName.ascii());

  Model* model = doc->getModel();

  if (model == NULL)
  {
    QMessageBox::critical(this, tr("Spatial UI"),
      tr("fatal errors while reading file %1\n%2.")
      .arg(fileName)
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
      .arg(fileName)
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

  displayItems.clear();
  ui->lstAssignments->clear();


  if (lstSpecies->count() > 0 && palettes
    .size() > 0)
  {
    displayItems.push_back(new DisplayItem(
      lstSpecies->itemText(0).ascii(),
      palettes[0]
    ));
    ui->lstAssignments->addItem(lstSpecies->itemText(0));
  }



  restart();

  ui->tblParameters->clear();
  ui->tblParameters->setColumnCount(1);
  ui->tblParameters->setHorizontalHeaderLabels(QStringList() << "Value");
  ui->tblParameters->setRowCount(model->getNumParameters());
  ui->tblParameters->horizontalHeader()->setResizeMode(0, QHeaderView::Stretch);
  ui->tblParameters->verticalHeader()->setResizeMode(QHeaderView::ResizeToContents);
  for (unsigned int i =0; i < model->getNumParameters(); i++)
  {
    const Parameter* param = model->getParameter(i);
    QTableWidgetItem* item = new QTableWidgetItem(param->getId().c_str()); 
    ui->tblParameters->setVerticalHeaderItem(i, item);
    ui->tblParameters->setItem(i, 0, new QTableWidgetItem( QString::number(param->getValue())));
  }

  updating = false;

  setCurrentFile(fileName);
  statusBar()->showMessage(tr("File loaded"), 2000);
}

void SpatialMainWindow::parameterChanged(int row, int column)
{
  if (updating) return;
  bool wasRunning = thread->mPaused;

  if (row < 0 || row > (int) doc->getModel()->getNumParameters())
    return;

  double newValue = ui->tblParameters->item(row, column)->text().toDouble();

  if (wasRunning) pause();

  thread->mpSimulator->setParameter(doc->getModel()->getParameter(row)->getId(), newValue);

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
  }

  thread->mpSimulator = new SpatialSimulator(doc, maxX, maxY);
  thread->mStep = txtStepsize->text().toDouble();
  thread->mTime = 0;
  thread->start();
  togglePlay();

}

void SpatialMainWindow::stepChanged(const QString &newStep)
{
  if (thread != NULL)
    thread->mStep = newStep.toDouble();
}

bool SpatialMainWindow::saveFile(const QString &fileName)
{
  if (thread != NULL)
  {
    thread->getImage().save(fileName);
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
  if (!thread->mPaused )
  {
    thread->mPaused = true;
    togglePlayAct->setText(tr("Start"));
    togglePlayAct->setIcon(QIcon(":/images/media-playback-start.png"));
  }
}

void SpatialMainWindow::play()
{
  if (thread->mPaused)
  {
    thread->mPaused = false;
    togglePlayAct->setText(tr("Pause"));
    togglePlayAct->setIcon(QIcon(":/images/media-playback-pause.png"));
  }
}

void SpatialMainWindow::togglePlay()
{
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

