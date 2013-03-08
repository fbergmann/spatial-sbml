#ifndef SPATIAL_MAIN_WIDGET
#define SPATIAL_MAIN_WIDGET

#include <QMainWindow>
#include <vector>
#include <string>

class QAction;
class QMenu;
class QGraphicsView;
class QGraphicsScene;
class QGraphicsPixmapItem;
class QImage;
class QTimer;
class SpatialSimulator;
class SBMLDocument;
class QLineEdit;
class QComboBox;
class ConcentrationPalette;
class DisplayItem;
class QResizeEvent;
class QShowEvent;
class SimulationThread;

namespace Ui {
  class SpatialMainWindow;
}

class SpatialMainWindow : public QMainWindow
{
  Q_OBJECT

public:
  SpatialMainWindow();

protected:
  void closeEvent(QCloseEvent *event);

public slots:
  void loadFile(const QString &fileName);
  void visualize(double time, const QImage& image);

private slots:;
  void newFile();
  void open();
  bool save();
  bool saveAs();
  void about();
  void togglePlay();
  void pause();
  void play();
  void stop();
  void restart();
  void editGeometry();

  void updatePosition(int x, int y);

  void stepChanged(const QString &newStep);
  void parameterChanged(int row, int column);

  void selectedSpeciesChanged(int row);
  void addSpecies();
  void removeSpecies();
  void palatteIndexChanged(int index);

public:
  virtual void resizeEvent( QResizeEvent *e );
  virtual void showEvent( QShowEvent *e );

private:
  void createActions();
  void createMenus();
  void createToolBars();
  void createStatusBar();
  void readSettings();
  void writeSettings();
  bool saveFile(const QString &fileName);

  void setCurrentFile(const QString &fileName);
  QString strippedName(const QString &fullFileName);

  QString curFile;
  QString lastDir;    

  QMenu *fileMenu;
  QMenu *editMenu;
  QMenu *helpMenu;
  QToolBar *fileToolBar;
  QToolBar *editToolBar;
  QAction *newAct;
  QAction *openAct;
  QAction *saveAct;
  QAction *saveAsAct;
  QAction *exitAct;
  QAction *togglePlayAct;
  QAction *stopAct;     
  QAction *restartAct;
  QAction *aboutAct;
  QAction *aboutQtAct;
  QAction *editGeometryAct;

  QLineEdit *txtStepsize;
  QLineEdit *txtUpdate;
  QComboBox* lstSpecies;

  SimulationThread *thread;
  SBMLDocument* doc;
  
  bool updating;

  int pickerX;
  int pickerY;

  int maxX;
  int maxY;

  std::vector<ConcentrationPalette*> palettes;
  std::vector<DisplayItem*> displayItems;

  Ui::SpatialMainWindow *ui;

};


#endif //SPATIAL_MAIN_WIDGET
