#ifndef SPATIAL_MAIN_WIDGET
#define SPATIAL_MAIN_WIDGET

#ifdef USE_SBW_INTEGRATION


# define WIN32_LEAN_AND_MEAN
# include <SBW/SBW.h>
# undef DELETE
# undef ERROR
# undef TRUE
# undef FALSE
# undef min
using namespace SystemsBiologyWorkbench;

# include <QApplication>
# include <QEvent>
# include <QMutex>
# include <QWaitCondition>
#else
class SBWListener;
#endif  // USE_SBW_INTEGRATION


#include <QMainWindow>
#include <QMap>
#include <QActionGroup>
#include <vector>
#include <string>




class QAction;
class QMenu;
class QGraphicsView;
class QGraphicsScene;
class QGraphicsPixmapItem;
class QListWidgetItem;
class QImage;
class QTimer;
class SpatialSimulator;
class SBMLDocument;
class Model;
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
#ifdef USE_SBW_INTEGRATION
  // A SBW listener can catch messages from SBW ... used here to allow COPASI to be shut down
  , public SBWListener
#endif // USE_SBW_INTEGRATION
{
  Q_OBJECT

public:
  SpatialMainWindow();
  bool loadFromDocument(SBMLDocument* doc);
protected:
  // SBW: handle the custom events
  void customEvent(QEvent *);

  void closeEvent(QCloseEvent *event);

public slots:
  void loadFile(const QString &fileName);
  void loadFromString(const std::string& sbml);
  void visualize(double time, const QImage& image);
  
private slots:;
  void newFile();
  void open();
  bool save();
  bool saveImage();
  bool saveAs();
  void about();
  void togglePlay();
  void pause();
  void play();
  void stop();
  void restart();
  void editGeometry();
  void exportConcentration();
  void importConcentration();
  void toggledHideBC(bool);
  void showAll ( );
  void hideAll ( );


  void itemChanged ( QListWidgetItem * );

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
  void fillParameters();
  void createActions();
  void createMenus();
  void createToolBars();
  void createStatusBar();
  void readSettings();
  void writeSettings();
  bool saveFile(const QString &fileName);
  bool saveImageFile(const QString &fileName);

  void initializeDisplay(Model* model);
  void saveDisplayToModel(Model* model);

  void setCurrentFile(const QString &fileName);
  QString strippedName(const QString &fullFileName);

  ConcentrationPalette* getPalette(const QString& name);
  ConcentrationPalette* getPalette(size_t index);

  DisplayItem* getDisplayItem(const QString& name);
  DisplayItem* getDisplayItem(size_t index);

  QListWidgetItem *getListItem(const QString& species);

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
  QAction *saveImageAct;
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


#ifdef USE_SBW_INTEGRATION
public:
  /**
   * This event is triggered by SBW asking COPASI to import an SBML document provided as a string
   */
  class QSBWSBMLEvent : public QEvent
  {
  public:
    /**
     * Constructor
     * @param const std::string & SBMLModel
     */
    QSBWSBMLEvent(const std::string & SBMLModel);

    /**
     * Retrieve the SBML model
     * @return const std::string & SBMLModel
     */
    const std::string & getSBMLModel() const;

  private:
    /**
     * A string holding the SBML model
     */
    std::string mSBML;
  };

  /**
   * This event is triggered by SBW asking COPASI shut down.
   */
  class QSBWShutdownEvent: public QEvent
  {
  public:
    QSBWShutdownEvent();
  };

  // We expose 2 methods to SBW, one to load an SBML file
  SystemsBiologyWorkbench::DataBlockWriter sbwAnalysis(SystemsBiologyWorkbench::Module from,
      SystemsBiologyWorkbench::DataBlockReader reader);

  // and another to return the SBML file COAPSI is currently working with
  SystemsBiologyWorkbench::DataBlockWriter sbwGetSBML(SystemsBiologyWorkbench::Module from,
      SystemsBiologyWorkbench::DataBlockReader reader);

  // This method must not be renamed as SBW calls it by name
  void registerMethods(SystemsBiologyWorkbench::MethodTable<SpatialMainWindow> & table);

  // as part of the SBWListener we tell SBW here, that we want to react on the shutdown event
  virtual void onShutdown();

private:
  /**
   * Connect to SBW
   */
  void sbwConnect();

  /**
   * Disconnect to SBW
   */
  void sbwDisconnect();

  /**
   * Register COPASI as a module ins SBW
   */
  void sbwRegister();

  /**
   * Unregister a module in SBW
   * @param const std::string & moduleName
   */
  void sbwUnregister(const std::string & moduleName) const;

  /**
   * Refresh the SBW menu.
   */
  void sbwRefreshMenu();

  /**
   * Retrieve the list of all services from the SBW broker
   * @param const std::string & category
   * @param const bool & recursive
   * @return std::vector< SystemsBiologyWorkbench::DataBlockReader > services
   */
  std::vector< SystemsBiologyWorkbench::DataBlockReader > sbwFindServices(const std::string & category,
      const bool & recursive);

protected slots:
  void sbwSlotMenuTriggered(QAction * pAction);
  void sbwSlotMenuTriggeredFinished(bool success);
  void sbwSlotGetSBMLFinished(bool success);

private:
  /**
   * The SBW module which handles the interaction with the SBW broker
   */
  SystemsBiologyWorkbench::ModuleImpl * mpSBWModule;

  /**
   * A list of SBW analyzer modules
   */
  QStringList mSBWAnalyzerModules;

  /**
   * A list of the corresponding SBW services
   */
  QStringList mSBWAnalyzerServices;

  /**
   * Map between actions and the index of SBW modules and services
   */
  QMap< QAction *, int > mSBWActionMap;

  /**
   * A group containing all actions of the SBW menu
   */
  QActionGroup * mpSBWActionGroup;

  /**
   * The SBW menu
   */
  QMenu * mpSBWMenu;

  /**
   * The SBW Action
   */
  QAction * mpSBWAction;

  /**
   * This variable indicates whether COPASI is to ignore SBW shutdown events
   */
  bool mSBWIgnoreShutdownEvent;

  QMutex mSBWMutex;

  QWaitCondition mSBWWaitSlot;

  bool mSBWCallFinished;

  bool mSBWSuccess;

  std::string mSBWDocumentString;

  QStringList::size_type mSBWActionId;

#endif // USE_SBW_INTEGRATION
  
};


#endif //SPATIAL_MAIN_WIDGET
