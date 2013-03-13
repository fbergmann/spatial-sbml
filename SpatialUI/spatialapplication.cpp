#include <QFileOpenEvent>
#include <QString>

#include "spatialapplication.h"
#include "spatialmainwidget.h"

SpatialApplication::SpatialApplication(int & argc, char ** argv):
  QApplication(argc, argv),
  mainWindow(NULL),
  file(),
  starting(true)
{
}

SpatialApplication::~SpatialApplication()
{
}

// virtual
bool SpatialApplication::event(QEvent * event)
{
  switch (event->type())
    {
      case QEvent::FileOpen:

        if (starting)
          {
            file = static_cast<QFileOpenEvent *>(event)->file();
          }
        else
          {
            // need to take the new file, otherwise whenever the application
            // is open we will re-open the first file that was supposed to be
            // opened.
            file = static_cast<QFileOpenEvent *>(event)->file();
            mainWindow->loadFile(file);
          }

        event->accept();
        return true;

        break;

      default:
        break;
    }

  return QApplication::event(event);
}

void SpatialApplication::setMainWindow(SpatialMainWindow * mainWindow)
{
  this->mainWindow = mainWindow;

  processEvents();

  starting = false;

  if (file.isNull() ||  file.isEmpty())
  return;
  mainWindow->loadFile(file);
}
