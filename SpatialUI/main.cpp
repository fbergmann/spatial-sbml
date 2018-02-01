
#include "spatialapplication.h"
#include "spatialmainwidget.h"

int main(int argc, char *argv[])
{
  Q_INIT_RESOURCE(spatialresources);

  SpatialApplication app(argc, argv);
  app.setOrganizationName("COPASI Team");
  app.setApplicationName("Spatial UI");
  SpatialMainWindow mainWin;
#if defined(Q_OS_SYMBIAN)
  mainWin.showMaximized();
#else
  mainWin.show();
  if (argc > 1)
    mainWin.loadFile(argv[1]);
#endif
  app.setMainWindow(&mainWin);
  return app.exec();
}