#include "simulationthread.h"
#include "spatialsimulator.h"
#include "displayitem.h"

#include <QDateTime>
#include <QImage>
#include <vector>

using namespace std;

SimulationThread::SimulationThread()
  : mpSimulator(NULL), mpDisplayItems(NULL), mPaused(false), lastUpdated(0), updateFrequency(200), image(NULL)
{
  QThread::setTerminationEnabled();
}

void SimulationThread::updateFrequencyChanged(const QString& number)
{
  updateFrequency = number.toInt();
}

SimulationThread::~SimulationThread()
{
  if (mpSimulator != NULL)
  {
    delete mpSimulator;
    mpSimulator = NULL;
  }

  if (image != NULL)
    {
      delete image; 
      image = NULL;
    }
}

void SimulationThread::stop()
{
  mStopped = true;
  QThread::quit();
  wait();
}

void SimulationThread::step()
{
  step(mStep);
}

void SimulationThread::step(double stepSize)
{
  if (mPaused || mpSimulator == NULL)
  {
    msleep(100);
  }
  else
  {
    mTime = mpSimulator->oneStep(mTime, stepSize);
  }  
  updateUI();
}

void SimulationThread::updateUI()
{
  long delta = QDateTime::currentMSecsSinceEpoch () - lastUpdated;
  if (delta < updateFrequency )
    return;

  //cout << "delta: " << delta << endl;

  // ascertain we have image of right size
  if (image == NULL || image->width() != mpSimulator->getXDim() || image->height() != mpSimulator->getYDim())
  {
    if (image != NULL)
    {
      delete image; 
      image = NULL;
    }
    image = new QImage(mpSimulator->getXDim(), mpSimulator->getYDim(), QImage::Format_RGB32);
  }

  
  updateImage();
  
  lastUpdated = QDateTime::currentMSecsSinceEpoch();
  
  emit visualize(mTime, *image);  

}

void SimulationThread::updateImage()
{
  if (image == NULL || mpSimulator == NULL || mpDisplayItems == NULL || mpDisplayItems->size() == 0)
    return;
  // fill image
  int length;
  double *xVal = mpSimulator->getX(length);
  double *yVal = mpSimulator->getY(length);
  
  double** allvalues = (double**) malloc(sizeof(double*)*mpDisplayItems->size());
  memset(allvalues, 0, sizeof(double*)*mpDisplayItems->size());
  for (unsigned int ind = 0; ind < mpDisplayItems->size(); ind++)
  {
    allvalues[ind] = mpSimulator->getVariable((*mpDisplayItems)[ind]->getId(), length);
  }
  
  int Xdiv =  mpSimulator->getXDim(); int Ydiv =  mpSimulator->getYDim();
  int Xindex = 2 * Xdiv - 1, Yindex = 2 * Ydiv - 1;
  int X, Y;
  int index;
  int stride = 2;
  for (Y = 0; Y < Yindex; Y += stride ) {
    for (X = 0; X < Xindex; X += stride ) {
      index = Y * Xindex + X;
      int x =  xVal[index];
      int y =  yVal[index];
      
      double current = allvalues[0][index];
      
      QRgb color = (*mpDisplayItems)[0]->getColorForConcentration(current);
      
      for (unsigned int ind = 1; ind < mpDisplayItems->size(); ind++)
      {
        current = allvalues[ind][index];
        color =  (*mpDisplayItems)[ind]->getBlendedColorForConcentration(color, current);
      }
      
      image->setPixel(x,Ydiv - 1-y,color );
    }
  }
  free(allvalues);
}

const std::vector<std::pair<std::string, double> > SimulationThread::getConcentrationsAt(int x, int y) const
{
  std::vector<std::pair<std::string, double> > result;

  if (mpDisplayItems == NULL || mpDisplayItems->size() == 0) return result;

  int index = mpSimulator->getIndexForPosition(x, y);
  if (index == -1) return result;
  
  int length;

  for (unsigned int ind = 0; ind < mpDisplayItems->size(); ind++)
  {
    string id = (*mpDisplayItems)[ind]->getId();
    std::pair<std::string, double> current (id,mpSimulator->getVariable(id, length)[index]);
    result.push_back(current);
  }

  return result;
}

void SimulationThread::applyDose(int x, int y,const QString &id, double value)
{
  if (mpDisplayItems == NULL || mpDisplayItems->size() == 0) return;

  bool wasRunning = mPaused;
  mPaused = true;
  int length;
  int index = mpSimulator->getIndexForPosition(x, y);
  if (index == -1) return;
  for (unsigned int ind = 0; ind < mpDisplayItems->size(); ind++)
  {
    if ((*mpDisplayItems)[ind]->getId() != id)
      continue;
    mpSimulator->getVariable(id.ascii(), length)[index] = value;
  }

  mPaused = wasRunning;
}

void SimulationThread::run()
{
  mStopped = false;
  while(!mStopped)
  {
    step();
  }
}

