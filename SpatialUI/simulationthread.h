#ifndef SIMULATION_THREAD_H
#define SIMULATION_THREAD_H

#include <QThread>
#include <QImage>
#include <vector>
#include <string>

class SpatialSimulator;
class DisplayItem; 

class SimulationThread : public QThread
{
  Q_OBJECT

public:

  SimulationThread();
  virtual ~SimulationThread();

  const std::vector<std::pair<std::string, double> > getConcentrationsAt(int x, int y) const;
  const QImage& getImage() const { return *image; }
  void run();

  SpatialSimulator* mpSimulator;
  std::vector<DisplayItem*>* mpDisplayItems;
  double mTime;
  double mStep;

  bool mPaused;

signals:
  void visualize(double time, const QImage& image);

public slots:
    void stop();
    void step();
    void step(double stepSize);
    void applyDose(int x, int y, const QString &id, double value);
    void updateFrequencyChanged(const QString& number);


private:
  void updateImage();
  void updateUI();
  qint64 lastUpdated;
  int updateFrequency; 
  bool mStopped;
  QImage *image;

};

#endif //SIMULATION_THREAD_H