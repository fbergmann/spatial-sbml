#include "concentrationpalette.h"

#include <QFile>
#include <QTextStream>
#include <QString>
#include <QRgb>

#include <vector>

using namespace std;

ConcentrationPalette::ConcentrationPalette(const QString& fileName, double maxValue)
  : mFileName(fileName), mMaxValue(maxValue)
{
  loadPalette(fileName);
}

QRgb ConcentrationPalette::blendColor(QRgb color1, QRgb color2)
{
  return qRgb(
    min(qRed(color1) + qRed(color2), 255),
    min(qGreen(color1) + qGreen(color2), 255),
    min(qBlue(color1) + qBlue(color2), 255)
    );

}

QRgb ConcentrationPalette::getBlendedColorForConcentration(QRgb color, double concentration) const
{
  return blendColor(color, getColorForConcentration(concentration));
}


void ConcentrationPalette::loadPalette(const QString& fileName)
{
  mFileName = fileName;
  mPalette.clear();
    QFile file(mFileName);

  file.open(QIODevice::ReadOnly);

  QTextStream stream( &file );

  QString line;
  while ((line = stream.readLine()) != "")
    mPalette.push_back(QColor(line));
  file.close();

}

const QString& 
  ConcentrationPalette::getFilename() const
{
  return mFileName;
}

size_t 
  ConcentrationPalette::getNumColors() const
{
  return mPalette.size();
}

double 
  ConcentrationPalette::getMaxValue() const
{
  return mMaxValue;
}

void 
  ConcentrationPalette::setMaxValue(double value)
{
  mMaxValue = value;
}

QRgb 
  ConcentrationPalette::getColorForConcentration(double concentration) const
{
  if (mPalette.empty()) return qRgb(0,0,0);

  double maxV = max(concentration, mMaxValue);

  int color = (int)
    ((concentration/maxV)*(mPalette.size() - 1));

  if (color < 0) color = 0;

  return mPalette[color].rgb();
}

