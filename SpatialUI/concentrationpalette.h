#ifndef CONCENTRATION_PALETTE_H
#define CONCENTRATION_PALETTE_H

#include <vector>
#include <QString>
#include <QColor>

class ConcentrationPalette 
{
public:

  ConcentrationPalette(const QString& fileName, double maxValue = 6.0);

  static QRgb blendColor(QRgb color1, QRgb color2);

  const QString& getFilename() const;
  QString getBasename() const;
  size_t getNumColors() const;

  double getMaxValue() const;
  void setMaxValue(double value);

  void loadPalette(const QString& fileName);

  QRgb getBlendedColorForConcentration(QRgb color, double concentration) const;
  QRgb getColorForConcentration(double concentration) const;

private:
  std::vector<QColor> mPalette;
  QString mFileName;
  double mMaxValue;
};


#endif // CONCENTRATION_PALETTE_H