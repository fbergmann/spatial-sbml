#ifndef DISPLAY_ITEM_H
#define DISPLAY_ITEM_H

#include <string>
#include <QString>
#include <QRgb>

class ConcentrationPalette;

class DisplayItem 
{
public: 
  DisplayItem( const std::string& sbmlId, const QString& fileName);
  DisplayItem( const std::string& sbmlId, ConcentrationPalette *palette);
  ~DisplayItem();
  const char* getId() const;
  const ConcentrationPalette* getPalette() const;
  ConcentrationPalette* getPalette();
  void setPalette(ConcentrationPalette* palette);
  QRgb getBlendedColorForConcentration(QRgb color, double concentration) const;
  QRgb getColorForConcentration(double concentration) const;

  bool isVisible() const;
  void setVisible(bool);

private:
  std::string mSBMLId;
  ConcentrationPalette *mPalette;
  bool needToDelete;
  bool mVisible;
};


#endif //DISPLAY_ITEM_H