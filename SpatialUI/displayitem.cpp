#include "displayitem.h"
#include "concentrationpalette.h"


DisplayItem::DisplayItem( const std::string& sbmlId, const QString& fileName)
  : mSBMLId(sbmlId), mPalette (NULL), needToDelete(true)
{
  mPalette = new ConcentrationPalette(fileName);
}

DisplayItem::DisplayItem( const std::string& sbmlId, ConcentrationPalette *palette)
  : mSBMLId(sbmlId), mPalette (palette), needToDelete(false)
{
  
}

DisplayItem::~DisplayItem()
{
  if (needToDelete && mPalette  != NULL)
  {
    delete mPalette ;
    mPalette  = NULL;
  }
}

QRgb 
DisplayItem::getColorForConcentration(double concentration) const
{
  return mPalette ->getColorForConcentration(concentration);
}

QRgb 
DisplayItem::getBlendedColorForConcentration(QRgb color,double concentration) const
{
  return mPalette ->getBlendedColorForConcentration(color,concentration);
}

const char* 
DisplayItem::getId() const 
{
  return mSBMLId.c_str(); 
}

void DisplayItem::setPalette(ConcentrationPalette* palette)
{
  if (needToDelete && mPalette  != NULL)
  {
    delete mPalette ;
    mPalette  = NULL;
  }
  needToDelete = false;
  mPalette = palette;
  
}

const ConcentrationPalette* 
DisplayItem::getPalette() const 
{ 
  return mPalette ; 
}

ConcentrationPalette* 
DisplayItem::getPalette() 
{ 
  return mPalette ; 
}
