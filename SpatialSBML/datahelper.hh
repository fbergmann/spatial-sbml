#ifndef DATAHELPER_HH
#define DATAHELPER_HH

#include <vector>
#include <string>
#include <fstream>
#include <math.h>

class DataHelper
{

private:
  std::vector<std::vector<double>> mData;
  int maxX;
  int maxY;
  bool valid;

public:
  DataHelper()
    : mData(0)
    , maxX(0)
    , maxY(0)
    , valid(false)
  {
  }

  DataHelper(const std::string& fileName)
    : mData(0)
    , maxX(0)
    , maxY(0)
    , valid(false)
  {
    initializeFromFile(fileName);
  }

  DataHelper(int numRows, int numCols, const double &value=double() ) 
    : mData(numRows)
    , maxX(numRows)
    , maxY(numCols)
    , valid(true)
  {
    if ( ((numCols==0) && (numRows!=0)) || ((numCols!=0) && (numRows==0)) )
    {
      maxX=0;
      maxY=0;
      mData.resize(maxX);
    }
    else
    {
      for (int i=0;i<maxX;++i)
        mData[i].resize(maxY);
      maxX=numRows;
      maxY=numCols;	    
    }
    for (int i=0;i<numRows;++i)
    {
      for (int j=0;j<numCols;++j)
        mData[i][j]=value;
    }

  }

  void resize(int numRows, int numCols, const double &value = double() )
  {
    mData.resize(numRows);
    for (size_t i=0;i<mData.size();++i)
    {
      mData[i].resize(numCols);
      for (size_t j=0;j<mData[i].size();++j)
        mData[i][j]=value;
    }
    maxX=numRows;
    maxY=numCols;
    valid = true;
  }


  void initializeFromFile(const std::string& fileName)
  {
    std::ifstream stream(fileName.c_str(), std::ios_base::binary);
    valid = stream.good();
    if (!valid) return;
    stream >> maxX >> maxY;
    mData.resize(maxX);
    for (int x = 0; x < maxX; ++x)
    {
      mData[x].resize(maxY);
      for (int y = 0; y < maxY; ++y)
      {
        stream >> mData[x][y];
      }
    }
    stream.close();
  }

  void writeToFile(const std::string& fileName)
  {
    std::ofstream data(fileName.c_str(), std::ios_base::binary);

    data << maxX << " ";
    data << maxY << std::endl;

    for (int x = 0; x < maxX; ++x)
    {
      for (int y = 0; y < maxY; ++y)
      {
        data << mData[x][y] << " ";
      }
      data << std::endl;
    }

    data.flush();
    data.close();
  }

  bool isValid() const { return valid; }

  double operator ()(double x, double y) const
  {
    int xpos = floor(x);
    if (xpos >= maxX) xpos = maxX - 1;
    if (xpos < 0) xpos = 0;
    int ypos = floor(y);
    if (ypos >= maxY) ypos = maxY - 1;
    if (ypos < 0) ypos = 0;

    return mData[xpos][ypos];
  }


  double &operator()(double x, double y)
  {
    int xpos = floor(x);
    if (xpos >= maxX) xpos = maxX - 1;
    if (xpos < 0) xpos = 0;
    int ypos = floor(y);
    if (ypos >= maxY) ypos = maxY - 1;
    if (ypos < 0) ypos = 0;
    return mData[xpos][ypos];
  }

  int Rows() const
  {
    return maxX;
  }

  int Cols() const
  {
    return maxY;
  }

};


#endif // DATAHELPER_HH
