#ifndef SETINFO_H_
#define SETINFO_H_

void setCompartmentInfo(Model *model, vector<variableInfo*> &varInfoList);
void setSpeciesInfo(SBMLDocument *doc, vector<variableInfo*> &varInfoList, unsigned int volDimension, unsigned int memDimension, int Xindex, int Yindex, int Zindex);
void setParameterInfo(SBMLDocument *doc, vector<variableInfo*> &varInfoList, int Xdiv, int Ydiv, int Zdiv, double &Xsize, double &Ysize, double &Zsize, double &deltaX, double &deltaY, double &deltaZ, char *&xaxis, char *&yaxis, char *&zaxis);
void setReactionInfo(Model *model, vector<variableInfo*> &varInfoList, vector<reactionInfo*> &rInfoList, vector<reactionInfo*> &fast_rInfoList, vector<double*> freeConstList, int numOfVolIndexes);
void setRateRuleInfo(Model *model, vector<variableInfo*> &varInfoList, vector<reactionInfo*> &rInfoList, vector<double*> freeConstList, int numOfVolIndexes);
normalUnitVector* setNormalAngle(vector<GeometryInfo*> &geoInfoList, double Xsize, double Ysize, double Zsize, int dimension, int Xindex, int Yindex, int Zindex, int numOfVolIndexes);
voronoiInfo* setVoronoiInfo(normalUnitVector *nuVec, variableInfo *xInfo, variableInfo *yInfo, variableInfo *zInfo, vector<GeometryInfo*> &geoInfoList, double Xsize, double Ysize, double Zsize, int dimension, int Xindex, int Yindex, int Zindex, int numOfVolIndexes);

#endif
