#ifndef CALCPDE_H_
#define CALCPDE_H_

void reversePolishInitial(vector<int> &indexList, reversePolishInfo *rpInfo, double *value, int numOfASTNodes, int Xindex, int Yindex, int Zindex, bool isAllArea);
void reversePolishRK(reactionInfo *rInfo, GeometryInfo *geoInfo, int Xindex, int Yindex, int Zindex, double dt, int m, int numOfReactants, bool isReaction);
//void reversePolishRK(reactionInfo *rInfo, GeometryInfo *geoi, Model *model, vector<variableInfo*> &varInfoList, int Xindex, int Yindex, int Zindex, double dt, int m);
void calcDiffusion(variableInfo *sInfo, double deltaX, double deltaY, double deltaZ, int Xindex, int Yindex, int Zindex, int m, double dt);
void cipCSLR(variableInfo *sInfo, double deltaX, double deltaY, double deltaZ, double dt, int Xindex, int Yindex, int Zindex, int dimension);
void calcBoundary(variableInfo *sInfo, double deltaX, double deltaY, double deltaZ, int Xindex, int Yindex, int Zindex, int m, int dimension);
void calcMemTransport(reactionInfo *rInfo, GeometryInfo *geoInfo, normalUnitVector *nuVec, int Xindex, int Yindex, int Zindex, double dt, int m, double deltaX, double deltaY, double deltaZ, int dimension, int numOfReactants);
void calcMemDiffusion(variableInfo *sInfo, voronoiInfo *vorI, int Xindex, int Yindex, int Zindex, int m, double dt, int dimension);

#endif
