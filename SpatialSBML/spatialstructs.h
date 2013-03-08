#ifndef SPATIAL_STRUCTS_H_
#define SPATIAL_STRUCTS_H_

using namespace std;
LIBSBML_CPP_NAMESPACE_USE;


typedef enum _materialTYpe {
  reactants = 0, products
}materialType;

typedef struct _reversePolishInfo {
  double **varList;
  double **deltaList;
  double **constList;
  int *opfuncList;
  int listNum;
}reversePolishInfo;

typedef struct _memIndexArray {
  int X;
  int Y;
}memIndexArray;

typedef struct _reactionInfo {
  const char* id;
  double *value;
  reversePolishInfo *rpInfo;
}reactionInfo;

typedef struct _analyticVolInfo {
  const char *compartmentId;
  const char *domainTypeId;
  const char *domainId;
  AdjacentDomains *adjacent0;
  AdjacentDomains *adjacent1;
  boundaryType *bt;
  bool isVol;
  bool implicit;
  double *isDomain;
  int *isBoundary;
  bool *isEdge;
  reversePolishInfo *rpInfo;
}analyticVolInfo;

typedef struct _boundaryCInfo {
  BoundaryCondition *bc;
  double *value;
  reversePolishInfo *rpInfo;
}boundaryCInfo;

typedef struct _bcInfoXYZ {
  const char *id;
  boundaryCInfo *Xp;
  boundaryCInfo *Xm;
  boundaryCInfo *Yp;
  boundaryCInfo *Ym;
  boundaryCInfo *Zp;
  boundaryCInfo *Zm;
}bcInfoXYZ;

typedef struct _diffCInfo {
  DiffusionCoefficient *diffc;
  double *value;
  reversePolishInfo *rpInfo;
}diffCInfo;

typedef struct _adCInfo {
  AdvectionCoefficient *adc;
  double *value;
  reversePolishInfo *rpInfo;
}adCInfo;

typedef struct _derivativeCIP {
  double *fx;
  double *fx_next;
  double *fy;
  double *fy_next;
  double *fz;
  double *fz_next;
}derivativeCIP;

typedef struct _variableInfo {
  Species *sp;
  Compartment *com;
  Parameter *para;
  const char* id;
  double *value;
  double *delta;
  bool inVol;
  bool hasAssignmentRule;
  boundaryType bType;
  reversePolishInfo *rpInfo;
  analyticVolInfo *avi;
  diffCInfo **dci;
  adCInfo **aci;
  bool isResolved;
  int dim;
  std::vector<_variableInfo*> dependence;
  derivativeCIP *dcip;
}variableInfo;

typedef struct _bcOfSpeciesInfo {
  const char* speciesId;
  variableInfo *bcXp;
  string bcXpType;
  variableInfo *bcXm;
  string bcXmType;
  variableInfo *bcYp;
  string bcYpType;
  variableInfo *bcYm;
  string bcYmType;
  variableInfo *bcZp;
  string bcZpType;
  variableInfo *bcZm;
  string bcZmType;
}bcOfSpeciesInfo;

#endif /* SPATIAL_STRUCTS_H_ */
