#ifndef MYSTRUCT_H_
#define MYSTRUCT_H_

using namespace std;

typedef enum _materialTYpe {
	reactants = 0, products
}materialType;

typedef enum _boundaryIndex {
	Xmax = 0, Xmin, Ymax, Ymin, Zmax, Zmin
}boundaryIndex;

typedef enum _preDirection {
	N = 0, NE, E, SE, S, SW, W, NW
}preDirection;

typedef struct _reversePolishInfo {
	double **varList;
	double **deltaList;
	double **constList;
	int *opfuncList;
	int listNum;
}reversePolishInfo;

typedef struct _boundaryType {
	bool isBofXp;
	bool isBofXm;
	bool isBofYp;
	bool isBofYm;
	bool isBofZp;
	bool isBofZm;
}boundaryType;

typedef struct _adjacentMemElement {
	int index_x;
	int index_y;
	int index_z;
}adjacentMemElement;

typedef struct _GeometryInfo {
	const char *compartmentId;
	const char *domainTypeId;
	const char *domainId;
	vector <const char *> domainIdList;
	AdjacentDomains *adjacent0;
	AdjacentDomains *adjacent1;
	_GeometryInfo *adjacentGeo1;
	_GeometryInfo *adjacentGeo2;
	boundaryType *bType;
	bool isVol;
	bool implicit;
	int *isDomain;
	int *isBoundary;
	reversePolishInfo *rpInfo;
	vector<int> domainIndex;
	vector<int> pseudoMemIndex;
	vector<int> boundaryIndex;
}GeometryInfo;

typedef struct _variableInfo {
	Species *sp;
	Compartment *com;
	Parameter *para;
	const char* id;
	double *value;
	double *delta;
	bool inVol;
	bool isUniform;
	bool hasAssignmentRule;
	//boundaryType bType;
	reversePolishInfo *rpInfo;
	GeometryInfo *geoi;
	_variableInfo **diffCInfo;
	_variableInfo **adCInfo;
	_variableInfo **boundaryInfo;
	bool isResolved;
	vector<_variableInfo*> dependence;
}variableInfo;

typedef struct _reactionInfo {
	const char* id;
	double *value;
	reversePolishInfo *rpInfo;
	bool isMemTransport;
	Reaction *reaction;
	vector<_variableInfo*> spRefList;
	vector<bool> isVariable;
	vector<double> srStoichiometry;
}reactionInfo;

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

typedef struct _normalUnitVector {
	double nx;
	double ny;
	double nz;
}normalUnitVector;

typedef struct _voronoiInfo {
	bool isAveXY[2];
	bool isAveYZ[2];
	bool isAveXZ[2];
	int adjacentIndexXY[2];
	int adjacentIndexYZ[2];
	int adjacentIndexXZ[2];
	double siXY[2];
	double siYZ[2];
	double siXZ[2];
	double diXY[2];
	double diYZ[2];
	double diXZ[2];
}voronoiInfo;

typedef struct _planeAdjacent {
	normalUnitVector XYcontour[2];
	normalUnitVector YZcontour[2];
	normalUnitVector XZcontour[2];
}planeAdjacent;

#endif /* MYSTRUCT_H_ */
