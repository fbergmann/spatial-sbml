#ifndef SEARCHFUNCTION_H_
#define SEARCHFUNCTION_H_

variableInfo* searchInfoById(vector<variableInfo*> &varInfoList, const char *varName);
GeometryInfo* searchAvolInfoByDomainType(vector<GeometryInfo*> &geoInfoList, const char *dtId);
GeometryInfo* searchAvolInfoByDomain(vector<GeometryInfo*> &geoInfoList, const char *dId);
GeometryInfo* searchAvolInfoByCompartment(vector<GeometryInfo*> &geoInfoList, const char *cId);

#endif
