#include "sbml/SBMLTypes.h"
#include "sbml/extension/SBMLExtensionRegistry.h"
#include "sbml/packages/spatial/common/SpatialExtensionTypes.h"
#include "sbml/packages/spatial/extension/SpatialModelPlugin.h"
#include "sbml/packages/spatial/extension/SpatialExtension.h"
#include <vector>
#include "mystruct.h"

variableInfo* searchInfoById(vector<variableInfo*> &varInfoList, const char *varName)
{
	vector<variableInfo*>::iterator it = varInfoList.begin();
	while (it != varInfoList.end()) {
		if (strcmp((*it)->id, varName) == 0) {
			return *it;
		}
		++it;
	}
	return 0;
}

GeometryInfo* searchAvolInfoByDomainType(vector<GeometryInfo*> &geoInfoList, const char *dtId)
{
	vector<GeometryInfo*>::iterator it = geoInfoList.begin();
	while (it != geoInfoList.end()) {
		if (strcmp((*it)->domainTypeId, dtId) == 0) {
			return *it;
		}
		++it;
	}
	return 0;
}

GeometryInfo* searchAvolInfoByDomain(vector<GeometryInfo*> &geoInfoList, const char *dId)
{
	vector<GeometryInfo*>::iterator it = geoInfoList.begin();
	while (it != geoInfoList.end()) {
		if (strcmp((*it)->domainId, dId) == 0) {
			return *it;
		}
		++it;
	}
	return 0;
}

GeometryInfo* searchAvolInfoByCompartment(vector<GeometryInfo*> &geoInfoList, const char *cId)
{
	vector<GeometryInfo*>::iterator it = geoInfoList.begin();
	while (it != geoInfoList.end()) {
		if (strcmp((*it)->compartmentId, cId) == 0) {
			return *it;
		}
		++it;
	}
	return 0;
}
