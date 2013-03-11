
#ifndef SPATIAL_SIMULATOR_H
#define SPATIAL_SIMULATOR_H

#include <sbml/common/common.h>
#include <string>
#include <vector>

#ifndef SWIG

LIBSBML_CPP_NAMESPACE_BEGIN
# endif
class SBMLDocument;
class Model;
class Species;
class Parameter;
class Compartment;
class AdjacentDomains;
class BoundaryCondition;
class DiffusionCoefficient;
class AdvectionCoefficient;
class ListOfSpecies;
class ListOfCompartments;
class ListOfReactions;
class ListOfParameters;
class ListOfRules;
#ifndef SWIG
LIBSBML_CPP_NAMESPACE_END


typedef enum _boundaryType {
  Xp = 0, Xm, Yp, Ym, Zp, Zm
} boundaryType;

# endif
#include "spatialstructs.h"


class SpatialSimulator 
{
public: 
  SpatialSimulator();
  void initFromFile(const std::string &fileName, int xdim, int ydim, int zdim=1);
  void initFromString(const std::string &sbml, int xdim, int ydim, int zdim=1);
#ifndef SWIG
  void initFromModel(SBMLDocument* doc, int xdim, int ydim, int zdim=1);
  SpatialSimulator(SBMLDocument* doc, int xdim, int ydim, int zdim=1);  
#endif
  virtual ~SpatialSimulator();

  void setGnuplotExecutable(const std::string &location);
  const std::string& getGnuplotExecutable() const;

  // run the simulator for a bit
  void run(double endTime, double step);

  // perform one step while providing output
  double oneStep(double initialTime, double step);

  void setParameterUniformly(const std::string &id, double value);
  void setParameter(const std::string &id, double value);
  int getIndexForPosition(double x, double y);
#ifndef SWIG
  void setParameterUniformly(variableInfo *species, double value);
#endif
#ifndef SWIG
  // return the values for the given variable (as 1d) flat array
  double* getVariable(const std::string &speciesId, int &length);
  double* getGeometry(const std::string &compartmentId, int &length);

  boundaryType* getBoundaryType(const std::string &compartmentId, int &length);

  int* getBoundary(const std::string &compartmentId, int &length);
  double* getX(int &length);
  double* getY(int &length);
  double* getZ(int &length);
#endif

  int getVariableLength() const;
  int getGeometryLength() const;
  int getXDim() const { return Xdiv;}
  int getYDim() const { return Ydiv;}
  int getZDim() const { return Zdiv;}

  double getVariableAt(const std::string& variable, int x, int y, int z=0);

#ifndef SWIG


  // for comparison allow the old stuff to run too
  static int runOldMain(int argc, const char* argv[]);


  const Model* getModel() const;
  void deleteValuesOutsideDomain(variableInfo *info );

#endif

  void deleteValuesOutsideDomain(const std::string& id);

private:

#ifndef SWIG
  Parameter* getDiffusionCoefficientForSpecies(const std::string& id, int index=0);

  // prints all values from the current point;
  void printValues();

  // perform one RK4 step
  void performStep(double initialTime,double step);

  // update values after RK4 step
  void updateValues(double step);

  void updateAssignmentRules();
  void updateSpecies(variableInfo *species, double step);

  std::string resultDirectory;
  std::string gnuplotExecutable;

  FILE* gp;

  int Xdiv, Ydiv, Zdiv,file_num ;
  double sim_time;

  unsigned int numOfSpecies      ;
  unsigned int numOfReactions    ;
  unsigned int numOfCompartments ;
  unsigned int numOfParameters   ;
  unsigned int numOfRules        ;

  std::vector<variableInfo*> varInfoList           ;
  std::vector<bcOfSpeciesInfo*> bcOfSpeciesInfoList;
  std::vector<analyticVolInfo*> avolInfoList       ;
  std::vector<boundaryCInfo*> bcInfoList           ;
  std::vector<reactionInfo*> rInfoList             ;
  std::vector<adCInfo*> acInfoList                 ;
  std::vector<const char*> memList                 ;
  unsigned int dimension                      ;

  variableInfo *t_info;
  variableInfo *xInfo, *yInfo, *zInfo;
  double deltaX, deltaY, deltaZ;
  char *xaxis, *yaxis, *zaxis;
  Model *model               ;
  unsigned int volDimension, memDimension ;
  
  int Xindex, Yindex , Zindex, numOfVolIndexes;

  ListOfSpecies *los      ;
  ListOfCompartments *loc ;
  ListOfReactions *lor    ;
  ListOfParameters *lop   ;
  ListOfRules *lorules    ;
#endif
};

#endif //SPATIAL_SIMULATOR_H 