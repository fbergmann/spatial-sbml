
#ifndef SPATIAL_SIMULATOR_H
#define SPATIAL_SIMULATOR_H

#include <sbml/common/common.h>
#include <string>
#include <vector>

LIBSBML_CPP_NAMESPACE_BEGIN

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

LIBSBML_CPP_NAMESPACE_END

typedef enum _boundaryType {
  Xp = 0, Xm, Yp, Ym, Zp, Zm
} boundaryType;

#include "spatialstructs.h"


class SpatialSimulator 
{
public: 
  SpatialSimulator(SBMLDocument* doc, int xdim, int ydim, int zdim=1);
  virtual ~SpatialSimulator();

  void setGnuplotExecutable(const std::string &location);
  const std::string& getGnuplotExecutable() const;

  // run the simulator for a bit
  void run(double endTime, double step);

  // perform one step while providing output
  double oneStep(double initialTime, double step);

  void setParameterUniformly(const std::string &id, double value);
  void setParameter(const std::string &id, double value);
  void setParameterUniformly(variableInfo *species, double value);
  int getIndexForPosition(double x, double y);
  // return the values for the given variable (as 1d) flat array
  double* getVariable(const std::string &speciesId, int &length);
  double* getGeometry(const std::string &compartmentId, int &length);
  boundaryType* getBoundaryType(const std::string &compartmentId, int &length);
  int* getBoundary(const std::string &compartmentId, int &length);
  double* getX(int &length);
  double* getY(int &length);
  double* getZ(int &length);

  int getXDim() const { return Xdiv;}
  int getYDim() const { return Ydiv;}
  int getZDim() const { return Zdiv;}

  // for comparison allow the old stuff to run too
  static int runOldMain(int argc, const char* argv[]);

  const Model* getModel() const;

  void deleteValuesOutsideDomain(variableInfo *info );
  void deleteValuesOutsideDomain(const std::string& id);

private:

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

};

#endif //SPATIAL_SIMULATOR_H 