/**
 * \file    spatialsimulator.h
 * \brief   Definition of SpatialSimulator the main entry point for SpatialSBML
 * \author  Frank T. Bergmann <fbergman@caltech.edu>
 * \author  Tatsuhiro Matsui 
 *
 */

#ifndef SPATIAL_SIMULATOR_H
#define SPATIAL_SIMULATOR_H

#include <sbml/common/common.h>
#include <string>
#include <vector>

#ifndef SWIG

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


# endif

#include "spatialstructs.h"


class SpatialSimulator 
{
public: 
  /**
   * Default constructor, constructing an uninitialized simulator, that needs to be 
   * initialized using one of the init functions
   */
  SpatialSimulator();

  /** 
   * Destructor, frees all elements that can be freed
   */
  virtual ~SpatialSimulator();

  /** 
   * Initializes this simulator, by reading the given SBML file and allocating 
   * the needed grid. 
   * 
   * @param fileName the SBML file to load
   * @param xdim the dimensions of the grid in the x direction
   * @param ydim the dimensions of the grid in the y direction
   * @param zdim optional dimensions of the grid in the z direction (defaults to 1)
   *
   */ 
  void initFromFile(const std::string &fileName, int xdim, int ydim, int zdim=1);

  /** 
   * Initializes this simulator, by reading the given SBML content and allocating 
   * the needed grid. 
   * 
   * @param sbml the SBML content to load
   * @param xdim the dimensions of the grid in the x direction
   * @param ydim the dimensions of the grid in the y direction
   * @param zdim optional dimensions of the grid in the z direction (defaults to 1)
   *
   */ 
  void initFromString(const std::string &sbml, int xdim, int ydim, int zdim=1);

  /** 
   * Perform one step with the specified stepsize. 
   * 
   * @param initialTime the current time
   * @param step the stepsize to be used
   * 
   * @return the end time (initialTime + step)
   */
  double oneStep(double initialTime, double step);

  /**
   * Looks through all parameters with given Id, and sets them 
   * uniformly (along the whole grid) to the given value. 
   * 
   * @param id the id of the variable to set
   * @param value the value to set
   */ 
  void setParameterUniformly(const std::string &id, double value);

  /** 
   * Sets the parameter with given id to the specified value. Prior to 
   * calling setParameterUniformly, this function also sets the parameter
   * of the internal SBML model.
   * 
   * @param id the id of the parameter to set
   * @param value the value to set
   */ 
  void setParameter(const std::string &id, double value);

  /**
   * Converts the given coordinates into a internal position 
   * 
   * @param x the x coordinate
   * @param y the y coordinate
   * @return the internal index for that position or -1 if non existent
   */ 
  int getIndexForPosition(double x, double y);

  /** 
   * Returns the internal length for 1d variable arrays. 
   * 
   * @return the internal length for 1d variable arrays. 
   */ 
  int getVariableLength() const;

  /** 
   * Returns the internal length for 1d geometry arrays. 
   * 
   * @return the internal length for 1d geometry arrays. 
   */ 
  int getGeometryLength() const;

  /** 
   * @return the x dimension
   */
  int getXDim() const { return Xdiv;}

  /** 
   * @return the y dimension
   */
  int getYDim() const { return Ydiv;}

  /**
   * @return the z dimension
   */
  int getZDim() const { return Zdiv;}

  /** 
   * Returns the value of the variable at the given position. 
   * 
   * @param id the id of the variable
   * @param x the index along the x coordinate
   * @param y the index along the y coordinate
   * @param z optionally the index along the z coordinate
   * 
   * @return the value of the variable at the specified coordinate
   */
  double getVariableAt(const std::string& id, int x, int y, int z=0);


#ifndef SWIG // API hidded from SWIG

  /** 
   * Constructs and initializes a new instance of the Spatial Simulator. 
   */
  SpatialSimulator(SBMLDocument* doc, int xdim, int ydim, int zdim=1);  

  /** 
   * Initializes this simulator, by reading the given SBML document and allocating 
   * the needed grid. 
   * 
   * @param doc the SBML document to load
   * @param xdim the dimensions of the grid in the x direction
   * @param ydim the dimensions of the grid in the y direction
   * @param zdim optional dimensions of the grid in the z direction (defaults to 1)
   *
   */ 
  void initFromModel(SBMLDocument* doc, int xdim, int ydim, int zdim=1);
  
  /**
   * @returns a pointer to the libSBML model object. 
   */ 
  const Model* getModel() const;

  /**
   * flips the order of analytic volumes
   */ 
  void flipVolumeOrder();
  


  /**
   * Sets the specified value uniformly on the provided variable info
   * structure. 
   * 
   * @param info the variable info element whose value to set
   * @param value the value to set
   */ 
  void setParameterUniformly(variableInfo *info, double value);

  /**
   * Returns the values of the variable with the specified id as (1d) flat array. 
   * to be indexed via getIndexForPosition.
   * 
   * @param id the id of the variable whose values to return
   * @param length output parameter that holds the length of the array
   * 
   * @return the value array, or NULL in case no variable with given id exists
   */
  double* getVariable(const std::string &id, int &length);

  /**
   * Returns the geometry array of the compartment with the specified id as (1d) flat array. 
   * to be indexed via getIndexForPosition.
   * 
   * @param id the id of the compartment whose geometry to return
   * @param length output parameter that holds the length of the array
   * 
   * @return the geometry array, or NULL in case no compartment with given id exists
   */
  int* getGeometry(const std::string &id, int &length);

  /**
   * Returns the boundary type array of the compartment with the specified id as (1d) flat array. 
   * to be indexed via getIndexForPosition.
   * 
   * @param id the id of the compartment whose geometry to return
   * @param length output parameter that holds the length of the array
   * 
   * @return the boundary type array, or NULL in case no compartment with given id exists
   */
  boundaryType* getBoundaryType(const std::string &id, int &length);

  /**
   * Returns the boundary array of the compartment with the specified id as (1d) flat array. 
   * to be indexed via getIndexForPosition.
   * 
   * @param id the id of the compartment whose geometry to return
   * @param length output parameter that holds the length of the array
   * 
   * @return the boundary array, or NULL in case no compartment with given id exists
   */
  int* getBoundary(const std::string &id, int &length);

  /** 
   * Returns the array of internal x coordinates.
   * 
   * @param length output parameter that holds the length of the array
   * 
   * @return the array of internal x coordinates, or NULL in case no initialization happened yet
   */
  double* getX(int &length);

  /** 
   * Returns the array of internal y coordinates.
   * 
   * @param length output parameter that holds the length of the array
   * 
   * @return the array of internal y coordinates, or NULL in case no initialization happened yet
   */
  double* getY(int &length);

  /** 
   * Returns the array of internal z coordinates.
   * 
   * @param length output parameter that holds the length of the array
   * 
   * @return the array of internal z coordinates, or NULL in case no initialization happened yet
   */
  double* getZ(int &length);


 /**
   * Sets the path to the GNUplot executable, that is called when 
   * the run method is called. 
   */
  void setGnuplotExecutable(const std::string &location);

  /** 
   * @return the currently used GNUplot executable. 
   */
  const std::string& getGnuplotExecutable() const;

  /** 
   * Runs the simulator and outputs results using GNUplot
   */
  void run(double endTime, double step);

  /** 
   * Runs the legacy main method.
   */
  static int runOldMain(int argc, const char* argv[]);
  

#endif


private:
  // end of public API
#ifndef SWIG

  void deleteValuesOutsideDomain(variableInfo *info );

  void deleteValuesOutsideDomain(const std::string& id);

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
  vector<GeometryInfo*> geoInfoList;    
  std::vector<reactionInfo*> rInfoList             ;
  std::vector<reactionInfo*> fast_rInfoList;
  std::vector<const char*> memList                 ;
  vector<double*> freeConstList;
  vector<string> spIdList;
  unsigned int dimension                      ;
  normalUnitVector * nuVec;
  voronoiInfo *vorI;
  vector<variableInfo*> orderedARule;

  int Xplus1, Xminus1, Yplus1, Yminus1, Zplus1, Zminus1, divIndex;

  variableInfo *t_info;
  variableInfo *xInfo, *yInfo, *zInfo;
  double deltaX, deltaY, deltaZ;
  char *xaxis, *yaxis, *zaxis;
  Model *model               ;
  Species* s;
  Reaction *r;
  SpeciesReference *sr;

  double Xsize, Ysize, Zsize;
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
