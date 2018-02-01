/**
* \file    spatialsimulator.cpp
* \brief   implementation of SpatialSimulator the main entry point for SpatialSBML
* \author  Frank T. Bergmann <fbergman@caltech.edu>
* \author  Tatsuhiro Matsui 
*
*/

#include "spatialsimulator.h"


#define INIT_DOUBLE(destination,size)\
  {\
  memset(destination, 0, sizeof(double)*size);\
  }

#define INIT_DOUBLE_WITH_VALUE(destination,size, value)\
  {\
  memset(destination, value, sizeof(double)*size);\
  }


#define INIT_INT(destination,size)\
  {\
  memset(destination, 0, sizeof(int)*size);\
  }


#include <cstdio>
#include <cstdlib>
#include <cstring>

#if WIN32
#define _USE_MATH_DEFINES
#ifndef LIBSBML_EXPORTS
#define LIBSBML_EXPORTS
#endif
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <sys/stat.h>

#include <sbml/SBMLTypes.h>

#include <sbml/common/extern.h>
#include <sbml/common/common.h>
#include <sbml/common/libsbml-namespace.h>

#include <sbml/extension/SBMLExtensionRegistry.h>

#include <sbml/packages/spatial/common/SpatialExtensionTypes.h>
#include <sbml/packages/spatial/extension/SpatialModelPlugin.h>
#include <sbml/packages/spatial/extension/SpatialExtension.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <stack>
#include <stdlib.h>
#include "spatialstructs.h"

#ifndef WIN32
#define MKDIR(file) (mkdir(file, 0755))
#else
#include <direct.h>
#if _MSC_VER > 1400
#define MKDIR(file) (_mkdir(file))
#else
#define MKDIR(file) (mkdir(file))
#endif
#endif

#ifdef WIN32
#define popen _popen
#define pclose _pclose
#define copysign _copysign
#include <math.h>
#endif

#define stackMax 25

using namespace std;

#include "initializeFunction.h"
#include "freeFunction.h"
#include "searchFunction.h"
#include "astFunction.h"
#include "calcPDE.h"
#include "setInfoFunction.h"
#include "boundaryFunction.h"
#include "outputFunction.h"
#include "checkStability.h"


#include "calcFastReaction.cxx"


std::string getFirstCompartmentId(const Model* model, const ASTNode *ast)
{

  for (unsigned int i = 0;i < ast->getNumChildren();++i)
  {
    string id = getFirstCompartmentId(model, ast->getChild(i));
    if (!id.empty()) return id;
  }

  if (!ast->isName()) return "";

  if (model->getCompartment(ast->getName()) != NULL)
    return ast->getName();
  
  return "";

}

bool isResolvedAll(vector<variableInfo*> &dependence)
{
  if (dependence.size() == 0) return true;
  vector<variableInfo*>::iterator it = dependence.begin();
  while (it != dependence.end()) {
    if ((*it) == NULL || !(*it)->isResolved) {
      return false;
    }
    ++it;
  }
  return true;
}


bool isNormalParameter(SpatialParameterPlugin *pPlugin ) 
{
  if (pPlugin == NULL) return true;

  if (pPlugin ->getSpatialSymbolReference() != NULL && !pPlugin ->getSpatialSymbolReference()->getId().empty())
    return false;

  if (pPlugin ->getDiffusionCoefficient() != NULL && !pPlugin ->getDiffusionCoefficient()->getVariable().empty())
    return false;

  if (pPlugin ->getBoundaryCondition() != NULL && !pPlugin ->getBoundaryCondition()->getVariable().empty())
    return false;

  if (pPlugin ->getAdvectionCoefficient() != NULL && !pPlugin ->getAdvectionCoefficient()->getVariable().empty())
    return false;

  return true;
}

int SpatialSimulator::getVariableLength() const
{
  return (Zindex -1) * Yindex * Xindex + (Yindex-1) * Xindex + Xindex-1;
}

int SpatialSimulator::runOldMain(int argc, const char* argv[])
{
  return -1;//old_main(argc, argv);
}

void SpatialSimulator::setParameterUniformly(const std::string& id, double value)
{
  setParameterUniformly(searchInfoById(varInfoList, id.c_str()), value);
}

void SpatialSimulator::setParameter(const std::string &id, double value)
{
  Parameter* param = model->getParameter(id);  
  if (param == NULL) return;
  param->setValue(value);
  SpatialParameterPlugin *pPlugin = static_cast<SpatialParameterPlugin*>(param->getPlugin("spatial"));
  if (pPlugin != NULL && !isNormalParameter(pPlugin))
  {

    DiffusionCoefficient *dc = pPlugin->getDiffusionCoefficient();
    if (dc != NULL) {
      variableInfo *sInfo = searchInfoById(varInfoList,dc->getVariable().c_str());
      if (sInfo != 0) {
        variableInfo *dci = sInfo->diffCInfo[dc->getCoordinateReference1()];
        if (dci != NULL)
        {
          if (dci->isUniform && dci->value != NULL)
          {
            *(dci->value) = value;
          }
          else if (dci->value != NULL)
          {
            for (int Z = 0; Z < Zindex; Z++) {
              for (int Y = 0; Y < Yindex; Y++) {
                for (int X = 0; X < Xindex; X++) {
                  dci->value[Z * Yindex * Xindex + Y * Xindex + X] = value;
                }
              }
            }
        }
        }
      }
    }
  }
  setParameterUniformly(id, value);
}

void SpatialSimulator::setParameterUniformly(variableInfo* info, double value)
{
  if (info == NULL) return;

  if (info->isUniform && info->value != NULL)
  {
    *(info->value) = value;
  }
  else if(info->value != NULL)
  {
  for (int Z = 0; Z < Zindex; Z++) {
#pragma omp parallel for
    for (int Y = 0; Y < Yindex; Y++) {
      for (int X = 0; X < Xindex; X++) {
        info->value[Z * Yindex * Xindex + Y * Xindex + X] = value;
      }
    }
  }
  }
}


Parameter* SpatialSimulator::getDiffusionCoefficientForSpecies(const std::string& id, int index)
{
  for(unsigned int i = 0; i < model->getNumParameters(); i++)
  {
    Parameter* current = model->getParameter(i);
    SpatialParameterPlugin *pPlugin = static_cast<SpatialParameterPlugin*>(current->getPlugin("spatial"));
    if (pPlugin == NULL) continue;
    DiffusionCoefficient* coefficient = pPlugin->getDiffusionCoefficient();
    if (coefficient == NULL || coefficient->getVariable() != id || coefficient->getCoordinateReference1() != index)
      continue;    
    return current;
  }
  return NULL;
}

SpatialSimulator::SpatialSimulator():
  Xdiv(100), Ydiv(100), Zdiv(1),  deltaX(0), deltaY(0), deltaZ(0), model(NULL),volDimension(0), memDimension(0)
{
}
void SpatialSimulator::initFromFile(const std::string &fileName, int xdim, int ydim, int zdim/*=1*/)
{
  SBMLDocument*doc = readSBMLFromFile(fileName.c_str());
  initFromModel(doc, xdim, ydim, zdim);
}
void SpatialSimulator::initFromString(const std::string &sbml, int xdim, int ydim, int zdim/*=1*/)
{
  SBMLDocument*doc = readSBMLFromString(sbml.c_str());
  initFromModel(doc, xdim, ydim, zdim);
}


SpatialSimulator::SpatialSimulator(SBMLDocument* doc, int xdim, int ydim, int zdim /*=1*/) :
  Xdiv(xdim), Ydiv(ydim), Zdiv(zdim), model(NULL), volDimension(0), memDimension(0)
{
  initFromModel(doc, xdim, ydim, zdim);
}


AnalyticVolume* getAnalyticVolumeForType(AnalyticGeometry* analyticGeo,const std::string& domainType )
{
  for (size_t i =0; i < analyticGeo->getNumAnalyticVolumes(); ++i)
  {
    AnalyticVolume* current = analyticGeo->getAnalyticVolume(i);
    if (current->getDomainType() == domainType)
      return current;
  }
  return NULL;
}

void SpatialSimulator::initFromModel(SBMLDocument* doc, int xdim, int ydim, int zdim/*=1*/)
{
  unsigned int i, j, k;
  int X = 0, Y = 0, Z = 0, index = 0;	  
  int numOfASTNodes = 0;

  Xdiv = xdim;
  Ydiv = ydim;
  Zdiv = zdim;

  //sbml core
  ASTNode *ast;
  model               = doc->getModel();

  los         = model->getListOfSpecies();
  loc    = model->getListOfCompartments();
  lor       = model->getListOfReactions();
  lop      = model->getListOfParameters();
  lorules       = model->getListOfRules();

  //sbml spatial package
  SpatialModelPlugin *spPlugin     = static_cast<SpatialModelPlugin*>(model->getPlugin("spatial"));
  Geometry *geometry               = spPlugin->getGeometry();
  //ListOfCoordinateComponents *locc = geometry->getListOfCoordinateComponents();

  //size of list
  numOfSpecies      = static_cast<unsigned int>(model->getNumSpecies());
  numOfReactions    = static_cast<unsigned int>(model->getNumReactions());
  numOfCompartments = static_cast<unsigned int>(model->getNumCompartments());
  numOfParameters   = static_cast<unsigned int>(model->getNumParameters());
  numOfRules        = static_cast<unsigned int>(model->getNumRules());

  varInfoList            = vector<variableInfo*>();
  varInfoList.reserve(numOfCompartments + numOfSpecies + numOfParameters);
  geoInfoList = vector<GeometryInfo*>();
  geoInfoList.reserve(geometry->getNumGeometryDefinitions());

  rInfoList              = vector<reactionInfo*>();
  rInfoList.reserve(numOfReactions + numOfRules);
  fast_rInfoList = vector<reactionInfo*>();
  fast_rInfoList.reserve(numOfReactions + numOfRules);
  memList = vector<const char*>();
  memList.reserve(numOfCompartments);
  dimension                       = geometry->getNumCoordinateComponents();

  spIdList = vector<string>();
  spIdList.reserve(numOfSpecies);
  for (i = 0; i < numOfSpecies; i++) {
    spIdList.push_back(model->getSpecies(i)->getId());
  }

  Xplus1 = 0, Xminus1 = 0, Yplus1 = 0, Yminus1 = 0, Zplus1 = 0, Zminus1 = 0;

  //div
  if (dimension <= 1) {
    Ydiv = 1;
    Zdiv = 1;
  }
  if (dimension <= 2) {
    Zdiv = 1;
  }




  /*bool isImageBased = false;
  for (i = 0; i < geometry->getNumGeometryDefinitions(); i++) {
    if (geometry->getGeometryDefinition(i)->isSampledFieldGeometry()) {
      //SampleFieldGeometry
      SampledFieldGeometry *sfGeo	= static_cast<SampledFieldGeometry*>(geometry->getGeometryDefinition(i));
      SampledField *samField = sfGeo->getSampledField();
      isImageBased = true;
      Xdiv = samField->getNumSamples1();
      Ydiv = samField->getNumSamples2();
      Zdiv = samField->getNumSamples3();
    }
  }*/


  Xindex = 2 * Xdiv - 1;
  Yindex = 2 * Ydiv - 1;
  Zindex = 2 * Zdiv - 1;

  numOfVolIndexes = Xindex * Yindex * Zindex;

  volDimension = 0; memDimension = 0;
  for (i = 0; i < numOfCompartments; i++) {
    Compartment *c = loc->get(i);
    if (volDimension == 0) {
      volDimension = c->getSpatialDimensions();
    } else {
      if (c->getSpatialDimensions() >= volDimension) {
        volDimension = c->getSpatialDimensions();
      } else {
        memDimension = c->getSpatialDimensions();
      }
    }
  }
  Xsize = 0.0;
  Ysize = 0.0;
  Zsize = 0.0;

  //set id and value
  //compartment
  setCompartmentInfo(model, varInfoList);
  //species
  setSpeciesInfo(doc, varInfoList, volDimension, memDimension, Xindex, Yindex, Zindex);
  //parameter
  setParameterInfo(doc, varInfoList, Xdiv, Ydiv, Zdiv, Xsize, Ysize, Zsize, deltaX, deltaY, deltaZ, xaxis, yaxis, zaxis);

  //time
  variableInfo *t_info = new variableInfo;
  InitializeVarInfo(t_info);
  varInfoList.push_back(t_info);
  t_info->id = (const char*)malloc(sizeof(char) * 1 + 1);
  strcpy(const_cast<char*>(t_info->id), "t");
  t_info->value = &sim_time;
  t_info->isResolved = true;
  t_info->isUniform = true;

  //volume index
  vector<int> volumeIndexList;
  for (Z = 0; Z < Zindex; Z += 2) 
  {
    for (Y = 0; Y < Yindex; Y += 2) 
    {
      for (X = 0; X < Xindex; X += 2) 
      {
        volumeIndexList.push_back(Z * Yindex * Xindex + Y * Xindex + X);
      }
    }
  }

  SpatialCompartmentPlugin* cPlugin = NULL;

  //geometryDefinition
  double *tmp_isDomain = new double[numOfVolIndexes];
  for (i = 0; i < geometry->getNumGeometryDefinitions(); i++) 
  {
    if (geometry->getGeometryDefinition(i)->isAnalyticGeometry()) 
    {
      //AnalyticVolumes
      AnalyticGeometry *analyticGeo = static_cast<AnalyticGeometry*>(geometry->getGeometryDefinition(i));
      //gather information of compartment, domainType, analyticVolume
      for (j = 0; j < numOfCompartments; j++) 
      {
        Compartment *c = loc->get(j);
        cPlugin = static_cast<SpatialCompartmentPlugin*>(c->getPlugin("spatial"));
        if (cPlugin != 0) {
          //          					if (AnalyticVolume *analyticVol = analyticGeo->getAnalyticVolume(cPlugin->getCompartmentMapping()->getDomainType())) {
          //
          AnalyticVolume *analyticVol =  getAnalyticVolumeForType( analyticGeo,cPlugin->getCompartmentMapping()->getDomainType() );
          if (analyticVol != NULL) {
            //  analyticGeo->getAnalyticVolume(cPlugin->getCompartmentMapping()->getDomainType())) {
            GeometryInfo *geoInfo = new GeometryInfo;
            InitializeAVolInfo(geoInfo);
            geoInfo->compartmentId = c->getId().c_str();
            geoInfo->domainTypeId = cPlugin->getCompartmentMapping()->getDomainType().c_str();
            geoInfo->domainId = 0;
            geoInfo->bType = new boundaryType[numOfVolIndexes];
            for (k = 0; k < (unsigned int)numOfVolIndexes; k++) {
              geoInfo->bType[k].isBofXp = false;
              geoInfo->bType[k].isBofXm = false;
              geoInfo->bType[k].isBofYp = false;
              geoInfo->bType[k].isBofYm = false;
              geoInfo->bType[k].isBofZp = false;
              geoInfo->bType[k].isBofZm = false;
            }
            geoInfo->isVol = true;
            ast = const_cast<ASTNode*>(analyticVol->getMath());
            if (ast == NULL) continue;
            rearrangeAST(ast);
            numOfASTNodes = 0;
            countAST(ast, numOfASTNodes);
            geoInfo->rpInfo = new reversePolishInfo();
            geoInfo->rpInfo->varList = new double*[numOfASTNodes];
            memset(geoInfo->rpInfo->varList, 0, numOfASTNodes*sizeof(double*));
            geoInfo->rpInfo->constList = new double*[numOfASTNodes];
            memset(geoInfo->rpInfo->constList, 0, numOfASTNodes*sizeof(double*));
            geoInfo->rpInfo->opfuncList = new int[numOfASTNodes];
            memset(geoInfo->rpInfo->opfuncList, 0, numOfASTNodes*sizeof(int));
            geoInfo->rpInfo->listNum = numOfASTNodes;
            geoInfo->isDomain = new int[numOfVolIndexes];
            memset(geoInfo->isDomain, 0, numOfASTNodes*sizeof(int));
            geoInfo->isBoundary = new int[numOfVolIndexes];
            memset(geoInfo->isBoundary, 0, numOfASTNodes*sizeof(int));
            geoInfo->adjacent0 = 0;
            geoInfo->adjacent1 = 0;
            parseAST(ast, geoInfo->rpInfo, varInfoList, numOfASTNodes, freeConstList);
            //judge if the coordinate point is inside the analytic volume
            memset(tmp_isDomain, 0, numOfVolIndexes*sizeof(double));
            reversePolishInitial(volumeIndexList, geoInfo->rpInfo, tmp_isDomain, numOfASTNodes, Xindex, Yindex, Zindex, false);
            for (k = 0; k < (unsigned int)numOfVolIndexes; k++) {
              //index = k;
              geoInfo->isDomain[k] = (int)tmp_isDomain[k];
              //Z = index / (Xindex * Yindex);
              //Y = (index - Z * Xindex * Yindex) / Xindex;
              //X = index - Z * Xindex * Yindex - Y * Xindex;
              //if (geoInfo->isDomain[k] == 1 && string(geoInfo->domainTypeId) == "nucleus") {
              //	//cerr << X << ", " << Y << ", " << Z << ", " << geoInfo->isDomain[k] << endl;
              //} else {
              //	//cout << geoInfo->domainTypeId << endl;
              //}
            }
            geoInfoList.push_back(geoInfo);
          }
        }
      }
    } else if (geometry->getGeometryDefinition(i)->isSampledFieldGeometry()) {
      //SampleFieldGeometry
      SampledFieldGeometry *sfGeo	= static_cast<SampledFieldGeometry*>(geometry->getGeometryDefinition(i));
      for (j = 0; j < numOfCompartments; j++) {
        Compartment *c = loc->get(j);
        if (c->getSpatialDimensions() == volDimension) 
        {
          cPlugin = static_cast<SpatialCompartmentPlugin*>(c->getPlugin("spatial"));
          if (cPlugin != 0) 
          {
            const std::string& domainType = cPlugin->getCompartmentMapping()->getDomainType();
            SampledField *samField = geometry->getSampledField(sfGeo->getSampledField());
            SampledVolume *samVol = 0;
            for (k = 0; k < sfGeo->getNumSampledVolumes(); k++) {
              if (sfGeo->getSampledVolume(k)->getDomainType() == domainType) {
                samVol = sfGeo->getSampledVolume(k);
              }
            }
            if (samVol == NULL) continue;
            unsigned int  uncomprLen = samField->getUncompressedLength();
            int *uncompr = (int*)calloc(sizeof(int), uncomprLen);
            samField->getUncompressed(uncompr);

            GeometryInfo *geoInfo = new GeometryInfo;
            InitializeAVolInfo(geoInfo);
            geoInfo->compartmentId = c->getId().c_str();
            geoInfo->domainTypeId = domainType.c_str();
            geoInfo->domainId = 0;
            geoInfo->bType = new boundaryType[numOfVolIndexes];
            for (k = 0; k < (unsigned int)numOfVolIndexes; k++) {
              geoInfo->bType[k].isBofXp = false;
              geoInfo->bType[k].isBofXm = false;
              geoInfo->bType[k].isBofYp = false;
              geoInfo->bType[k].isBofYm = false;
              geoInfo->bType[k].isBofZp = false;
              geoInfo->bType[k].isBofZm = false;
            }
            geoInfo->isVol = true;
            geoInfo->isDomain = new int[numOfVolIndexes];
            memset(geoInfo->isDomain, 0, numOfVolIndexes*sizeof(int));
            geoInfo->isBoundary = new int[numOfVolIndexes];
            memset(geoInfo->isBoundary, 0, numOfVolIndexes*sizeof(int));
            geoInfoList.push_back(geoInfo);
            /*
            for (Y = 0; Y <= samField->getNumSamples2(); Y++) {
            for (X = 0; X <= samField->getNumSamples1(); X++) {
            sam_ofs << (unsigned int)uncompr[Y * samField->getNumSamples1() + X] << ", ";
            }
            sam_ofs << endl;
            }
            sam_ofs.close();
            */
            for (Z = 0; Z < Zindex; Z += 2) {
              for (Y = 0; Y < Yindex; Y += 2) {
                for (X = 0; X < Xindex; X += 2) {
                  int currentZ = Z / 2 * samField->getNumSamples2() * samField->getNumSamples1();
                  int currentY = (samField->getNumSamples2() - Y / 2 - 1) * samField->getNumSamples1();
                  int currentX = X / 2;
                  if (static_cast<unsigned int>(uncompr[currentZ + currentY + currentX]) == static_cast<unsigned int>(samVol->getSampledValue())) {
                    geoInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                    geoInfo->domainIndex.push_back(Z * Yindex * Xindex + Y * Xindex + X);
                  } else {
                    geoInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 0;
                  }
                  //sam_ofs << (unsigned int)uncompr[X * samField->getNumSamples1() + Y] << ", ";
                }
              }
            }
            for (k = 0; k < geoInfo->domainIndex.size(); k++) {
              index = geoInfo->domainIndex[k];
              Z = index / (Xindex * Yindex);
              Y = (index - Z * Xindex * Yindex) / Xindex;
              X = index - Z * Xindex * Yindex - Y * Xindex;
              if ((dimension == 2 && (X == 0 || X == Xindex - 1 || Y == 0 || Y == Yindex - 1)) ||
                (dimension == 3 && (X == 0 || X == Xindex - 1 || Y == 0 || Y == Yindex - 1 || Z == 0 || Z == Zindex - 1))) {
                  geoInfo->isBoundary[index] = 1;
                  geoInfo->boundaryIndex.push_back(index);
                  if (dimension >= 2) {
                    if (X == 0) geoInfo->bType[index].isBofXm = true;
                    if (X == Xindex - 1) geoInfo->bType[index].isBofXp = true;
                    if (Y == 0) geoInfo->bType[index].isBofYm = true;
                    if (Y == Yindex - 1) geoInfo->bType[index].isBofYp = true;
                  }
                  if (dimension == 3) {
                    if (Z == 0) geoInfo->bType[index].isBofZm = true;
                    if (Z == Zindex - 1) geoInfo->bType[index].isBofZp = true;
                  }
              } else {
                if (dimension >= 2) {
                  if (geoInfo->isDomain[Z * Xindex * Yindex + Y * Xindex + (X + 2)] == 0) {
                    geoInfo->isBoundary[index] = 1;
                    geoInfo->bType[index].isBofXp = true;
                  }
                  if (geoInfo->isDomain[Z * Xindex * Yindex + Y * Xindex + (X - 2)] == 0) {
                    geoInfo->isBoundary[index] = 1;
                    geoInfo->bType[index].isBofXm = true;
                  }
                  if (geoInfo->isDomain[Z * Xindex * Yindex + (Y + 2) * Xindex + X] == 0) {
                    geoInfo->isBoundary[index] = 1;
                    geoInfo->bType[index].isBofYp = true;
                  }
                  if (geoInfo->isDomain[Z * Xindex * Yindex + (Y - 2) * Xindex + X] == 0) {
                    geoInfo->isBoundary[index] = 1;
                    geoInfo->bType[index].isBofYm = true;
                  }
                }
                if (dimension == 3) {
                  if (geoInfo->isDomain[(Z + 2) * Xindex * Yindex + Y * Xindex + X] == 0) {
                    geoInfo->isBoundary[index] = 1;
                    geoInfo->bType[index].isBofZp = true;
                  }
                  if (geoInfo->isDomain[(Z - 2) * Xindex * Yindex + Y * Xindex + X] == 0) {
                    geoInfo->isBoundary[index] = 1;
                    geoInfo->bType[index].isBofZm = true;
                  }
                }
              }
            }
            free(uncompr);
          }
        }
      }
    } else if (geometry->getGeometryDefinition(i)->isCSGeometry()){
      //CSGeometry
    } else if (geometry->getGeometryDefinition(i)->isParametricGeometry()) {
      //ParametricGeometry
    }
  }
  delete[] tmp_isDomain;
  //merge external and internal analytic volumes and get boundary points of the geometry
  for (i = 0; i < geometry->getNumGeometryDefinitions(); i++) {
    AnalyticGeometry *analyticGeo = static_cast<AnalyticGeometry*>(geometry->getGeometryDefinition(i));
    if (geometry->getGeometryDefinition(i)->isAnalyticGeometry()) {
      for (j = 0; j < analyticGeo->getNumAnalyticVolumes(); j++) {
        AnalyticVolume *analyticVolEx = analyticGeo->getAnalyticVolume(j);
        GeometryInfo *geoInfoEx = searchAvolInfoByDomainType(geoInfoList, analyticVolEx->getDomainType().c_str());
        if (geoInfoEx == NULL)
          continue;

        for (k = 0; k < analyticGeo->getNumAnalyticVolumes(); k++) {
          AnalyticVolume *analyticVolIn = analyticGeo->getAnalyticVolume(k);
          if (analyticVolEx->getOrdinal() < analyticVolIn->getOrdinal()) {//ex's ordinal is smaller than in's ordinal
            GeometryInfo *geoInfoIn = searchAvolInfoByDomainType(geoInfoList, analyticVolIn->getDomainType().c_str());
            if (geoInfoIn == NULL)
              continue;

            //merge
            for (Z = 0; Z < Zindex; Z += 2) {
              for (Y = 0; Y < Yindex; Y += 2) {
                for (X = 0; X < Xindex; X += 2) {
                  index = Z * Yindex * Xindex + Y * Xindex + X;
                  //if external == true && internal == true, then external = false;
                  if (geoInfoEx->isDomain[index] == 1 && geoInfoIn->isDomain[index] == 1) {
                    geoInfoEx->isDomain[index] = 0;
                  }
                }
              }
            }
          }
        }
        //domainの位置(analytic geo)
        for (Z = 0; Z < Zindex; Z += 2) {
          for (Y = 0; Y < Yindex; Y += 2) {
            for (X = 0; X < Xindex; X += 2) {
              index = Z * Xindex * Yindex + Y * Xindex + X;
              if (geoInfoEx->isDomain[index] == 1) {
                geoInfoEx->domainIndex.push_back(index);
              }
            }
          }
        }
        //boundary
        switch(dimension) {
        case 1://1D
          for (X = 0; X < Xindex; X += 2) {
            geoInfoEx->isBoundary[X] = 0;//initialize
            if (geoInfoEx->isDomain[X] == 1) {
              if (X == 0 || X == Xindex - 1) {//bounary of the domain
                geoInfoEx->isBoundary[X] = 1;
                geoInfoEx->boundaryIndex.push_back(X);
              } else if (geoInfoEx->isDomain[X + 2] == 0 || geoInfoEx->isDomain[X - 2] == 0) {
                geoInfoEx->isBoundary[Y * Xindex + X] = 1;
                geoInfoEx->boundaryIndex.push_back(X);
              }
            }
          }
          break;
        case 2://2D
          for (Y = 0; Y < Yindex; Y += 2) {
            for (X = 0; X < Xindex; X += 2) {
              geoInfoEx->isBoundary[Y * Xindex + X] = 0;//initialize
              if (geoInfoEx->isDomain[Y * Xindex + X] == 1) {
                if (X == 0 || X == Xindex - 1 || Y == 0 || Y == Yindex - 1) {//boundary of the domain
                  geoInfoEx->isBoundary[Y * Xindex + X] = 1;
                  geoInfoEx->boundaryIndex.push_back(Y * Xindex + X);
                } else if (geoInfoEx->isDomain[Y * Xindex + (X + 2)] == 0
                  || geoInfoEx->isDomain[Y * Xindex + (X - 2)] == 0
                  || geoInfoEx->isDomain[(Y + 2) * Xindex + X] == 0
                  || geoInfoEx->isDomain[(Y - 2) * Xindex + X] == 0) {
                    geoInfoEx->isBoundary[Y * Xindex + X] = 1;
                    geoInfoEx->boundaryIndex.push_back(Y * Xindex + X);
                }
              }
            }
          }
          break;
        case 3://3D
          for (Z = 0; Z < Zindex; Z += 2) {
            for (Y = 0; Y < Yindex; Y += 2) {
              for (X = 0; X < Xindex; X += 2) {
                geoInfoEx->isBoundary[Z * Yindex * Xindex + Y * Xindex + X] = 0;//initialize
                if (geoInfoEx->isDomain[Z * Yindex * Xindex + Y * Xindex + X] == 1) {
                  if (X == 0 || X == Xindex - 1 || Y == 0 || Y == Yindex - 1 || Z == 0 || Z == Zindex - 1) {//boundary of th domain
                    geoInfoEx->isBoundary[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                    geoInfoEx->boundaryIndex.push_back(Z * Yindex * Xindex + Y * Xindex + X);
                  } else if (geoInfoEx->isDomain[Z * Yindex * Xindex + Y * Xindex + (X + 2)] == 0
                    || geoInfoEx->isDomain[Z * Yindex * Xindex + Y * Xindex + (X - 2)] == 0
                    || geoInfoEx->isDomain[Z * Yindex * Xindex + (Y + 2) * Xindex + X] == 0
                    || geoInfoEx->isDomain[Z * Yindex * Xindex + (Y - 2) * Xindex + X] == 0
                    || geoInfoEx->isDomain[(Z + 2) * Yindex * Xindex + Y * Xindex + X] == 0
                    || geoInfoEx->isDomain[(Z - 2) * Yindex * Xindex + Y * Xindex + X] == 0) {
                      geoInfoEx->isBoundary[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      geoInfoEx->boundaryIndex.push_back(Z * Yindex * Xindex + Y * Xindex + X);
                  }
                }
              }
            }
          }
          break;
        default:
          break;
        }
      }
    }
  }

  //membrane compartment
  //adjacentとかdomainとか使って隣接関係を割り出す。
  for (i = 0; i < numOfCompartments; i++) {
    SpatialCompartmentPlugin *cPlugin = static_cast<SpatialCompartmentPlugin*>(loc->get(i)->getPlugin("spatial"));
    if (cPlugin == 0) continue;
    DomainType *dType = geometry->getDomainType(cPlugin->getCompartmentMapping()->getDomainType());
    GeometryInfo *geoInfo = searchAvolInfoByCompartment(geoInfoList, loc->get(i)->getId().c_str());
    if (geoInfo == NULL) {
      geoInfo = new GeometryInfo;
      InitializeAVolInfo(geoInfo);
      geoInfo->compartmentId = loc->get(i)->getId().c_str();
      geoInfo->domainTypeId = cPlugin->getCompartmentMapping()->getDomainType().c_str();
      //geoInfo->bt = new boundaryType[numOfVolIndexes];
      if (dType == NULL) {
        geoInfo->isVol = true;
      }
      else if (dType->getSpatialDimensions() == volDimension) {
        geoInfo->isVol = true;
      } else if (dType->getSpatialDimensions() == memDimension) {
        geoInfo->isVol = false;
      }
      geoInfoList.push_back(geoInfo);
    }
    for (j = 0; j < geometry->getNumDomains(); j++) {//Domain
      Domain *domain = geometry->getDomain(j);
      if (domain->getDomainType() == geoInfo->domainTypeId) {
        geoInfo->domainIdList.push_back(domain->getId().c_str());
      }
    }
  }

  for (i = 0; i < geometry->getNumAdjacentDomains(); i++) {//AdjacentDomain
    AdjacentDomains *adDomain = geometry->getAdjacentDomains(i);
    Domain *ad1 = geometry->getDomain(adDomain->getDomain1());
    Domain *ad2 = geometry->getDomain(adDomain->getDomain2());

    if (ad1 == NULL || ad2 == NULL)
      continue;

    GeometryInfo *geoi1 = searchAvolInfoByDomainType(geoInfoList, ad1->getDomainType().c_str());
    GeometryInfo *geoi2 = searchAvolInfoByDomainType(geoInfoList, ad2->getDomainType().c_str());
    if (geoi1 == NULL || geoi2 == NULL)
      continue;

    if (geoi1->adjacentGeo1 == 0) {
      geoi1->adjacentGeo1 = geoi2;
    } else if (geoi1->adjacentGeo2 == 0 && geoi1->adjacentGeo1 != geoi2) {
      geoi1->adjacentGeo2 = geoi2;
    }
    if (geoi2->adjacentGeo1 == 0) {
      geoi2->adjacentGeo1 = geoi1;
    } else if (geoi2->adjacentGeo2 == 0 && geoi2->adjacentGeo1 != geoi1) {
      geoi2->adjacentGeo2 = geoi1;
    }
  }

  //membrane position
  for (i = 0; i < geoInfoList.size(); i++) {
    GeometryInfo *geoInfo = geoInfoList[i];
    if (geoInfo->isVol == false) {//avol is membrane
      geoInfo->isDomain = new int[numOfVolIndexes];
      memset(geoInfo->isDomain, 0, numOfVolIndexes*sizeof(int));
      geoInfo->isBoundary = new int[numOfVolIndexes];
      memset(geoInfo->isBoundary, 0, numOfVolIndexes*sizeof(int));
      geoInfo->bType = new boundaryType[numOfVolIndexes];
      for (j = 0; j < (unsigned int)numOfVolIndexes; j++) {
        geoInfo->bType[j].isBofXp = false;
        geoInfo->bType[j].isBofXm = false;
        geoInfo->bType[j].isBofYp = false;
        geoInfo->bType[j].isBofYm = false;
        geoInfo->bType[j].isBofZp = false;
        geoInfo->bType[j].isBofZm = false;
      }

      if (geoInfo->adjacentGeo1 == NULL || geoInfo->adjacentGeo2 == NULL)
        continue;

      switch (dimension) {
      case 1:
        for (X = 0; X < Xindex; X++) {
          if (X % 2 != 0) {
            Xplus1 = Y * Xindex + (X + 1);
            Xminus1 = Y * Xindex + (X - 1);
            if (X != 0 && X != (Xindex - 1)) {
              if ((geoInfo->adjacentGeo1->isDomain[Xplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Xminus1] == 1) ||
                (geoInfo->adjacentGeo1->isDomain[Xminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Xplus1] == 1)) {
                  geoInfo->isDomain[X] = 1;
                  geoInfo->domainIndex.push_back(X);
              }
            }
          }
        }
        break;
      case 2:
        for (Y = 0; Y < Yindex; Y++) {
          for (X = 0; X < Xindex; X++) {
            if ((Y * Xindex + X) % 2 != 0) {
              Xplus1 = Y * Xindex + (X + 1);
              Xminus1 = Y * Xindex + (X - 1);
              Yplus1 = (Y + 1) * Xindex + X;
              Yminus1 = (Y - 1) * Xindex + X;
              if (X != 0 && X != (Xindex - 1) && Y != 0 && Y != (Yindex - 1)) {
                if ((geoInfo->adjacentGeo1->isDomain[Xplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Xminus1] == 1) ||
                  (geoInfo->adjacentGeo1->isDomain[Xminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Xplus1] == 1)) {
                    geoInfo->isDomain[Y * Xindex + X] = 1;
                    geoInfo->domainIndex.push_back(Y * Xindex + X);
                    geoInfo->bType[Y * Xindex + X].isBofXp = true;
                    geoInfo->bType[Y * Xindex + X].isBofXm = true;
                } else if ((geoInfo->adjacentGeo1->isDomain[Yplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Yminus1] == 1) ||
                  (geoInfo->adjacentGeo1->isDomain[Yminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Yplus1] == 1)) {
                    geoInfo->isDomain[Y * Xindex + X] = 1;
                    geoInfo->domainIndex.push_back(Y * Xindex + X);
                    geoInfo->bType[Y * Xindex + X].isBofYp = true;
                    geoInfo->bType[Y * Xindex + X].isBofYm = true;
                }
              } else if (X == 0 || X == Xindex - 1) {
                if ((geoInfo->adjacentGeo1->isDomain[Yplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Yminus1] == 1) ||
                  (geoInfo->adjacentGeo1->isDomain[Yminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Yplus1] == 1)) {
                    geoInfo->isDomain[Y * Xindex + X] = 1;
                    geoInfo->domainIndex.push_back(Y * Xindex + X);
                    geoInfo->bType[Y * Xindex + X].isBofYp = true;
                    geoInfo->bType[Y * Xindex + X].isBofYm = true;
                }
              } else if (Y == 0 || Y == Yindex - 1) {
                if ((geoInfo->adjacentGeo1->isDomain[Xplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Xminus1] == 1) ||
                  (geoInfo->adjacentGeo1->isDomain[Xminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Xplus1] == 1)) {
                    geoInfo->isDomain[Y * Xindex + X] = 1;
                    geoInfo->domainIndex.push_back(X);
                    geoInfo->bType[Y * Xindex + X].isBofXp = true;
                    geoInfo->bType[Y * Xindex + X].isBofXm = true;
                }
              }
            }
          }
        }
        //pseudo membrane
        for (Y = 0; Y < Yindex; Y++) {
          for (X = 0; X < Xindex; X++) {
            //if (!(X % 2 == 0 && Y % 2 == 0)) {
            if ((X % 2 != 0 && Y % 2 != 0)) {
              Xplus1 = Y * Xindex + (X + 1);
              Xminus1 = Y * Xindex + (X - 1);
              Yplus1 = (Y + 1) * Xindex + X;
              Yminus1 = (Y - 1) * Xindex + X;
              if (X != 0 && X != Xindex - 1 && Y != 0 && Y != Yindex - 1){
                if ((geoInfo->isDomain[Xplus1] == 1 && geoInfo->isDomain[Xminus1] == 1) ||
                  (geoInfo->isDomain[Yplus1] == 1 && geoInfo->isDomain[Yminus1] == 1)) {
                    geoInfo->isDomain[Y * Xindex + X] = 2;
                    geoInfo->pseudoMemIndex.push_back(Y * Xindex + X);
                }
                if ((geoInfo->isDomain[Xplus1] == 1 && geoInfo->isDomain[Yplus1] == 1) ||
                  (geoInfo->isDomain[Xplus1] == 1 && geoInfo->isDomain[Yminus1] == 1) ||
                  (geoInfo->isDomain[Yplus1] == 1 && geoInfo->isDomain[Xminus1] == 1) ||
                  (geoInfo->isDomain[Yminus1] == 1 && geoInfo->isDomain[Xminus1] == 1)) {
                    geoInfo->isDomain[Y * Xindex + X] = 2;
                    geoInfo->pseudoMemIndex.push_back(Y * Xindex + X);
                }
              } else if (X == 0 || X == Xindex - 1) {
                if (geoInfo->isDomain[Yplus1] == 1 && geoInfo->isDomain[Yminus1] == 1) {
                  geoInfo->isDomain[Y * Xindex + X] = 2;
                  geoInfo->pseudoMemIndex.push_back(Y * Xindex + X);
                }
              } else if (Y == 0 || Y == Yindex - 1) {
                if (geoInfo->isDomain[Xplus1] == 1 && geoInfo->isDomain[Xminus1] == 1) {
                  geoInfo->isDomain[Y * Xindex + X] = 2;
                  geoInfo->pseudoMemIndex.push_back(Y * Xindex + X);
                }
              }
            }
          }
        }
        break;
      case 3:
        for (Z = 0; Z < Zindex; Z++) {
          for (Y = 0; Y < Yindex; Y++) {
            for (X = 0; X < Xindex; X++) {
              if ((Z * Yindex * Xindex + Y * Xindex + X) % 2 != 0) {
                Xplus1 = Z * Yindex * Xindex + Y * Xindex + (X + 1);
                Xminus1 = Z * Yindex * Xindex + Y * Xindex + (X - 1);
                Yplus1 = Z * Yindex * Xindex + (Y + 1) * Xindex + X;
                Yminus1 = Z * Yindex * Xindex + (Y - 1) * Xindex + X;
                Zplus1 = (Z + 1) * Yindex * Xindex + Y * Xindex + X;
                Zminus1 = (Z - 1) * Yindex * Xindex + Y * Xindex + X;
                if ((X * Y * Z) != 0 && X != (Xindex - 1) && Y != (Yindex - 1) && Z != (Zindex - 1)) {//not at the edge of simulation space
                  if ((geoInfo->adjacentGeo1->isDomain[Xplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Xminus1] == 1) ||
                    (geoInfo->adjacentGeo1->isDomain[Xminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Xplus1] == 1)) {//X
                      geoInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      geoInfo->domainIndex.push_back(Z * Yindex * Xindex + Y * Xindex + X);
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofXp = true;
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofXm = true;
                  } else if ((geoInfo->adjacentGeo1->isDomain[Yplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Yminus1] == 1) ||
                    (geoInfo->adjacentGeo1->isDomain[Yminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Yplus1] == 1)) {//Y
                      geoInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      geoInfo->domainIndex.push_back(Z * Yindex * Xindex + Y * Xindex + X);
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofYp = true;
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofYm = true;
                  } else if ((geoInfo->adjacentGeo1->isDomain[Zplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Zminus1] == 1) ||
                    (geoInfo->adjacentGeo1->isDomain[Zminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Zplus1] == 1)) {//Z
                      geoInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      geoInfo->domainIndex.push_back(Z * Yindex * Xindex + Y * Xindex + X);
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofZp = true;
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofZm = true;
                  }
                } else if (X == 0 || X == Xindex - 1) {
                  if ((geoInfo->adjacentGeo1->isDomain[Yplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Yminus1] == 1) ||
                    (geoInfo->adjacentGeo1->isDomain[Yminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Yplus1] == 1)) {//Y
                      geoInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      geoInfo->domainIndex.push_back(Z * Yindex * Xindex + Y * Xindex + X);
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofYp = true;
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofYm = true;
                  } else if ((geoInfo->adjacentGeo1->isDomain[Zplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Zminus1] == 1) ||
                    (geoInfo->adjacentGeo1->isDomain[Zminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Zplus1] == 1)) {//Z
                      geoInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      geoInfo->domainIndex.push_back(Z * Yindex * Xindex +Y * Xindex + X);
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofZp = true;
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofZm = true;
                  }
                } else if (Y == 0 || Y == Yindex - 1) {
                  if ((geoInfo->adjacentGeo1->isDomain[Xplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Xminus1] == 1) ||
                    (geoInfo->adjacentGeo1->isDomain[Xminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Xplus1] == 1)) {//X
                      geoInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      geoInfo->domainIndex.push_back(Z * Yindex * Xindex + Y * Xindex + X);
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofXp = true;
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofXm = true;
                  } else if ((geoInfo->adjacentGeo1->isDomain[Zplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Zminus1] == 1) ||
                    (geoInfo->adjacentGeo1->isDomain[Zminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Zplus1] == 1)) {//Y
                      geoInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      geoInfo->domainIndex.push_back(Z * Yindex * Xindex +Y * Xindex + X);
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofYp = true;
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofYm = true;
                  }
                } else if (Z == 0 || Z == Zindex - 1) {
                  if ((geoInfo->adjacentGeo1->isDomain[Xplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Xminus1] == 1) ||
                    (geoInfo->adjacentGeo1->isDomain[Xminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Xplus1] == 1)) {//X
                      geoInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      geoInfo->domainIndex.push_back(Z * Yindex * Xindex + Y * Xindex + X);
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofXp = true;
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofXm = true;
                  } else if ((geoInfo->adjacentGeo1->isDomain[Yplus1] == 1 && geoInfo->adjacentGeo2->isDomain[Yminus1] == 1) ||
                    (geoInfo->adjacentGeo1->isDomain[Yminus1] == 1 && geoInfo->adjacentGeo2->isDomain[Yplus1] == 1)) {//Y
                      geoInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      geoInfo->domainIndex.push_back(Z * Yindex * Xindex +Y * Xindex + X);
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofYp = true;
                      geoInfo->bType[Z * Yindex * Xindex + Y * Xindex + X].isBofYm = true;
                  }
                }
              }
            }
          }
        }
        //pseudo membrane
        for (Z = 0; Z < Zindex; Z++) {
          for (Y = 0; Y < Yindex; Y++) {
            for (X = 0; X < Xindex; X++) {
              if ((X % 2 != 0 && Y % 2 != 0) || (Y % 2 != 0 && Z % 2 != 0) || (Z % 2 != 0 && X % 2 != 0)) {
                index = Z * Yindex * Xindex + Y * Xindex + X;
                Xplus1 = Z * Yindex * Xindex + Y * Xindex + (X + 1);
                Xminus1 = Z * Yindex * Xindex + Y * Xindex + (X - 1);
                Yplus1 = Z * Yindex * Xindex + (Y + 1) * Xindex + X;
                Yminus1 = Z * Yindex * Xindex + (Y - 1) * Xindex + X;
                Zplus1 = (Z + 1) * Yindex * Xindex + Y * Xindex + X;
                Zminus1 = (Z - 1) * Yindex * Xindex + Y * Xindex + X;
                if (X != 0 && X != Xindex - 1 && Y != 0 && Y != Yindex - 1 && Z != 0 && Z != Zindex - 1){
                  if ((geoInfo->isDomain[Xplus1] == 1 && geoInfo->isDomain[Xminus1] == 1) ||
                    (geoInfo->isDomain[Yplus1] == 1 && geoInfo->isDomain[Yminus1] == 1) ||
                    (geoInfo->isDomain[Zplus1] == 1 && geoInfo->isDomain[Zminus1] == 1)) {
                      geoInfo->isDomain[index] = 2;
                      geoInfo->pseudoMemIndex.push_back(index);
                  }
                  if ((geoInfo->isDomain[Xplus1] == 1 && geoInfo->isDomain[Yplus1] == 1) ||
                    (geoInfo->isDomain[Xplus1] == 1 && geoInfo->isDomain[Yminus1] == 1) ||
                    (geoInfo->isDomain[Xplus1] == 1 && geoInfo->isDomain[Zplus1] == 1) ||
                    (geoInfo->isDomain[Xplus1] == 1 && geoInfo->isDomain[Zminus1] == 1) ||
                    (geoInfo->isDomain[Xminus1] == 1 && geoInfo->isDomain[Yplus1] == 1) ||
                    (geoInfo->isDomain[Xminus1] == 1 && geoInfo->isDomain[Yminus1] == 1) ||
                    (geoInfo->isDomain[Xminus1] == 1 && geoInfo->isDomain[Zplus1] == 1) ||
                    (geoInfo->isDomain[Xminus1] == 1 && geoInfo->isDomain[Zminus1] == 1) ||
                    (geoInfo->isDomain[Yplus1] == 1 && geoInfo->isDomain[Zplus1] == 1) ||
                    (geoInfo->isDomain[Yplus1] == 1 && geoInfo->isDomain[Zminus1] == 1) ||
                    (geoInfo->isDomain[Yminus1] == 1 && geoInfo->isDomain[Zplus1] == 1) ||
                    (geoInfo->isDomain[Yminus1] == 1 && geoInfo->isDomain[Zminus1] == 1)) {
                      geoInfo->isDomain[index] = 2;
                      geoInfo->pseudoMemIndex.push_back(index);
                  }
                } else if (X == 0 || X == Xindex - 1) {
                  if ((geoInfo->isDomain[Yplus1] == 1 && geoInfo->isDomain[Yminus1] == 1) ||
                    (geoInfo->isDomain[Zplus1] == 1 && geoInfo->isDomain[Zminus1] == 1) ||
                    (geoInfo->isDomain[Yplus1] == 1 && geoInfo->isDomain[Zplus1] == 1) ||
                    (geoInfo->isDomain[Yplus1] == 1 && geoInfo->isDomain[Zminus1] == 1) ||
                    (geoInfo->isDomain[Yminus1] == 1 && geoInfo->isDomain[Zplus1] == 1) ||
                    (geoInfo->isDomain[Yminus1] == 1 && geoInfo->isDomain[Zminus1] == 1)) {
                      geoInfo->isDomain[index] = 2;
                      geoInfo->pseudoMemIndex.push_back(index);
                  }
                } else if (Y == 0 || Y == Yindex - 1) {
                  if ((geoInfo->isDomain[Xplus1] == 1 && geoInfo->isDomain[Xminus1] == 1) ||
                    (geoInfo->isDomain[Zplus1] == 1 && geoInfo->isDomain[Zminus1] == 1) ||
                    (geoInfo->isDomain[Yplus1] == 1 && geoInfo->isDomain[Zplus1] == 1) ||
                    (geoInfo->isDomain[Yplus1] == 1 && geoInfo->isDomain[Zminus1] == 1) ||
                    (geoInfo->isDomain[Yminus1] == 1 && geoInfo->isDomain[Zplus1] == 1) ||
                    (geoInfo->isDomain[Yminus1] == 1 && geoInfo->isDomain[Zminus1] == 1)) {
                      geoInfo->isDomain[index] = 2;
                  }
                } else if (Z == 0 || Z == Zindex - 1) {
                  if ((geoInfo->isDomain[Xplus1] == 1 && geoInfo->isDomain[Xminus1] == 1) ||
                    (geoInfo->isDomain[Yplus1] == 1 && geoInfo->isDomain[Yminus1] == 1)) {
                      geoInfo->isDomain[index] = 2;
                      geoInfo->pseudoMemIndex.push_back(index);
                  }
                }
              }
            }
          }
        }
        break;
      default:
        break;
      }
    }
  }


  //draw geometries
  if (dimension >= 1) {
    xInfo = searchInfoById(varInfoList, xaxis);
  }
  if (dimension >= 2) {
    yInfo = searchInfoById(varInfoList, yaxis);
  }
  if (dimension >= 3) {
    zInfo = searchInfoById(varInfoList, zaxis);
  }
  vector<variableInfo*> notOrderedInfo;
  for (i = 0; i < varInfoList.size(); ++i) {
    variableInfo *info = varInfoList[i];
    ast = 0;
    if (model->getInitialAssignment(info->id) != 0) {//initial assignment
      if (info->value == 0) {//value is not set yet
        info->value = new double[numOfVolIndexes];
        memset(info->value, 0, numOfVolIndexes*sizeof(double));
        if (info->sp != 0 && (!info->sp->isSetConstant() || !info->sp->getConstant())) {
          info->delta = new double[4 * numOfVolIndexes];
          memset(info->delta, 0, 4*numOfVolIndexes*sizeof(double));
        }
      }
      ast = const_cast<ASTNode*>((model->getInitialAssignment(info->id))->getMath());
    } else if (model->getRule(info->id) != 0 && model->getRule(info->id)->isAssignment()) {//assignment rule
      info->hasAssignmentRule = true;
      if (info->value == 0) {//value is not set yet
        info->value = new double[numOfVolIndexes];
        memset(info->value, 0, numOfVolIndexes*sizeof(double));        
        if (info->sp != 0 && (!info->sp->isSetConstant() || !info->sp->getConstant())) {
          //the species is variable
          info->delta = new double[4 * numOfVolIndexes];
          memset(info->delta, 0, 4*numOfVolIndexes*sizeof(double));
        }
      }
      ast = const_cast<ASTNode*>(((AssignmentRule*)model->getRule(info->id))->getMath());
    }


     //if (info->sp != 0) 
     //{
     //  if (info->geoi == NULL)
     //    info->geoi = searchAvolInfoByCompartment(geoInfoList, info->sp->getCompartment().c_str());
     //   if (info->sp != 0 && info->delta == NULL && (!info->sp->isSetConstant() || !info->sp->getConstant())) {
     //     //the species is variable
     //     info->delta = new double[4 * numOfVolIndexes];
     //     fill_n(info->delta, 4 * numOfVolIndexes, 0);
     //   }
     //}
    
    if (ast != 0) {
      rearrangeAST(ast);
      numOfASTNodes = 0;
      countAST(ast, numOfASTNodes);
      info->rpInfo = new reversePolishInfo();
      info->rpInfo->varList = new double*[numOfASTNodes];      
      memset(info->rpInfo->varList, 0, numOfASTNodes*sizeof(double*));
      info->rpInfo->deltaList = 0;
      info->rpInfo->constList = new double*[numOfASTNodes];
      memset(info->rpInfo->constList, 0, numOfASTNodes*sizeof(double*));
      info->rpInfo->opfuncList = new int[numOfASTNodes];
      memset(info->rpInfo->opfuncList, 0, numOfASTNodes*sizeof(int));
      info->rpInfo->listNum = numOfASTNodes;
      info->isResolved = false;
      parseDependence(ast, info->dependence, varInfoList);
      notOrderedInfo.push_back(info);
    }
  }
  //dependency of symbols
  unsigned int resolved_count = 0;

  if (notOrderedInfo.size() != 0) {
    for (i = 0;i < notOrderedInfo.size();++i) {
      variableInfo *info = notOrderedInfo[i];
      if (isResolvedAll(info->dependence) && info->isResolved == false) {
        std::string formula; 
        if (model->getInitialAssignment(info->id) != 0) {//initial assignment
          ast = const_cast<ASTNode*>((model->getInitialAssignment(info->id))->getMath());
          formula = SBML_formulaToString((model->getInitialAssignment(info->id))->getMath());
        } else if (model->getRule(info->id) != 0 && model->getRule(info->id)->isAssignment()) {//assignment rule
          ast = const_cast<ASTNode*>(((AssignmentRule*)model->getRule(info->id))->getMath());
          formula = SBML_formulaToString(((AssignmentRule*)model->getRule(info->id))->getMath());
        }
        parseAST(ast, info->rpInfo, varInfoList, info->rpInfo->listNum, freeConstList);
        std::string formula2 = SBML_formulaToString(ast); 
        bool isAllArea = (info->sp != 0)? false: true;
        if (info->sp != 0) info->geoi = searchAvolInfoByCompartment(geoInfoList, info->sp->getCompartment().c_str());
        else if (info->com != NULL) info->geoi = searchAvolInfoByCompartment(geoInfoList, info->com->getId().c_str());
        else
        {
          string id = getFirstCompartmentId(model, ast);
          if (!id.empty())
            info->geoi = searchAvolInfoByCompartment(geoInfoList, id.c_str());          
        }
        if (info->geoi == NULL)
        {
          vector<int> vTemp;
          reversePolishInitial(vTemp, info->rpInfo, info->value, info->rpInfo->listNum, Xindex, Yindex, Zindex, isAllArea);
        }
        else
          reversePolishInitial(info->geoi->domainIndex, info->rpInfo, info->value, info->rpInfo->listNum, Xindex, Yindex, Zindex, isAllArea);
        info->isResolved = true;
        if (info->hasAssignmentRule) orderedARule.push_back(info);
      }
      if (i == notOrderedInfo.size() - 1) {
        for (j = 0, resolved_count = 0; j < notOrderedInfo.size(); j++) {
          if (notOrderedInfo[j]->isResolved == true) resolved_count++;
        }
        if (resolved_count == notOrderedInfo.size()) break;
        else i = 0;
      }
    }
  }

  //calc normal unit vector of membrane (for mem diffusion and mem transport)
  nuVec = 0;
  vorI = 0;
  if (dimension >= 2) {
    //calc normalUnitVector at membrane
    nuVec = setNormalAngle(geoInfoList, Xsize, Ysize, Zsize, dimension, Xindex, Yindex, Zindex, numOfVolIndexes);
    //calc voronoi at membrane
    vorI = setVoronoiInfo(nuVec, xInfo, yInfo, zInfo, geoInfoList, Xsize, Ysize, Zsize, dimension, Xindex, Yindex, Zindex, numOfVolIndexes);
  }

  //set boundary type
  setBoundaryType(model, varInfoList, geoInfoList, Xindex, Yindex, Zindex, dimension);
  //double dt = 0.1;
  ////numerical stability analysis of diffusion and advection
  //double min_dt = dt;
  //for (i = 0; i < numOfSpecies; i++) {
  //  variableInfo *sInfo = searchInfoById(varInfoList, los->get(i)->getId().c_str());
  //  if (sInfo==NULL) continue;
  //  //volume diffusion
  //  if (sInfo->diffCInfo != 0 && sInfo->geoi->isVol) {
  //    min_dt = min(min_dt, checkDiffusionStab(sInfo, deltaX, deltaY, deltaZ, Xindex, Yindex, dt));
  //  }
  //  //membane diffusion
  //  if (sInfo->diffCInfo != 0 && !sInfo->geoi->isVol) {
  //    min_dt = min(min_dt, checkMemDiffusionStab(sInfo, vorI, Xindex, Yindex, dt, dimension));
  //  }
  //  //advection
  //  if (sInfo->adCInfo != 0) {
  //    min_dt = min(min_dt, checkAdvectionStab(sInfo, deltaX, deltaY, deltaZ, dt, Xindex, Yindex, dimension));
  //  }
  //}
  //cout << "finished" << endl;
  //if (dt > min_dt) {
  //  cout << "dt must be less than " << min_dt << endl;
  //  exit(1);
  //}



  //reaction information
  setReactionInfo(model, varInfoList, rInfoList, fast_rInfoList, freeConstList, numOfVolIndexes);

  //rate rule information
  setRateRuleInfo(model, varInfoList, rInfoList, freeConstList, numOfVolIndexes);

  //output geometries
  int *geo_edge = new int[numOfVolIndexes];
  memset(geo_edge, 0, numOfVolIndexes*sizeof(int));  
  vector<GeometryInfo*> memInfoList = vector<GeometryInfo*>();
  for (i = 0; i < geoInfoList.size(); i++) {
    GeometryInfo *geoInfo = geoInfoList[i];
    if (geoInfo->isVol == false) {//avol is membrane
      memList.push_back(geoInfo->domainTypeId);
      memInfoList.push_back(geoInfo);
      for (Z = 0; Z < Zindex; Z++) {
        for (Y = 0; Y < Yindex; Y++) {
          for (X = 0; X < Xindex; X++) {
            if (geoInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] > 0) {
              if ((Z * Yindex * Xindex + Y * Xindex + X) % 2 != 0) geo_edge[Z * Yindex * Xindex + Y * Xindex + X] = 1;
              else geo_edge[Z * Yindex * Xindex + Y * Xindex + X] = 2;
            }
          }
        }
      }
    }
  }


}

SpatialSimulator::~SpatialSimulator()
{
  freeVarInfo(varInfoList);
  freeAvolInfo(geoInfoList);
  freeRInfo(rInfoList);
  for (size_t i = 0; i < freeConstList.size(); ++i) {
    delete freeConstList[i];
    freeConstList[i] = NULL;
  }  
}

double SpatialSimulator::oneStep(double initialTime, double step)
{
  //printValues();
  performStep(initialTime, step);
  updateValues(step);
  return initialTime + step;
}

double* SpatialSimulator::getX(int &length)
{
  if (xInfo == NULL)
  {
    length = 0; return NULL;
  }
  length = (Zindex -1) * Yindex * Xindex + (Yindex-1) * Xindex + Xindex-1;
  return xInfo->value;

}

double* SpatialSimulator::getY(int &length)
{
  if (yInfo == NULL)
  {
    length = 0; return NULL;
  }
  length = (Zindex -1) * Yindex * Xindex + (Yindex-1) * Xindex + Xindex-1;
  return yInfo->value;
}

double* SpatialSimulator::getZ(int &length)
{
  if (zInfo == NULL)
  {
    length = 0; return NULL;
  }
  length = (Zindex -1) * Yindex * Xindex + (Yindex-1) * Xindex + Xindex-1;
  return zInfo->value;
}

int* SpatialSimulator::getGeometry(const std::string &compartmentId, int &length)
{
  GeometryInfo* info = searchAvolInfoByCompartment(geoInfoList, compartmentId.c_str());
  if (info == NULL) 
  {
    length = 0;
    return NULL;
  }
  length = numOfVolIndexes;
  return info->isDomain;
}

int SpatialSimulator::getGeometryLength() const
{
  return numOfVolIndexes;
}

boundaryType* SpatialSimulator::getBoundaryType(const std::string &compartmentId, int &length)
{
  GeometryInfo* info = searchAvolInfoByCompartment(geoInfoList, compartmentId.c_str());
  if (info == NULL)
  {
    length = 0;
    return NULL;
  }
  length = numOfVolIndexes;
  return info->bType;
}

int* SpatialSimulator::getBoundary(const std::string &compartmentId, int &length)
{
  GeometryInfo* info = searchAvolInfoByCompartment(geoInfoList, compartmentId.c_str());
   if (info == NULL)
  {
    length = 0;
    return NULL;
  }
  length = numOfVolIndexes;
  return info->isBoundary;
}

double SpatialSimulator::getVariableAt(const std::string& variable, int x, int y, int z)
{  
  variableInfo *sInfo = searchInfoById(varInfoList, variable.c_str());
  if (sInfo == NULL)
  {
    return 0.0;
  }

  int index = getIndexForPosition(x, y);
  if (index == -1)
    return 0.0;

  if (sInfo->geoi == NULL || sInfo->geoi->isDomain[index] == 0)
    return 0.0;

  return sInfo->value[index];
}

double* SpatialSimulator::getVariable(const std::string &speciesId, int &length)
{
  variableInfo *sInfo = searchInfoById(varInfoList, speciesId.c_str());
  if (sInfo == NULL)
  {
    length = 0; return NULL;
  }
  length = (Zindex -1) * Yindex * Xindex + (Yindex-1) * Xindex + Xindex-1;
  return sInfo->value;
}

void SpatialSimulator::printValues() 
{

}

void SpatialSimulator::performStep(double t, double dt)
{
  unsigned int /* i,*/ j , m;
  //int X = 0, Y = 0, Z = 0;
  //int index = 0;

  //int Xplus = 0, Xminus = 0, Yplus = 0, Yminus = 0, Zplus = 0, Zminus = 0;
  //double end_time=1000;
  sim_time = t * dt;
  //output
  //if (count % out_step == 0) {
  //  // funa
  //	outputTimeCource(gp, model, varInfoList, memList, xInfo, yInfo, zInfo, &sim_time, end_time, dt, range_max, dimension, Xindex, Yindex, Zindex, Xsize, Ysize, Zsize, file_num, fname);
  //	file_num++;
  //}

  //calculation

  //runge-kutta
  for (m = 0; m < 4; m++) {

    //reaction
    //for (i = 0; i < numOfReactions; i++) {
    //slow reaction
    for (size_t i = 0; i < rInfoList.size(); i++) {
      //Reaction *r = model->getReaction(i);
      Reaction *r = rInfoList[i]->reaction;
      if (r == NULL) continue;
      if (!rInfoList[i]->isMemTransport) 
      {//normal reaction

        GeometryInfo* gInfo = NULL;
        // get first non-null geometry
        for (size_t k = 0; k < rInfoList[i]->spRefList.size(); ++k)
          if (rInfoList[i]->spRefList[k]->geoi != NULL)
          {
            gInfo = rInfoList[i]->spRefList[k]->geoi;
            break;
          }

          reversePolishRK(rInfoList[i], gInfo, Xindex, Yindex, Zindex, dt, m, r->getNumReactants(), true);
      } 
      else 
      {//membrane transport
        GeometryInfo *reactantGeo = searchInfoById(varInfoList, r->getReactant(0)->getSpecies().c_str())->geoi;
        GeometryInfo *productGeo = searchInfoById(varInfoList, r->getProduct(0)->getSpecies().c_str())->geoi;

        if (reactantGeo == NULL)
        {
          const std::string& speciesId = r->getReactant(0)->getSpecies();
          const std::string& compartment = model->getSpecies(speciesId)->getCompartment();
          for (size_t i = 0; i < model->getNumSpecies(); ++i)
          {
            const Species* current = model->getSpecies(i);
            if (current->getId() != speciesId && current->getCompartment() == compartment)
            {
              reactantGeo = searchInfoById(varInfoList, current->getId().c_str())->geoi;
              if (reactantGeo != NULL) 
                break;
            }
          }
        }

        if (productGeo == NULL)
        {
          const std::string& speciesId = r->getProduct(0)->getSpecies();
          const std::string& compartment = model->getSpecies(speciesId)->getCompartment();
          for (size_t i = 0; i < model->getNumSpecies(); ++i)
          {
            const Species* current = model->getSpecies(i);
            if (current->getId() != speciesId && current->getCompartment() == compartment)
            {
              productGeo = searchInfoById(varInfoList, current->getId().c_str())->geoi;
              if (productGeo != NULL) 
                break;
            }
          }
        }

        for (j = 0; j < geoInfoList.size(); j++) {
          if (!geoInfoList[j]->isVol) {
            if ((geoInfoList[j]->adjacentGeo1 == reactantGeo && geoInfoList[j]->adjacentGeo2 == productGeo)
              || (geoInfoList[j]->adjacentGeo1 == productGeo && geoInfoList[j]->adjacentGeo2 == reactantGeo)) {//mem transport
                //cout << geoInfoList[j]->domainTypeId << endl;
                calcMemTransport(rInfoList[i], geoInfoList[j], nuVec, Xindex, Yindex, Zindex, dt, m, deltaX, deltaY, deltaZ, dimension, r->getNumReactants());
                break;
            }
          }
        }

        bool haveReactantVol = reactantGeo != NULL && reactantGeo->isVol;
        bool haveProductVol = productGeo != NULL && productGeo->isVol;

        if (haveReactantVol ^ haveProductVol) {
          if (!haveReactantVol) {
            calcMemTransport(rInfoList[i], reactantGeo, nuVec, Xindex, Yindex, Zindex, dt, m, deltaX, deltaY, deltaZ, dimension, r->getNumReactants());
          }
          if (!haveProductVol) {
            calcMemTransport(rInfoList[i], productGeo, nuVec, Xindex, Yindex, Zindex, dt, m, deltaX, deltaY, deltaZ, dimension, r->getNumReactants());
          }
        }
      }
    }

    //rate rule
    for (size_t i = 0; i < numOfRules; ++i) {
      if (model->getRule(i)->isRate()) {
        RateRule *rrule = (RateRule*)model->getRule(i);
        variableInfo *sInfo = searchInfoById(varInfoList, rrule->getVariable().c_str());
        reversePolishRK(rInfoList[i], sInfo->geoi, Xindex, Yindex, Zindex, dt, m, 1, false);
      }
    }

    //diffusion
    for (size_t i = 0; i < numOfSpecies; ++i) {
      variableInfo *sInfo = searchInfoById(varInfoList, los->get(i)->getId().c_str());
      if (sInfo == NULL) continue;
      //volume diffusion
      if (sInfo->diffCInfo != 0 && sInfo->geoi->isVol) {
        calcDiffusion(sInfo, deltaX, deltaY, deltaZ, Xindex, Yindex, Zindex, m, dt, dimension);
      }
      //membane diffusion
      if (sInfo->diffCInfo != 0 && !sInfo->geoi->isVol) {
        calcMemDiffusion(sInfo, vorI, Xindex, Yindex, Zindex, m, dt, dimension);
      }
      //boundary condition
      if (sInfo->boundaryInfo != 0 && sInfo->geoi->isVol) {
        calcBoundary(sInfo, deltaX, deltaY, deltaZ, Xindex, Yindex, Zindex, m, dimension);
      }
    }

  }//end of runge-kutta

  //advection
  for (size_t i = 0; i < numOfSpecies; i++) {
    variableInfo *sInfo = searchInfoById(varInfoList, los->get(i)->getId().c_str());
    if (sInfo == NULL) continue;
    //advection
    if (sInfo->adCInfo != 0) {
      cipCSLR(sInfo, deltaX, deltaY, deltaZ, dt, Xindex, Yindex, Zindex, dimension);
    }//end of advection
  }
}

void SpatialSimulator::updateValues(double dt)
{
  if (los == NULL) return;
#pragma region     // update values
  //update values
#pragma omp parallel for 
  for (int i = 0; i < (int)numOfSpecies; i++) 
  {
    s = los->get(i);
    if (s == NULL) continue;
    if (s->isSetConstant() && s->getConstant())
      continue;

    variableInfo *sInfo = searchInfoById(varInfoList, s->getId().c_str());
    updateSpecies(sInfo, dt);
  }

  updateAssignmentRules();

#pragma endregion 


}

void SpatialSimulator::updateSpecies(variableInfo *sInfo, double dt)
{
  if (sInfo == NULL || sInfo->geoi == NULL) 
    return;
  int X = 0, Y = 0, Z = 0, index = 0;
  if (sInfo->sp->isSetConstant() && sInfo->sp->getConstant()) 
    return;

  //update values (advection, diffusion, slow reaction)
  for (size_t j = 0; j < sInfo->geoi->domainIndex.size(); ++j) 
  {
    index = sInfo->geoi->domainIndex[j];
    Z = index / (Xindex * Yindex);
    Y = (index - Z * Xindex * Yindex) / Xindex;
    X = index - Z * Xindex * Yindex - Y * Xindex;
    divIndex = (Z / 2) * Ydiv * Xdiv + (Y / 2) * Xdiv + (X / 2);
    //update values for the next time
    sInfo->value[index] += dt * (sInfo->delta[index] + 2.0 * sInfo->delta[numOfVolIndexes + index] + 2.0 * sInfo->delta[2 * numOfVolIndexes + index] + sInfo->delta[3 * numOfVolIndexes + index]) / 6.0;
    for (size_t k = 0; k < 4; ++k) sInfo->delta[k * numOfVolIndexes + index] = 0.0;
  }

  //boundary condition
  if (sInfo->boundaryInfo != 0) {
    calcBoundary(sInfo, deltaX, deltaY, deltaZ, Xindex, Yindex, Zindex, 0, dimension);
  }    

}

void SpatialSimulator::deleteValuesOutsideDomain(const std::string& id)
{
  deleteValuesOutsideDomain(searchInfoById(varInfoList, id.c_str()));
}

void SpatialSimulator::deleteValuesOutsideDomain(variableInfo *info )
{
  if (info == NULL) return;
  for (int Z = 0; Z < Zindex; Z += 2) {
    for (int Y = 0; Y < Yindex; Y += 2) {
      for (int X = 0; X < Xindex; X += 2) {
        int index = Z * Yindex * Xindex + Y * Xindex + X;
        if (static_cast<int>(info->geoi->isDomain[index]) == 0) {
          info->value[index] = 0.0;
        }
      }
    }
  }
}

void SpatialSimulator::updateAssignmentRules()
{

  int X = 0, Y = 0, Z = 0, index = 0;
  //fast reaction
  // 		for (i = 0; i < fast_rInfoList.size(); i++) {
  // 			Reaction *r = fast_rInfoList[i]->reaction;
  // 		}
  //assignment rule
  for (size_t i = 0; i < orderedARule.size(); ++i) {
    variableInfo *info = orderedARule[i];
    bool isAllArea = (info->sp != 0)? false: true;
    if (info->sp != 0) {
      info->geoi = searchAvolInfoByCompartment(geoInfoList, info->sp->getCompartment().c_str());
    }
    reversePolishInitial(info->geoi->domainIndex, info->rpInfo, info->value, info->rpInfo->listNum, Xindex, Yindex, Zindex, isAllArea);
    /*
    double tmp_x, tmp_y, tmp_z, tmp_len;
    double tmp_phi, tmp_theta;
    for (j = 0; j < info->geoi->domainIndex.size(); j++) {
    tmp_x = xInfo->value[info->geoi->domainIndex[j]];
    tmp_y = yInfo->value[info->geoi->domainIndex[j]];
    tmp_z = zInfo->value[info->geoi->domainIndex[j]];
    tmp_phi = atan2(tmp_y, tmp_x);
    tmp_theta = acos(tmp_z);
    //tmp_x = sin(tmp_theta) * cos(tmp_phi);
    //tmp_y = sin(tmp_theta) * sin(tmp_phi);
    //tmp_z = cos(tmp_theta);
    //tmp_len = 1.0 - sqrt(pow(tmp_x, 2) + pow(tmp_y, 2) + pow(tmp_z, 2));
    //tmp_x += tmp_len * nuVec[info->geoi->domainIndex[j]].nx;
    //tmp_y += tmp_len * nuVec[info->geoi->domainIndex[j]].ny;
    //tmp_z += tmp_len * nuVec[info->geoi->domainIndex[j]].nz;
    //info->value[info->geoi->domainIndex[j]] = 2.0 + 0.5 * exp(-2.0 * (*sim_time)) * tmp_x + exp(-6.0 * (*sim_time)) * (pow(tmp_x, 2) - pow(tmp_y, 2)) + exp(-12.0 * (*sim_time)) * (pow(tmp_x, 3) - 3 * tmp_x * pow(tmp_y, 2));
    //if ((tmp_x - 5.0 * (*sim_time) >= 10.0) && (tmp_x - 5.0 * (*sim_time) <= 30.0) && (tmp_y - 5.0 * (*sim_time) >= 10.0) && (tmp_y - 5.0 * (*sim_time) <= 30.0)) {
    // 				if ((tmp_x - 5.0 * (*sim_time) >= 10.0) && (tmp_x - 5.0 * (*sim_time) <= 30.0) && (tmp_y - 5.0 * (*sim_time) >= 10.0) && (tmp_y - 5.0 * (*sim_time) <= 30.0)) {
    // 					info->value[info->geoi->domainIndex[j]] = 5.0;
    // 				} else {
    // 					info->value[info->geoi->domainIndex[j]] = 0.0;
    // 				}

    //info->value[info->geoi->domainIndex[j]] = 5.0 * exp(-pow(0.1 * (-30.0 + tmp_x - 5.0 *sim_time), 2) + pow(0.1 * (-30.0 + tmp_y + t * t * -1), 2) * -1);
    //cout << info->value[info->geoi->domainIndex[j]] << endl;
    //cout << tmp_len * nuVec[index].nx << endl;
    //cout << sqrt(pow(tmp_x, 2) + pow(tmp_y, 2) + pow(tmp_z, 2)) << endl;
    }
    */
  }
  //pseudo membrane
  for (size_t i = 0; i < numOfSpecies; ++i) {
    s = los->get(i);
    variableInfo *sInfo = searchInfoById(varInfoList, s->getId().c_str());
    if (sInfo == NULL|| sInfo->geoi == NULL) continue;
    if (!sInfo->geoi->isVol) {
      for (size_t j = 0; j < sInfo->geoi->pseudoMemIndex.size(); ++j) {
        index = sInfo->geoi->pseudoMemIndex[j];
        Z = index / (Xindex * Yindex);
        Y = (index - Z * Xindex * Yindex) / Xindex;
        X = index - Z * Xindex * Yindex - Y * Xindex;
        Xplus1 = Z * Yindex * Xindex + Y * Xindex + (X + 1);
        Xminus1 = Z * Yindex * Xindex + Y * Xindex + (X - 1);
        Yplus1 = Z * Yindex * Xindex + (Y + 1) * Xindex + X;
        Yminus1 = Z * Yindex * Xindex + (Y - 1) * Xindex + X;
        Zplus1 = (Z + 1) * Yindex * Xindex + Y * Xindex + X;
        Zminus1 = (Z - 1) * Yindex * Xindex + Y * Xindex + X;
        if (sInfo->geoi->isDomain[Xplus1] == 1 && sInfo->geoi->isDomain[Xminus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Xminus1]);
        } else if (sInfo->geoi->isDomain[Yplus1] == 1 && sInfo->geoi->isDomain[Yminus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Yplus1], sInfo->value[Yminus1]);
        } else if (dimension == 3 && sInfo->geoi->isDomain[Zplus1] == 1 && sInfo->geoi->isDomain[Zminus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Zplus1], sInfo->value[Zminus1]);
        } else if (sInfo->geoi->isDomain[Xplus1] == 1 && sInfo->geoi->isDomain[Yplus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Yplus1]);
        } else if (sInfo->geoi->isDomain[Xplus1] == 1 && sInfo->geoi->isDomain[Yminus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Yminus1]);
        } else if (dimension == 3 && sInfo->geoi->isDomain[Xplus1] == 1 && sInfo->geoi->isDomain[Zplus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Zplus1]);
        } else if (dimension == 3 && sInfo->geoi->isDomain[Xplus1] == 1 && sInfo->geoi->isDomain[Zminus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Zminus1]);
        } else if (sInfo->geoi->isDomain[Xminus1] == 1 && sInfo->geoi->isDomain[Yplus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Xminus1], sInfo->value[Yplus1]);
        } else if (sInfo->geoi->isDomain[Xminus1] == 1 && sInfo->geoi->isDomain[Yminus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Xminus1], sInfo->value[Yminus1]);
        } else if (dimension == 3 && sInfo->geoi->isDomain[Xminus1] == 1 && sInfo->geoi->isDomain[Zplus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Xminus1], sInfo->value[Zplus1]);
        } else if (dimension == 3 && sInfo->geoi->isDomain[Xminus1] == 1 && sInfo->geoi->isDomain[Zminus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Xminus1], sInfo->value[Zminus1]);
        } else if (dimension == 3 && sInfo->geoi->isDomain[Yplus1] == 1 && sInfo->geoi->isDomain[Zplus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Yplus1], sInfo->value[Zplus1]);
        } else if (dimension == 3 && sInfo->geoi->isDomain[Yplus1] == 1 && sInfo->geoi->isDomain[Zminus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Yplus1], sInfo->value[Zminus1]);
        } else if (dimension == 3 && sInfo->geoi->isDomain[Yminus1] == 1 && sInfo->geoi->isDomain[Zplus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Yminus1], sInfo->value[Zplus1]);
        } else if (dimension == 3 && sInfo->geoi->isDomain[Yminus1] == 1 && sInfo->geoi->isDomain[Zminus1] == 1) {
          sInfo->value[index] = min(sInfo->value[Yminus1], sInfo->value[Zminus1]);
        }
      }
    }
  }

}

void SpatialSimulator::run(double endTime, double step)
{
  double end_time = endTime; 
  double dt = step;
  double out_time = 1.0;
  int t = 0, count = 0,  percent = 0;	

  file_num = 0;

  ofstream ofs;
  string voidString;
#if WIN32
  voidString = "> NUL";
#endif
  gp = popen((gnuplotExecutable + " -persist " + voidString).c_str(),"w");
  for (t = 0; t <= static_cast<int>(end_time / dt); t++) {
    //cerr << static_cast<int>(100.0 * static_cast<double>(t) / (end_time / dt)) << endl;
    sim_time = t * dt;


    //output
    if (count % static_cast<int>(out_time / dt) == 0) {
      printValues();
    }

    count++;

    oneStep(t, dt);


    if (t == (static_cast<int>(end_time / dt) / 10) * percent) {
      cout << percent * 10 << "% finished" << endl;
      percent++;
    }
  }
  pclose(gp);
}

void  SpatialSimulator::setGnuplotExecutable(const std::string &location)
{
  gnuplotExecutable = location;
}

const std::string& SpatialSimulator::getGnuplotExecutable() const
{
  return gnuplotExecutable;
}

const Model* SpatialSimulator::getModel() const
{
  return model;
}

struct sort_pair
{
  bool operator() (std::pair<int, int> i, std::pair<int, int> j) 
  {
    return (i.second < j.second);
  }
} sorter;


void flipOrder(AnalyticGeometry& geometry)
{
  if (geometry.getNumAnalyticVolumes() == 0) return;

  std::vector< std::pair<int, int> > list;
  for (size_t i = 0; i < geometry.getNumAnalyticVolumes(); ++i)
  {
    AnalyticVolume* current = geometry.getAnalyticVolume(i);
    list.push_back(std::pair<int, int>((int)i, (int)current->getOrdinal()));
  }

  std::sort(list.begin(), list.end(), sorter);

  for (std::vector< std::pair<int, int> >::iterator it=list.begin(); it!=list.end(); ++it)
  {
    AnalyticVolume* current = geometry.getAnalyticVolume(it->first);
    current->setOrdinal((geometry.getNumAnalyticVolumes()-1)-it->second);
  }

}

void SpatialSimulator::flipVolumeOrder()
{    
  if (model == NULL)
    return;
  SpatialModelPlugin* plugin = dynamic_cast<SpatialModelPlugin*>(model->getPlugin("spatial"));
  if (plugin == NULL)
    return;

  Geometry* geometry = plugin->getGeometry();

  for (size_t i = 0; i < geometry->getNumGeometryDefinitions(); ++i)
  {
    AnalyticGeometry* definition = dynamic_cast<AnalyticGeometry*>(geometry->getGeometryDefinition(i));
    if (definition == NULL) 
      continue;
    flipOrder(*definition);    
  }
}

int SpatialSimulator::getIndexForPosition(double x, double y)
{
  for (int Y = 0; Y < Yindex; Y += 2 ) {
    for (int X = 0; X < Xindex; X += 2 ) {
      int index = Y * Xindex + X;
      double currentX =  xInfo->value[index];
      double currentY =  yInfo->value[index];

      if (currentX == x && currentY == y)
        return index;

    }
  }
  return -1;
}
