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

#include <cmath>
#include <ctime>
#include <sys/stat.h>

#include <sbml/SBMLTypes.h>

#include <sbml/common/extern.h>
#include <sbml/common/common.h>
#include <sbml/common/libsbml-namespace.h>

#include <sbml/extension/SBMLExtensionRegistry.h>

#include <sbml/packages/req/common/RequiredElementsExtensionTypes.h>
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


void spatialSimulator(SBMLDocument *doc, int argc, const char *argv[]);
void rearrangeAST(ASTNode *ast);
void countAST(ASTNode *ast, int &numOfASTNode);
void parseAST(ASTNode *ast, double **&varList, double **&deltaList, double **&constList, int *&opfuncList, vector<variableInfo*> &varInfoList, int index_max);
void InitializeVarInfo(variableInfo *&varInfo);
void InitializeBcOfSpeciesInfo(bcOfSpeciesInfo *&bcos);
void InitializeAVolInfo(analyticVolInfo *&avolInfo);
variableInfo* searchInfoById(vector<variableInfo*> &varInfoList, const char *varName);
bcOfSpeciesInfo* searchBCOSInfoBySpeId(vector<bcOfSpeciesInfo*> &bcOfSpeciesInfoList, const char *varName);
analyticVolInfo* searchAvolInfoByDomainType(vector<analyticVolInfo*> &avolInfoList, const char *dtId);
analyticVolInfo* searchAvolInfoByDomain(vector<analyticVolInfo*> &avolInfoList, const char *dId);
analyticVolInfo* searchAvolInfoByCompartment(vector<analyticVolInfo*> &avolInfoList, const char *cId);
void reversePolishInitial(double **&varList, double **&constList, int *&opfuncList, double *&value, int numOfASTNodes, int Xindex, int Yindex, int Zindex, bool inVol);
void reversePolishRK(reactionInfo *&rInfo, variableInfo *&sInfo, int Xindex, int Ydiv, int Zdiv, double dt, int m, materialType mType, bool inVol, diffCInfo** dci);
void freeVarInfo(vector<variableInfo*> &varInfoList, double *sim_time_p);
void freeBcOfSpeciesInfo(vector<bcOfSpeciesInfo*> &bcOfSpeciesInfoList);
void freeAvolInfo(vector<analyticVolInfo*> &avolInfoList);
void freeRInfo(vector<reactionInfo*> &avolInfoList);
void freeAvolInfo(vector<analyticVolInfo*> &avolInfoList);
void parseDependence(const ASTNode *ast, vector<variableInfo*> &dependence, vector<variableInfo*> &varInfoList);
bool isResolvedAll(vector<variableInfo*> &dependence);

// moved old version out of the way
#include "spatialsimulator_legacy.cxx"

void InitializeVarInfo(variableInfo *&varInfo)
{
  varInfo->sp = NULL;
  varInfo->para = NULL;
  varInfo->com = NULL;
  varInfo->avi = NULL;
  varInfo->delta = NULL;
  varInfo->id = NULL;
  varInfo->dcip = NULL;
  varInfo->inVol = true;
  varInfo->hasAssignmentRule = false;
  //varInfo->next = 0;
  varInfo->rpInfo = NULL;
  varInfo->value = NULL;
  varInfo->dci = NULL;
  varInfo->aci = NULL;
  varInfo->isResolved = false;
  varInfo->dim = 0;  
}

void InitializeBcOfSpeciesInfo(bcOfSpeciesInfo *&bcos)
{
  bcos->bcXp = NULL;
  bcos->bcXm = NULL;
  bcos->bcYp = NULL;
  bcos->bcYm = NULL;
  bcos->bcZp = NULL;
  bcos->bcZm = NULL;
}

void InitializeAVolInfo(analyticVolInfo *&avolInfo)
{
  avolInfo->compartmentId = NULL;
  avolInfo->domainTypeId = NULL;
  avolInfo->domainId = NULL;
  avolInfo->adjacent0 = NULL;
  avolInfo->adjacent1 = NULL;
  avolInfo->bt = NULL;
  avolInfo->isVol = true;
  avolInfo->implicit = true;
  avolInfo->isDomain = NULL;
  avolInfo->isBoundary = NULL;
  avolInfo->isEdge = NULL;
  avolInfo->rpInfo = NULL;
}

void parseAST(ASTNode *ast, double **&varList, double **&deltaList, double **&constList, int *&opfuncList, vector<variableInfo*> &varInfoList, int index_max)
{
  static int index = 0;
  unsigned int i;
  for (i = 0; i < ast->getNumChildren(); i++) {
    parseAST(ast->getChild(i), varList, deltaList, constList, opfuncList, varInfoList, index_max);
  }
  if (ast->isFunction() || ast->isOperator() || ast->isRelational() || ast->isLogical()) {
    //ast is function, operator, relational or logical
    varList[index] = 0;
    constList[index] = 0;
    opfuncList[index] = ast->getType();
    if (deltaList != 0) {
      deltaList[index] = 0;
    }
  } else if (ast->isReal()) {//ast is real number
    varList[index] = 0;
    constList[index] = new double(ast->getReal());
    opfuncList[index] = 0;
    if (deltaList != 0) {
      deltaList[index] = 0;
    }
  } else if (ast->isInteger()) {//ast is integer
    varList[index] = 0;
    constList[index] = new double((double)ast->getInteger());
    opfuncList[index] = 0;
    if (deltaList != 0) {
      deltaList[index] = 0;
    }
  } else if (ast->isConstant()) {//ast is constant
    ASTNodeType_t type = ast->getType();
    if (type == AST_CONSTANT_E) {
      varList[index] = 0;
      constList[index] = new double(M_E);
      opfuncList[index] = 0;
      if (deltaList != 0) {
        deltaList[index] = 0;
      }
    } else if (type == AST_CONSTANT_PI) {
      varList[index] = 0;
      constList[index] = new double(M_PI);
      opfuncList[index] = 0;
      if (deltaList != 0) {
        deltaList[index] = 0;
      }
    }  else if (type == AST_CONSTANT_FALSE) {
      varList[index] = 0;
      constList[index] = new double(0.0);
      opfuncList[index] = 0;
      if (deltaList != 0) {
        deltaList[index] = 0;
      }
    }  else if (type == AST_CONSTANT_TRUE) {
      varList[index] = 0;
      constList[index] = new double(1.0);
      opfuncList[index] = 0;
      if (deltaList != 0) {
        deltaList[index] = 0;
      }
    }
  } else {//variable
    variableInfo *info = searchInfoById(varInfoList, ast->getName());
    if (info != 0) {
      if (strcmp(ast->getName(), "t") != 0) {
        varList[index] = info->value;
        constList[index] = 0;
        opfuncList[index] = 0;
        if (deltaList != 0) {
          deltaList[index] = info->delta;
        }
      } else {
        varList[index] = 0;
        constList[index] = info->value;
        opfuncList[index] = 0;
        if (deltaList != 0) {
          deltaList[index] = 0;
        }
      }
    } else {
      //cout << "couldn't find value of \"" << ast->getName() << "\"" << endl;
    }
  }
  index++;
  if (index == index_max) {
    index = 0;
  }
}

variableInfo* searchInfoById(vector<variableInfo*> &varInfoList, const char *varName)
{
  vector<variableInfo*>::iterator it = varInfoList.begin();
  while (it != varInfoList.end()) {
    if (strcmp((*it)->id, varName) == 0) {
      return *it;
    }
    it++;
  }
  return NULL;
}

bcOfSpeciesInfo* searchBCOSInfoBySpeId(vector<bcOfSpeciesInfo*> &bcOfSpeciesInfoList, const char *varName)
{
  vector<bcOfSpeciesInfo*>::iterator it = bcOfSpeciesInfoList.begin();
  while (it != bcOfSpeciesInfoList.end()) {
    if (strcmp((*it)->speciesId, varName) == 0) {
      return *it;
    }
    it++;
  }
  return NULL;
}

analyticVolInfo* searchAvolInfoByDomainType(vector<analyticVolInfo*> &avolInfoList, const char *dtId)
{
  vector<analyticVolInfo*>::iterator it = avolInfoList.begin();
  while (it != avolInfoList.end()) {
    if (strcmp((*it)->domainTypeId, dtId) == 0) {
      return *it;
    }
    it++;
  }
  return NULL;
}

analyticVolInfo* searchAvolInfoByDomain(vector<analyticVolInfo*> &avolInfoList, const char *dId)
{
  vector<analyticVolInfo*>::iterator it = avolInfoList.begin();
  while (it != avolInfoList.end()) {
    if (strcmp((*it)->domainId, dId) == 0) {
      return *it;
    }
    it++;
  }
  return NULL;
}

analyticVolInfo* searchAvolInfoByCompartment(vector<analyticVolInfo*> &avolInfoList, const char *cId)
{
  vector<analyticVolInfo*>::iterator it = avolInfoList.begin();
  while (it != avolInfoList.end()) {
    if (strcmp((*it)->compartmentId, cId) == 0) {
      return *it;
    }
    it++;
  }
  return NULL;

}

void countAST(ASTNode *ast, int &numOfASTNodes)
{
  for (unsigned int i = 0; i < ast->getNumChildren(); i++) {
    countAST(ast->getChild(i), numOfASTNodes);
  }
  numOfASTNodes++;
}

void rearrangeAST(ASTNode *ast)
{
  unsigned int i, j;
  ASTNodeType_t type = ast->getType();
  if (type == AST_MINUS && ast->getNumChildren() == 1) {//minus which has one child
    //"-a" to "-1.0 * a"
    ast->setType(AST_TIMES);
    ASTNode *ast_minus_one = new ASTNode(AST_REAL);
    ast_minus_one->setValue(-1.0);
    ast->addChild(ast_minus_one);
  } else if (type == AST_PLUS && ast->getNumChildren() == 1) {//plus which has one child
    //"+a" to "1.0 * a"
    ast->setType(AST_TIMES);
    ASTNode *ast_plus_one = new ASTNode(AST_REAL);
    ast_plus_one->setValue(1.0);
    ast->addChild(ast_plus_one);
  } else if (type == AST_FUNCTION_PIECEWISE) {//piecewise
    //remove all children
    vector<ASTNode*> astChildrenList;
    vector<ASTNode*> astBooleanList;
    unsigned int noc = ast->getNumChildren();
    ASTNode *astOtherwise = new ASTNode;

    //"piece" to boolean * expression
    for (j = 0; j < noc / 2; j++) {
      ASTNode *ast_times = new ASTNode(AST_TIMES);
      astBooleanList.push_back(ast->getChild(1)->deepCopy());
      ast_times->addChild(ast->getChild(0));
      ast_times->addChild(ast->getChild(1));
      astChildrenList.push_back(ast_times);
      ast->removeChild(0);
      ast->removeChild(0);
    }
    //"otherwise" to nand
    if (noc % 2 != 0) {
      astOtherwise->setType(AST_TIMES);
      ASTNode *otherwiseExpression = ast->getChild(0);
      ast->removeChild(0);
      ASTNode *ast_and = new ASTNode(AST_LOGICAL_AND);
      vector<ASTNode*>::iterator it = astBooleanList.begin();
      while (it != astBooleanList.end()) {
        ASTNode *ast_not = new ASTNode(AST_LOGICAL_NOT);
        ast_not->addChild(*it);
        ast_and->addChild(ast_not);
        it++;
      }
      astOtherwise->addChild(ast_and);
      astOtherwise->addChild(otherwiseExpression);
    } else {
      astOtherwise->setType(AST_INTEGER);
      astOtherwise->setValue(0);
    }
    ast->setType(AST_PLUS);
    ASTNode *ast_next = ast;
    ast_next->addChild(astOtherwise);

    for (j = 0; j < astChildrenList.size() - 1; j++) {
      ast_next->addChild(new ASTNode(AST_PLUS));
      ast_next = ast_next->getChild(1);
      ast_next->addChild(astChildrenList[j]);
    }
    ast_next->addChild(astChildrenList[astChildrenList.size() - 1]);
  } else if (ast->isLogical()) {//logical
    if (type != AST_LOGICAL_NOT) {//except for ast_logical_not
      if (ast->getNumChildren() == 1) {
        ASTNode *ast_one = new ASTNode(AST_INTEGER);
        ast_one->setValue(1);
        ast->addChild(ast_one);
      } else {
        ast->reduceToBinary();
      }
    } else {//not
      // 			vector<ASTNode*> astChildrenList;
      // 			unsigned int noc = ast->getNumChildren();
      // 			for (j = 0; j < noc; j++) {
      // 				astChildrenList.push_back(ast->getChild(0));
      // 				ast->removeChild(0);
      // 			}
      // 			ast->addChild(new ASTNode(AST_LOGICAL_AND));
      // 			vector<ASTNode*>::iterator it = astChildrenList.begin();
      // 			while (it != astChildrenList.end()) {
      // 				ast->getChild(0)->addChild(*it);
      // 				it++;
      // 			}
    }
  }
  for (i = 0; i < ast->getNumChildren(); i++) {
    rearrangeAST(ast->getChild(i));
  }
}

void reversePolishInitial(double **&varList, 
                          double **&constList, 
                          int *&opfuncList, 
                          double *&value, 
                          int numOfASTNodes, 
                          int Xindex, 
                          int Yindex, 
                          int Zindex, 
                          bool inVol)
{
  int X, Y, Z, i;
  int index = 0;
  double rpStack[stackMax] = {0};
  int increment = 0;
  if (inVol == true) {
    increment = 2;
  } else {
    increment = 1;
  }
  /*
  //volume
  for (Y = 0; Y < Yindex; Y += 2) {
  for (X = 0; X < Xindex; X += 2) {
  cout << (Y * Xindex + X) << " ";
  }
  }
  cout << endl;
  //membrane
  for (Y = 0; Y < Yindex; Y++) {
  for (X = 0; X < Xindex; X++) {
  if (!(X % 2 == 0 && Y % 2 == 0)) {
  cout << (Y * Xindex + X) << " ";
  }
  }
  }
  cout << endl;
  */

  for (Z = 0; Z < Zindex; Z += increment) {
    for (Y = 0; Y < Yindex; Y += increment) {
      for (X = 0; X < Xindex; X += increment) {
        if (inVol == true || (inVol == false && (Z * Yindex * Xindex + Y * Xindex + X) % 2 != 0)) {
          index = 0;
          for (i = 0; i < numOfASTNodes; i++) {
            if (varList[i] != 0) {//set variable into the stack
              rpStack[index] = varList[i][Z * Yindex * Xindex + Y * Xindex + X];
              index++;
            } else if (constList[i] != 0) {//set const into the stack
              rpStack[index] = *(constList[i]);
              index++;
            } else {//operation
              index--;
              switch (opfuncList[i]) {
              case AST_PLUS:
                rpStack[index - 1] += rpStack[index];
                break;
              case AST_MINUS:
                rpStack[index - 1] -= rpStack[index];
                break;
              case AST_TIMES:
                rpStack[index - 1] *= rpStack[index];
                break;
              case AST_DIVIDE:
                rpStack[index - 1] /= rpStack[index];
                break;
              case AST_POWER:
              case AST_FUNCTION_POWER:
                rpStack[index - 1] = pow(rpStack[index - 1], rpStack[index]);
                break;
              case AST_FUNCTION_ABS:
                rpStack[index] = fabs(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_ARCCOS:
                rpStack[index] = acos(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_ARCCOSH:
#if WIN32
                rpStack[index] = log(rpStack[index]-sqrt(pow(rpStack[index],2)-1));
#else
                rpStack[index] = acosh(rpStack[index]);
#endif
                index++;
                break;
              case AST_FUNCTION_ARCCSC:
                rpStack[index] = asin(1.0 / rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_ARCCSCH:
                break;
              case AST_FUNCTION_ARCSEC:
                break;
              case AST_FUNCTION_ARCSECH:
                break;
              case AST_FUNCTION_ARCSIN:
                rpStack[index] = asin(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_ARCSINH:
                break;
              case AST_FUNCTION_ARCTAN:
                rpStack[index] = atan(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_ARCTANH:
                break;
              case AST_FUNCTION_CEILING:
                break;
              case AST_FUNCTION_COS:
                rpStack[index] = cos(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_COSH:
                rpStack[index] = cosh(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_COT:
                break;
              case AST_FUNCTION_COTH:
                break;
              case AST_FUNCTION_CSC:
                rpStack[index] = 1.0 / sin(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_CSCH:
                break;
              case AST_FUNCTION_DELAY:
                break;
              case AST_FUNCTION_EXP:
                rpStack[index] = exp(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_FACTORIAL:
                break;
              case AST_FUNCTION_FLOOR:
                break;
              case AST_FUNCTION_LN:
                rpStack[index] = log(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_LOG:
                rpStack[index - 1] = log(rpStack[index]) / log(rpStack[index - 1]);
                break;
              case AST_FUNCTION_PIECEWISE:
                break;
              case AST_FUNCTION_ROOT:
                rpStack[index] = sqrt(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_SEC:
                break;
              case AST_FUNCTION_SECH:
                break;
              case AST_FUNCTION_SIN:
                rpStack[index] = sin(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_SINH:
                rpStack[index] = sinh(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_TAN:
                rpStack[index] = tan(rpStack[index]);
                index++;
                break;
              case AST_FUNCTION_TANH:
                rpStack[index] = tanh(rpStack[index]);
                index++;
                break;
              case AST_LAMBDA:
                break;
              case AST_LOGICAL_AND:
                rpStack[index - 1] = (static_cast<int>(rpStack[index - 1]) == 1 && static_cast<int>(rpStack[index] == 1)) ? 1.0: 0.0;
                break;
              case AST_LOGICAL_NOT:
                rpStack[index] = (static_cast<int>(rpStack[index]) == 0) ? 1.0: 0.0;
                index++;
                break;
              case AST_LOGICAL_OR:
                rpStack[index - 1] = (static_cast<int>(rpStack[index - 1]) == 1 || static_cast<int>(rpStack[index] == 1)) ? 1.0: 0.0;
                break;
              case AST_LOGICAL_XOR:
                rpStack[index - 1] = ((static_cast<int>(rpStack[index - 1]) == 1 && static_cast<int>(rpStack[index] == 0))
                  || (static_cast<int>(rpStack[index - 1]) == 0 && static_cast<int>(rpStack[index] == 1))) ? 1.0: 0.0;
                break;
              case AST_RATIONAL:
                break;
              case AST_RELATIONAL_EQ:
                rpStack[index - 1] = (rpStack[index - 1] == rpStack[index]) ? 1.0: 0.0;
                break;
              case AST_RELATIONAL_GEQ:
                rpStack[index - 1] = (rpStack[index - 1] >= rpStack[index]) ? 1.0: 0.0;
                break;
              case AST_RELATIONAL_GT:
                rpStack[index - 1] = (rpStack[index - 1] > rpStack[index]) ? 1.0: 0.0;
                break;
              case AST_RELATIONAL_LEQ:
                rpStack[index - 1] = (rpStack[index - 1] <= rpStack[index]) ? 1.0: 0.0;
                break;
              case AST_RELATIONAL_LT:
                rpStack[index - 1] = (rpStack[index - 1] < rpStack[index]) ? 1.0: 0.0;
                break;
              case AST_RELATIONAL_NEQ:
                rpStack[index - 1] = (rpStack[index - 1] != rpStack[index]) ? 1.0: 0.0;
                break;
              default:
                break;
              }
            }
          }
          index--;
          value[Z * Yindex * Xindex + Y * Xindex + X] = rpStack[index];
        }
      }
    }
  }
}

void freeVarInfo(vector<variableInfo*> &varInfoList, double *sim_time_p)
{
  for (size_t i = 0; i < varInfoList.size(); i++) {
    variableInfo *info = varInfoList[i];
    if (strcmp(info->id, "t") != 0) {
      //value
      delete[] info->value;
      info->value = 0;
      //delta
      delete[] info->delta;
      info->delta = 0;
      //rpinfo
      if (info->rpInfo != 0) {
        //varList
        delete[] info->rpInfo->varList;
        info->rpInfo->varList = 0;
        //deltaList
        delete[] info->rpInfo->deltaList;
        info->rpInfo->deltaList = 0;
        //constList
        if (info->rpInfo->constList != 0) {//constList
          for (int j = 0; j < info->rpInfo->listNum; j++) {
            if (info->rpInfo->constList[j] != sim_time_p) {
              delete info->rpInfo->constList[j];
            }
          }
        }
        delete[] info->rpInfo->constList;
        info->rpInfo->constList = 0;
        //opfuncList
        delete[] info->rpInfo->opfuncList;
        info->rpInfo->opfuncList = 0;
        //rpInfo
        delete info->rpInfo;
        info->rpInfo = 0;
      }
      //diffCInfo
      if (info->dci != 0) {
        for (int j = 0; j < 3; j++) {
          delete info->dci[j];
          info->dci[j] = NULL;
        }
        delete[] info->dci;
        info->dci = 0;
      }
      //adCInfo
      if (info->aci != 0) {
        for (int j = 0; j < 3; j++) {
          delete info->aci[j];
        }
        delete[] info->aci;
        info->aci = 0;
      }
      //derivativeCIP
      if (info->dcip != 0) {
        delete[] info->dcip->fx;
        info->dcip->fx = 0;
        delete[] info->dcip->fx_next;
        info->dcip->fx_next = 0;
        delete[] info->dcip->fy;
        info->dcip->fy = 0;
        delete[] info->dcip->fy_next;
        info->dcip->fy_next = 0;
        delete[] info->dcip->fz;
        info->dcip->fz_next = 0;
        delete[] info->dcip->fz_next;
        info->dcip->fz_next = 0;
        delete info->dcip;
        info->dcip = 0;
      }
      //info
      delete info;
      info = 0;
    } else {
      // const char* does not need to be freed!
      //free(info->id);
      //info->id = 0;
      delete info;
      info = 0;
    }
  }
}

void freeBcOfSpeciesInfo(vector<bcOfSpeciesInfo*> &bcOfSpeciesInfoList)
{
  for (size_t i = 0; i < bcOfSpeciesInfoList.size(); i++) {
    bcOfSpeciesInfo *bcos = bcOfSpeciesInfoList[i];
    delete bcos;
    bcos = 0;
  }
}

void freeAvolInfo(vector<analyticVolInfo*> &avolInfoList)
{
  for (size_t i = 0; i < avolInfoList.size(); i++) {
    analyticVolInfo *avolInfo = avolInfoList[i];
    if (avolInfo->rpInfo != 0) {
      //varList
      delete[] avolInfo->rpInfo->varList;
      avolInfo->rpInfo->varList = 0;
      //deltaList
      delete[] avolInfo->rpInfo->deltaList;
      avolInfo->rpInfo->deltaList = 0;
      //constList
      if (avolInfo->rpInfo->constList != 0) {//constList
        for (int j = 0; j < avolInfo->rpInfo->listNum; j++) {
          delete avolInfo->rpInfo->constList[j];
        }
      }
      delete[] avolInfo->rpInfo->constList;
      avolInfo->rpInfo->constList = 0;
      //opfuncList
      delete[] avolInfo->rpInfo->opfuncList;
      avolInfo->rpInfo->opfuncList = 0;
      //rpInfo
      delete avolInfo->rpInfo;
      avolInfo->rpInfo = 0;
    }
    //isDomain
    delete[] avolInfo->isDomain;
    avolInfo->isDomain = 0;
    //isBoundary
    delete[] avolInfo->isBoundary;
    avolInfo->isBoundary = 0;
    //isEdge
    // 		delete avolInfo->isEdge;
    // 		avolInfo->isEdge = 0;
    //boundaryType
    delete[] avolInfo->bt;
    avolInfo->bt = 0;
    //avolinfo
    delete avolInfo;
    avolInfo = 0;
  }
}

bool isNormalParameter(SpatialParameterPlugin *pPlugin ) 
{
  if (pPlugin == NULL) return true;

  if (pPlugin ->getSpatialSymbolReference() != NULL && !pPlugin ->getSpatialSymbolReference()->getSpatialId().empty())
    return false;

  if (pPlugin ->getDiffusionCoefficient() != NULL && !pPlugin ->getDiffusionCoefficient()->getVariable().empty())
    return false;

  if (pPlugin ->getBoundaryCondition() != NULL && !pPlugin ->getBoundaryCondition()->getVariable().empty())
    return false;

  if (pPlugin ->getAdvectionCoefficient() != NULL && !pPlugin ->getAdvectionCoefficient()->getVariable().empty())
    return false;

  return true;
}

void freeRInfo(vector<reactionInfo*> &rInfoList)
{
  for (size_t i = 0; i < rInfoList.size(); i++) {
    reactionInfo *rInfo = rInfoList[i];
    //value
    delete[] rInfo->value;
    rInfo->value = 0;
    if (rInfo->rpInfo != 0) {
      //varList
      delete[] rInfo->rpInfo->varList;
      rInfo->rpInfo->varList = 0;
      //deltaList
      delete[] rInfo->rpInfo->deltaList;
      rInfo->rpInfo->deltaList = 0;
      //constList
      if (rInfo->rpInfo->constList != 0) {//constList
        for (int j = 0; j < rInfo->rpInfo->listNum; j++) {
          delete rInfo->rpInfo->constList[j];
        }
      }
      delete[] rInfo->rpInfo->constList;
      rInfo->rpInfo->constList = 0;
      //opfuncList
      delete[] rInfo->rpInfo->opfuncList;
      rInfo->rpInfo->opfuncList = 0;
      //rpInfo
      delete rInfo->rpInfo;
      rInfo->rpInfo = 0;
    }
    //rInfo
    delete rInfo;
    rInfo = 0;
  }
}

int SpatialSimulator::getVariableLength() const
{
return (Zindex -1) * Yindex * Xindex + (Yindex-1) * Xindex + Xindex-1;
}

void reversePolishRK(reactionInfo *&rInfo, variableInfo *&sInfo, int Xindex, int Yindex, int Zindex, double dt, int m, materialType mType, bool inVol, diffCInfo **dci)
{
  int X, Y, Z, i;
  int st_index = 0, index = 0, numOfVolIndexes = Xindex * Yindex * Zindex;
  double rpStack[stackMax] = {0};
  double rk[4] = {0, 0.5, 0.5, 1.0};
  double **variable = rInfo->rpInfo->varList;
  double **constant = rInfo->rpInfo->constList;
  double **d = rInfo->rpInfo->deltaList;
  int *operation = rInfo->rpInfo->opfuncList;
  int numOfASTNodes = rInfo->rpInfo->listNum;
  analyticVolInfo *avolInfo = sInfo->avi;
  int increment = 0;
  if (inVol == true) {
    increment = 2;
  } else {
    increment = 1;
  }

  for (Z = 0; Z < Zindex; Z += increment) {
    for (Y = 0; Y < Yindex; Y += increment) {
      for (X = 0; X < Xindex; X += increment) {
        if (inVol == true || (inVol == false && (Z * Yindex * Xindex + Y * Xindex + X) % 2 != 0)) {
          st_index = 0;
          index = Z * Yindex * Xindex + Y * Xindex + X;
          //if (avolInfo->isDomain[index] == 1 && ((avolInfo->isBoundary[index] == 0 && sInfo->dci != 0) || sInfo->dci == 0)) {
          //if (avolInfo->isDomain[index] == 1) {
          if (avolInfo->isDomain[index] == 1 && (dci == 0 || (dci != 0 && avolInfo->isBoundary[index] == 0))) {
            for (i = 0; i < numOfASTNodes; i++) {
              if (variable[i] != 0) {//set variable into the stack
                //rpStack[st_index] = variable[i][index];
                if (d != 0 && d[i] != 0) {
                  if (m == 0) {
                    rpStack[st_index] = variable[i][index];
                  } else {
                    //rpStack[st_index] = variable[i][index];
                    rpStack[st_index] = variable[i][index] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + index];
                  }
                } else {
                  rpStack[st_index] = variable[i][index];
                }
                st_index++;
              } else if (constant[i] != 0) {//set const into the stack
                rpStack[st_index] = *(constant[i]);
                st_index++;
              } else {//operation
                st_index--;
                switch (operation[i]) {
                case AST_PLUS:
                  rpStack[st_index - 1] += rpStack[st_index];
                  break;
                case AST_MINUS:
                  rpStack[st_index - 1] -= rpStack[st_index];
                  break;
                case AST_TIMES:
                  rpStack[st_index - 1] *= rpStack[st_index];
                  break;
                case AST_DIVIDE:
                  rpStack[st_index - 1] /= rpStack[st_index];
                  break;
                case AST_POWER:
                case AST_FUNCTION_POWER:
                  rpStack[st_index - 1] = pow(rpStack[st_index - 1], rpStack[st_index]);
                  break;
                case AST_FUNCTION_ABS:
                  rpStack[st_index] = fabs(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_ARCCOS:
                  rpStack[st_index] = acos(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_ARCCOSH:
#if WIN32
                  rpStack[index] = log(rpStack[index]-sqrt(pow(rpStack[index],2)-1));
#else
                  rpStack[index] = acosh(rpStack[index]);
#endif
                  st_index++;
                  break;
                case AST_FUNCTION_ARCCSC:
                  rpStack[st_index] = asin(1.0 / rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_ARCCSCH:
                  break;
                case AST_FUNCTION_ARCSEC:
                  break;
                case AST_FUNCTION_ARCSECH:
                  break;
                case AST_FUNCTION_ARCSIN:
                  rpStack[st_index] = asin(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_ARCSINH:
                  break;
                case AST_FUNCTION_ARCTAN:
                  rpStack[st_index] = atan(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_ARCTANH:
                  break;
                case AST_FUNCTION_CEILING:
                  break;
                case AST_FUNCTION_COS:
                  rpStack[st_index] = cos(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_COSH:
                  rpStack[st_index] = cosh(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_COT:
                  break;
                case AST_FUNCTION_COTH:
                  break;
                case AST_FUNCTION_CSC:
                  rpStack[st_index] = 1.0 / sin(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_CSCH:
                  break;
                case AST_FUNCTION_DELAY:
                  break;
                case AST_FUNCTION_EXP:
                  rpStack[st_index] = exp(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_FACTORIAL:
                  break;
                case AST_FUNCTION_FLOOR:
                  break;
                case AST_FUNCTION_LN:
                  rpStack[st_index] = log(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_LOG:
                  rpStack[st_index - 1] = log(rpStack[st_index]) / log(rpStack[st_index - 1]);
                  break;
                case AST_FUNCTION_PIECEWISE:
                  break;
                case AST_FUNCTION_ROOT:
                  rpStack[st_index] = sqrt(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_SEC:
                  break;
                case AST_FUNCTION_SECH:
                  break;
                case AST_FUNCTION_SIN:
                  rpStack[st_index] = sin(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_SINH:
                  rpStack[st_index] = sinh(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_TAN:
                  rpStack[st_index] = tan(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_FUNCTION_TANH:
                  rpStack[st_index] = tanh(rpStack[st_index]);
                  st_index++;
                  break;
                case AST_LAMBDA:
                  break;
                case AST_LOGICAL_AND:
                  rpStack[st_index - 1] = (static_cast<int>(rpStack[st_index - 1]) == 1 && static_cast<int>(rpStack[st_index] == 1)) ? 1.0: 0.0;
                  break;
                case AST_LOGICAL_NOT:
                  rpStack[st_index] = (static_cast<int>(rpStack[st_index]) == 0) ? 1.0: 0.0;
                  st_index++;
                  break;
                case AST_LOGICAL_OR:
                  rpStack[st_index - 1] = (static_cast<int>(rpStack[st_index - 1]) == 1 || static_cast<int>(rpStack[st_index] == 1)) ? 1.0: 0.0;
                  break;
                case AST_LOGICAL_XOR:
                  rpStack[st_index - 1] = ((static_cast<int>(rpStack[st_index - 1]) == 1 && static_cast<int>(rpStack[st_index] == 0))
                    || (static_cast<int>(rpStack[st_index - 1]) == 0 && static_cast<int>(rpStack[st_index] == 1))) ? 1.0: 0.0;
                  break;
                case AST_RATIONAL:
                  break;
                case AST_RELATIONAL_EQ:
                  rpStack[st_index - 1] = (rpStack[st_index - 1] == rpStack[st_index]) ? 1.0: 0.0;
                  break;
                case AST_RELATIONAL_GEQ:
                  rpStack[st_index - 1] = (rpStack[st_index - 1] >= rpStack[st_index]) ? 1.0: 0.0;
                  break;
                case AST_RELATIONAL_GT:
                  rpStack[st_index - 1] = (rpStack[st_index - 1] > rpStack[st_index]) ? 1.0: 0.0;
                  break;
                case AST_RELATIONAL_LEQ:
                  rpStack[st_index - 1] = (rpStack[st_index - 1] <= rpStack[st_index]) ? 1.0: 0.0;
                  break;
                case AST_RELATIONAL_LT:
                  rpStack[st_index - 1] = (rpStack[st_index - 1] < rpStack[st_index]) ? 1.0: 0.0;
                  break;
                case AST_RELATIONAL_NEQ:
                  rpStack[st_index - 1] = (rpStack[st_index - 1] != rpStack[st_index]) ? 1.0: 0.0;
                  break;
                default:
                  break;
                }
              }
            }
            st_index--;
            if (mType == products) {//products
              sInfo->delta[m * numOfVolIndexes + index] += rpStack[st_index];
            } else if (mType == reactants) {//reactants
              sInfo->delta[m * numOfVolIndexes + index] -= rpStack[st_index];
            }
          }
        }
      }
    }
  }
}

void parseDependence(const ASTNode *ast, vector<variableInfo*> &dependence, vector<variableInfo*> &varInfoList)
{
  unsigned int i;
  for (i = 0; i < ast->getNumChildren(); i++) {
    parseDependence(ast->getChild(i), dependence, varInfoList);
  }
  if (ast->isName()) {
    if (searchInfoById(dependence, ast->getName()) == 0) {
      dependence.push_back(searchInfoById(varInfoList, ast->getName()));
    }
  }
}

bool isResolvedAll(vector<variableInfo*> &dependence)
{
  vector<variableInfo*>::iterator it = dependence.begin();
  while (it != dependence.end()) {
    if (!(*it)->isResolved) {
      return false;
    }
    it++;
  }
  return true;
}

int SpatialSimulator::runOldMain(int argc, const char* argv[])
{
  return old_main(argc, argv);
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

    if (pPlugin->getDiffusionCoefficient() != 0) {
      variableInfo *sInfo = searchInfoById(varInfoList, pPlugin->getDiffusionCoefficient()->getVariable().c_str());
      if (sInfo != 0) {
        DiffusionCoefficient *dc = pPlugin->getDiffusionCoefficient();
        for (int Z = 0; Z < Zindex; Z++) {
          for (int Y = 0; Y < Yindex; Y++) {
            for (int X = 0; X < Xindex; X++) {
              sInfo->dci[dc->getCoordinateIndex()]->value[Z * Yindex * Xindex + Y * Xindex + X] = value;
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

  for (int Z = 0; Z < Zindex; Z++) {
#pragma omp parallel for
    for (int Y = 0; Y < Yindex; Y++) {
      for (int X = 0; X < Xindex; X++) {
        info->value[Z * Yindex * Xindex + Y * Xindex + X] = value;
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
    if (coefficient == NULL || coefficient->getVariable() != id || coefficient->getCoordinateIndex() != index)
      continue;    
    return current;
  }
  return NULL;
}

  SpatialSimulator::SpatialSimulator():
  Xdiv(100), Ydiv(100), Zdiv(1), model(NULL), volDimension(0), memDimension(0)
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

  XMLNamespaces *xns         = doc->getNamespaces();
  string spatialPrefix       = xns->getPrefix("http://www.sbml.org/sbml/level3/version1/spatial/version1");
  string reqPrefix           = xns->getPrefix("http://www.sbml.org/sbml/level3/version1/requiredElements/version1");
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
  bcOfSpeciesInfoList = vector<bcOfSpeciesInfo*>();
  avolInfoList        = vector<analyticVolInfo*>();
  bcInfoList            = vector<boundaryCInfo*>();
  rInfoList              = vector<reactionInfo*>();
  acInfoList                  = vector<adCInfo*>();
  memList                  = vector<const char*>();
  dimension                       = geometry->getNumCoordinateComponents();

  int Xplus = 0, Xminus = 0, Yplus = 0, Yminus = 0, Zplus = 0, Zminus = 0;


  Xindex = 2 * Xdiv - 1;
  Yindex = 2 * Ydiv - 1;
  Zindex = 2 * Zdiv - 1;

  numOfVolIndexes = Xindex * Yindex * Zindex;
  //unit

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

#pragma region // read compartment
  //set id and value
  //compartment
  for (i = 0; i < numOfCompartments; i++) {
    Compartment *c = loc->get(i);
    variableInfo *info = new variableInfo;
    InitializeVarInfo(info);
    varInfoList.push_back(info);
    info->com = c;
    info->id = c->getId().c_str();
    info->value = new double[numOfVolIndexes];
    INIT_DOUBLE(info->value, numOfVolIndexes);
    if (c->isSetSize()) {
#ifdef PRINT_DEBUG
      cout << c->getId() << endl;
#endif
      info->isResolved = true;
      if (c->getSpatialDimensions() == volDimension) {//volume
        for (Z = 0; Z < Zindex; Z += 2) {
          for (Y = 0; Y < Yindex; Y += 2) {
            for (X = 0; X < Xindex; X += 2) {
              info->value[Z * Yindex * Xindex + Y * Xindex + X] = c->getSize();
            }
          }
        }
      } else if (c->getSpatialDimensions() == memDimension) {//membrane
        for (Z = 0; Z < Zindex; Z++) {
          for (Y = 0; Y < Yindex; Y++) {
            for (X = 0; X < Xindex; X++) {
              if (!(X % 2 == 0 && Y % 2 == 0 && Z % 2 == 0)) {
                info->value[Z * Yindex * Xindex + Y * Xindex + X] = c->getSize();
              }
            }
          }
        }
      }
    }
  }
#pragma endregion

#pragma region // read species
  //species
  for (i = 0; i < numOfSpecies; i++) {
    Species *s = los->get(i);
    RequiredElementsSBasePlugin* reqplugin = static_cast<RequiredElementsSBasePlugin*>(s->getPlugin("req"));
    SpatialSpeciesRxnPlugin* splugin = static_cast<SpatialSpeciesRxnPlugin*>(s->getPlugin("spatial"));
    //species have spatial extension
    if (reqplugin->getMathOverridden() == spatialPrefix) {
      bcOfSpeciesInfo *bcos = new bcOfSpeciesInfo;
      bcos->speciesId = s->getId().c_str();
      InitializeBcOfSpeciesInfo(bcos);
      bcOfSpeciesInfoList.push_back(bcos);
      variableInfo *info = new variableInfo;
      InitializeVarInfo(info);
      varInfoList.push_back(info);
      info->sp = s;
      info->com = model->getCompartment(s->getCompartment());
      info->id = s->getId().c_str();
      if (info->com->getSpatialDimensions() == volDimension) {
        info->inVol = true;
      } else if (info->com->getSpatialDimensions() == memDimension) {
        info->inVol = false;
      }
      //species value is specified by initial amount, initial value, rule or initial assignment
      //species is spatially defined
      if (reqplugin->getCoreHasAlternateMath() && splugin->getIsSpatial()) {
        if (s->isSetInitialAmount() || s->isSetInitialConcentration()) {//Initial Amount or Initial Concentration
          info->value = new double[numOfVolIndexes];
          INIT_DOUBLE(info->value, numOfVolIndexes);
          if (!s->isSetConstant() || !s->getConstant()) {
            info->delta = new double[4 * numOfVolIndexes];
            INIT_DOUBLE(info->delta, 4*numOfVolIndexes);
          }
          if (s->isSetInitialAmount()) {//Initial Amount
            info->isResolved = true;
            for (Z = 0; Z < Zindex; Z++) {
              for (Y = 0; Y < Yindex; Y++) {
                for (X = 0; X < Xindex; X++) {
                  info->value[Z * Yindex * Xindex + Y * Xindex + X] = s->getInitialAmount();
                }
              }
            }
          } else if (s->isSetInitialConcentration()) {//Initial Concentration
            for (Z = 0; Z < Zindex; Z++) {
              for (Y = 0; Y < Yindex; Y++) {
                for (X = 0; X < Xindex; X++) {
                  info->value[Z * Yindex * Xindex + Y * Xindex + X] = s->getInitialConcentration();
                }
              }
            }
          }
        }
      }
    }
  }
#pragma endregion


  
  for (i = 0; i < varInfoList.size(); i++)
  {
    variableInfo* current = varInfoList[i];
    Parameter* diffX = getDiffusionCoefficientForSpecies(current->id);
    if (diffX == NULL)
      continue;
    Parameter* diffY = getDiffusionCoefficientForSpecies(current->id, 1);
    if (diffY != NULL)
      continue;

    diffY = diffX->clone();
    diffY->setId(diffY->getId() + "_1");
    string id = diffY->getId();
    model->addParameter(diffY);
    diffY = model->getParameter(id);
    SpatialParameterPlugin *plugin = (SpatialParameterPlugin*)diffY->getPlugin("spatial");
    plugin->getDiffusionCoefficient()->setCoordinateIndex(1);
   
  }
  
  numOfParameters   = static_cast<unsigned int>(model->getNumParameters());  
  lop      = model->getListOfParameters();

#pragma region // read Parameter
   
  //parameter
  for (i = 0; i < numOfParameters; i++) {
    Parameter *p = lop->get(i);
    SpatialParameterPlugin *pPlugin = static_cast<SpatialParameterPlugin*>(p->getPlugin("spatial"));
    variableInfo *info = new variableInfo;
    InitializeVarInfo(info);
    varInfoList.push_back(info);
    info->para = p;
    info->id = p->getId().c_str();
    info->value = new double[numOfVolIndexes];
    INIT_DOUBLE(info->value, numOfVolIndexes);
    if (isNormalParameter(pPlugin ) ) {//normal parameter
      if (p->isSetValue()) {
        info->isResolved = true;
        setParameterUniformly(info, p->getValue());        
      }
    } else {//spatial parameter plugin
      //なぜかnullにならない
      // 			//normal
      // 			if (p->isSetValue()) {
      // 				info->isResolved = true;
      // 				for (Z = 0; Z < Zindex; Z++) {
      // 					for (Y = 0; Y < Yindex; Y++) {
      // 						for (X = 0; X < Xindex; X++) {
      // 							info->value[Z * Yindex * Xindex + Y * Xindex + X] = p->getValue();
      // 						}
      // 					}
      // 				}
      // 			}
      //diffusion coefficient
      if (pPlugin->getDiffusionCoefficient() != 0) {
        variableInfo *sInfo = searchInfoById(varInfoList, pPlugin->getDiffusionCoefficient()->getVariable().c_str());
        if (sInfo != 0) {
          DiffusionCoefficient *dc = pPlugin->getDiffusionCoefficient();
          if (sInfo->dci == 0) {
            sInfo->dci = new diffCInfo*[3];
            for (j = 0; j < 3; j++) {
              sInfo->dci[j] = 0;
            }
          }

          sInfo->dci[dc->getCoordinateIndex()] = new diffCInfo;
          sInfo->dci[dc->getCoordinateIndex()]->diffc = dc;
          sInfo->dci[dc->getCoordinateIndex()]->value = new double[numOfVolIndexes];
          for (Z = 0; Z < Zindex; Z++) {
              for (Y = 0; Y < Yindex; Y++) {
                for (X = 0; X < Xindex; X++) {
                  sInfo->dci[dc->getCoordinateIndex()]->value[Z * Yindex * Xindex + Y * Xindex + X]  =  p->getValue();
                }
              }
            }
          //cout << sInfo->id << " " << sInfo->dci[dc->getCoordinateIndex()]->value << " " << dc->getCoordinateIndex() << " " << sInfo->dci << endl;
          sInfo->dci[dc->getCoordinateIndex()]->rpInfo = 0;
        }
      }
      //advection coefficient
      if (pPlugin->getAdvectionCoefficient() != 0) {
        variableInfo *sInfo = searchInfoById(varInfoList, pPlugin->getAdvectionCoefficient()->getVariable().c_str());
        if (sInfo != 0) {
          AdvectionCoefficient *ac = pPlugin->getAdvectionCoefficient();
          if (sInfo->aci == 0) {
            sInfo->aci = new adCInfo*[3];
            for (j = 0; j < 3; j++) {
              sInfo->aci[j] = 0;
            }
          }
          sInfo->aci[ac->getCoordinateIndex()] = new adCInfo;
          sInfo->aci[ac->getCoordinateIndex()]->adc = ac;
          sInfo->aci[ac->getCoordinateIndex()]->value = new double[numOfVolIndexes];
                    for (Z = 0; Z < Zindex; Z++) {
              for (Y = 0; Y < Yindex; Y++) {
                for (X = 0; X < Xindex; X++) {
                  sInfo->aci[ac->getCoordinateIndex()]->value[Z * Yindex * Xindex + Y * Xindex + X]  =  p->getValue();
                }
              }
            }

          sInfo->aci[ac->getCoordinateIndex()]->rpInfo = 0;
        }
      }
      //get coordinate information
      //Cartesian
      if (pPlugin->getSpatialSymbolReference() != 0) {
        if (pPlugin->getSpatialSymbolReference()->getType() == "coordinateComponent") {
          CoordinateComponent *cc = geometry->getCoordinateComponent(p->getId());
          double min = cc->getBoundaryMin()->getValue();
          double max = cc->getBoundaryMax()->getValue();
          if (cc->getComponentType() == "cartesianX") {
            xaxis = const_cast<char*>(p->getId().c_str());
            deltaX = (max - min) / (Xdiv - 1);
            info->isResolved = true;
            for (Z = 0; Z < Zindex; Z++) {
              for (Y = 0; Y < Yindex; Y++) {
                for (X = 0; X < Xindex; X++) {
                  info->value[Z * Yindex * Xindex + Y * Xindex + X] = min + (double)X * deltaX / 2.0;
                }
              }
            }
          } else if (cc->getComponentType() == "cartesianY") {
            yaxis = const_cast<char*>(p->getId().c_str());
            deltaY = (max - min) / (Ydiv - 1);
            info->isResolved = true;
            for (Z = 0; Z < Zindex; Z++) {
              for (Y = 0; Y < Yindex; Y++) {
                for (X = 0; X < Xindex; X++) {
                  info->value[Z * Yindex * Xindex + Y * Xindex + X] = min + (double)Y * deltaY / 2.0;
                }
              }
            }
          } else if (cc->getComponentType() == "cartesianZ") {
            zaxis = const_cast<char*>(p->getId().c_str());
            deltaZ = (max - min) / (Zdiv - 1);
            info->isResolved = true;
            for (Z = 0; Z < Zindex; Z++) {
              for (Y = 0; Y < Yindex; Y++) {
                for (X = 0; X < Xindex; X++) {
                  info->value[Z * Yindex * Xindex + Y * Xindex + X] = min + (double)Z * deltaZ / 2.0;
                }
              }
            }
          }
        }
      }
      //boundary condition
      if (pPlugin->getBoundaryCondition() != 0) {
        bcOfSpeciesInfo *bcos = searchBCOSInfoBySpeId(bcOfSpeciesInfoList, pPlugin->getBoundaryCondition()->getVariable().c_str());
        BoundaryCondition *bcon = pPlugin->getBoundaryCondition();
        if (bcos != 0 && bcon != 0) {
          for (Z = 0; Z < Zindex; Z += 2) {
            for (Y = 0; Y < Yindex; Y += 2) {
              for (X = 0; X < Xindex; X += 2) {
                info->value[Z * Yindex * Xindex + Y * Xindex + X] = p->getValue();
              }
            }
          }
          if (bcon->getCoordinateBoundary() == "Xmax") {
            bcos->bcXp = info;
            bcos->bcXpType = bcon->getType();
          }
          if (bcon->getCoordinateBoundary() == "Xmin") {
            bcos->bcXm = info;
            bcos->bcXmType = bcon->getType();
          }
          if (bcon->getCoordinateBoundary() == "Ymax") {
            bcos->bcYp = info;
            bcos->bcYpType = bcon->getType();
          }
          if (bcon->getCoordinateBoundary() == "Ymin") {
            bcos->bcYm = info;
            bcos->bcYmType = bcon->getType();
          }
          if (bcon->getCoordinateBoundary() == "Zmax") {
            bcos->bcZp = info;
            bcos->bcZpType = bcon->getType();
          }
          if (bcon->getCoordinateBoundary() == "Zmin") {
            bcos->bcZm = info;
            bcos->bcZmType = bcon->getType();
          }
        }
      }
    }
  }
#pragma endregion



  //time
  t_info = new variableInfo;
  InitializeVarInfo(t_info);
  varInfoList.push_back(t_info);
  t_info->id = (const char*)malloc(sizeof(char) * 1 + 1);
  strcpy(const_cast<char*>(t_info->id), "t");
  t_info->value = &sim_time;
  t_info->isResolved = true;

#pragma region // read geometry definitions
  //geometryDefinition
  for (i = 0; i < geometry->getNumGeometryDefinitions(); i++) {
    if (geometry->getGeometryDefinition(i)->isAnalyticGeometry()) {
      //AnalyticVolumes
      AnalyticGeometry *analyticGeo = static_cast<AnalyticGeometry*>(geometry->getGeometryDefinition(i));
      //gather information of compartment, domainType, analyticVolume
      for (j = 0; j < numOfCompartments; j++) {
        Compartment *c = loc->get(j);
        if (SpatialCompartmentPlugin *cPlugin = static_cast<SpatialCompartmentPlugin*>(c->getPlugin("spatial"))) 
        {
          if (AnalyticVolume *analyticVol = analyticGeo->getAnalyticVolume(cPlugin->getCompartmentMapping()->getDomainType())) 
          {
            analyticVolInfo *avolInfo = new analyticVolInfo;
            InitializeAVolInfo(avolInfo);
            avolInfo->compartmentId = c->getId().c_str();
            avolInfo->domainTypeId = cPlugin->getCompartmentMapping()->getDomainType().c_str();
            avolInfo->domainId = 0;
            avolInfo->bt = new boundaryType[numOfVolIndexes];
            avolInfo->isVol = true;
            ast = const_cast<ASTNode*>(analyticVol->getMath());
#ifdef PRINT_DEBUG
            cout << "before: " << SBML_formulaToString(ast) << endl;
#endif
            rearrangeAST(ast);
#ifdef PRINT_DEBUG
            cout << "after: " << SBML_formulaToString(ast) << endl;
#endif
            numOfASTNodes = 0;
            countAST(ast, numOfASTNodes);
            avolInfo->rpInfo = new reversePolishInfo();
            avolInfo->rpInfo->varList = new double*[numOfASTNodes];
            avolInfo->rpInfo->constList = new double*[numOfASTNodes];
            avolInfo->rpInfo->opfuncList = new int[numOfASTNodes];
            avolInfo->isDomain = new double[numOfVolIndexes];
            INIT_DOUBLE(avolInfo->isDomain, numOfVolIndexes);
            avolInfo->isBoundary = new int[numOfVolIndexes];
            INIT_INT(avolInfo->isBoundary, numOfVolIndexes);
            avolInfo->rpInfo->listNum = numOfASTNodes;
            avolInfo->adjacent0 = 0;
            avolInfo->adjacent1 = 0;
            parseAST(ast, avolInfo->rpInfo->varList, avolInfo->rpInfo->deltaList, avolInfo->rpInfo->constList, avolInfo->rpInfo->opfuncList, varInfoList, numOfASTNodes);
            /*
            for (k = 0; k < geometry->getNumAdjacentDomains(); k++) {//AdjacentDomain
            AdjacentDomains *adDomain = geometry->getAdjacentDomains(k);
            if (adDomain->getDomain1() == domain->getSpatialId() || adDomain->getDomain2() == domain->getSpatialId()) {
            if (avolInfo->adjacent0 == 0) {
            avolInfo->adjacent0 = adDomain;
            } else if (avolInfo->adjacent1 == 0) {
            avolInfo->adjacent1 = adDomain;
            }
            }
            }
            */
            //parseAST(ast, avolInfo->rpInfo->varList, avolInfo->rpInfo->constList, avolInfo->rpInfo->opfuncList, varInfoList, numOfASTNodes);
            //judge if the coordinate point is inside the analytic volume
            reversePolishInitial(
              avolInfo->rpInfo->varList, 
              avolInfo->rpInfo->constList, 
              avolInfo->rpInfo->opfuncList, 
              avolInfo->isDomain, 
              numOfASTNodes, 
              Xindex, 
              Yindex, 
              Zindex, true);
            avolInfoList.push_back(avolInfo);
          }
        }
      }
    } else if (geometry->getGeometryDefinition(i)->isSampledFieldGeometry()) {
      //SampleFieldGeometry
    } else if (geometry->getGeometryDefinition(i)->isCSGeometry()){
      //CSGeometry
    } else if (geometry->getGeometryDefinition(i)->isParametricGeometry()) {
      //ParametricGeometry
    }
  }
#pragma endregion

#pragma region // merge geometries
  //merge external and internal analytic volumes and get boundary points of the geometry
  for (i = 0; i < geometry->getNumGeometryDefinitions(); i++) {
    AnalyticGeometry *analyticGeo = static_cast<AnalyticGeometry*>(geometry->getGeometryDefinition(i));
    if (geometry->getGeometryDefinition(i)->isAnalyticGeometry()) {
      for (j = 0; j < analyticGeo->getNumAnalyticVolumes(); j++) {
        AnalyticVolume *analyticVolEx = analyticGeo->getAnalyticVolume(j);
        analyticVolInfo *avolInfoEx = searchAvolInfoByDomainType(avolInfoList, analyticVolEx->getDomainType().c_str());
        if (avolInfoEx == NULL)
          continue;
        for (k = 0; k < analyticGeo->getNumAnalyticVolumes(); k++) {
          AnalyticVolume *analyticVolIn = analyticGeo->getAnalyticVolume(k);
          if (analyticVolEx->getOrdinal() > analyticVolIn->getOrdinal()) {//ex's ordinal is larger than in's ordinal
            analyticVolInfo *avolInfoIn = searchAvolInfoByDomainType(avolInfoList, analyticVolIn->getDomainType().c_str());
            if (avolInfoIn != NULL)
            {
            //merge
            for (Z = 0; Z < Zindex; Z += 2) {
              for (Y = 0; Y < Yindex; Y += 2) {
                for (X = 0; X < Xindex; X += 2) {
                  index = Z * Yindex * Xindex + Y * Xindex + X;
                  //if external == true && internal == true, then external = false;
                  if (avolInfoEx->isDomain[index] == 1 && avolInfoIn->isDomain[index] == 1) {
                    avolInfoEx->isDomain[index] = 0;
                  }
                }
              }
            }
            }
          }
        }
        //boundary
        switch(dimension) {
        case 1://1D
          for (X = 0; X < Xindex; X++) {
            avolInfoEx->isBoundary[X] = 0;//initialize
            if (static_cast<int>(avolInfoEx->isDomain[X]) == 1) {
              if (X == 0 || X == Xindex - 1) {//bounary of the domain
                avolInfoEx->isBoundary[X] = 1;
              } else if (static_cast<int>(avolInfoEx->isDomain[X + 2]) == 0
                || static_cast<int>(avolInfoEx->isDomain[X - 2]) == 0) {
                  avolInfoEx->isBoundary[Y * Xindex + X] = 1;
              }
            }
          }
          break;
        case 2://2D
          for (Y = 0; Y < Yindex; Y += 2) {
            for (X = 0; X < Xindex; X += 2) {
              avolInfoEx->isBoundary[Y * Xindex + X] = 0;//initialize
              if (static_cast<int>(avolInfoEx->isDomain[Y * Xindex + X]) == 1) {
                if (X == 0 || X == Xindex - 1 || Y == 0 || Y == Yindex - 1) {//boundary of the domain
                  avolInfoEx->isBoundary[Y * Xindex + X] = 1;
                } else if (static_cast<int>(avolInfoEx->isDomain[Y * Xindex + (X + 2)]) == 0
                  || static_cast<int>(avolInfoEx->isDomain[Y * Xindex + (X - 2)]) == 0
                  || static_cast<int>(avolInfoEx->isDomain[(Y + 2) * Xindex + X]) == 0
                  || static_cast<int>(avolInfoEx->isDomain[(Y - 2) * Xindex + X]) == 0) {
                    avolInfoEx->isBoundary[Y * Xindex + X] = 1;
                }
              }
            }
          }
          break;
        case 3://3D
          for (Z = 0; Z < Zindex; Z += 2) {
            for (Y = 0; Y < Yindex; Y += 2) {
              for (X = 0; X < Xindex; X += 2) {
                avolInfoEx->isBoundary[Z * Yindex * Xindex + Y * Xindex + X] = 0;//initialize
                if (static_cast<int>(avolInfoEx->isDomain[Z * Yindex * Xindex + Y * Xindex + X]) == 1) {
                  if (X == 0 || X == Xindex - 1 || Y == 0 || Y == Yindex - 1 || Z == 0 || Z == Zindex - 1) {//boundary of th domain
                    avolInfoEx->isBoundary[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                  } else if (static_cast<int>(avolInfoEx->isDomain[Z * Yindex * Xindex + Y * Xindex + (X + 2)]) == 0
                    || static_cast<int>(avolInfoEx->isDomain[Z * Yindex * Xindex + Y * Xindex + (X - 2)]) == 0
                    || static_cast<int>(avolInfoEx->isDomain[Z * Yindex * Xindex + (Y + 2) * Xindex + X]) == 0
                    || static_cast<int>(avolInfoEx->isDomain[Z * Yindex * Xindex + (Y - 2) * Xindex + X]) == 0
                    || static_cast<int>(avolInfoEx->isDomain[(Z + 2) * Yindex * Xindex + Y * Xindex + X]) == 0
                    || static_cast<int>(avolInfoEx->isDomain[(Z - 2) * Yindex * Xindex + Y * Xindex + X]) == 0) {
                      avolInfoEx->isBoundary[Z * Yindex * Xindex + Y * Xindex + X] = 1;
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
#pragma endregion 
  //membrane compartment
  //adjacentとかdomainとか使って隣接関係を割り出す。
  for (i = 0; i < numOfCompartments; i++) {
    if (SpatialCompartmentPlugin *cPlugin = static_cast<SpatialCompartmentPlugin*>(loc->get(i)->getPlugin("spatial"))) {
      DomainType *dType = geometry->getDomainType(cPlugin->getCompartmentMapping()->getDomainType());
      for (j = 0; j < geometry->getNumDomains(); j++) {//Domain
        if (geometry->getDomain(j)->getDomainType() == cPlugin->getCompartmentMapping()->getDomainType()) {
          analyticVolInfo *avolInfo = 0;
          if (!(avolInfo = searchAvolInfoByCompartment(avolInfoList, loc->get(i)->getId().c_str()))) {
            avolInfo = new analyticVolInfo;
            InitializeAVolInfo(avolInfo);
            avolInfo->compartmentId = loc->get(i)->getId().c_str();
            avolInfo->domainTypeId = cPlugin->getCompartmentMapping()->getDomainType().c_str();
            avolInfo->domainId = geometry->getDomain(j)->getSpatialId().c_str();
            avolInfo->bt = new boundaryType[numOfVolIndexes];
            if (dType->getSpatialDimensions() == volDimension) {
              avolInfo->isVol = true;
            } else if (dType->getSpatialDimensions() == memDimension) {
              avolInfo->isVol = false;
            }
            avolInfo->adjacent0 = 0;
            avolInfo->adjacent1 = 0;
            avolInfo->isDomain = 0;
            avolInfo->isBoundary = 0;
            avolInfo->rpInfo = 0;
            avolInfoList.push_back(avolInfo);
          } else {
            avolInfo->domainId = geometry->getDomain(j)->getSpatialId().c_str();
          }
          for (k = 0; k < geometry->getNumAdjacentDomains(); k++) {//AdjacentDomain
            AdjacentDomains *adDomain = geometry->getAdjacentDomains(k);
            if (adDomain->getDomain1() == geometry->getDomain(j)->getSpatialId() || adDomain->getDomain2() == geometry->getDomain(j)->getSpatialId()) {
              if (avolInfo->adjacent0 == 0) {
                avolInfo->adjacent0 = adDomain;
              } else if (avolInfo->adjacent1 == 0) {
                avolInfo->adjacent1 = adDomain;
              }
            }
          }
        }
      }
    }
  }

  //膜の位置
  for (i = 0; i < avolInfoList.size(); i++) {
    analyticVolInfo *avolInfo = avolInfoList[i];
    if (avolInfo->isVol == false) {//avol is membrane
      analyticVolInfo *ad0 = NULL, *ad1 = NULL;
      //adjacent0
      if (strcmp(avolInfo->adjacent0->getDomain1().c_str(), avolInfo->domainId) == 0) {
        ad0 = searchAvolInfoByDomain(avolInfoList, avolInfo->adjacent0->getDomain2().c_str());
      } else if (strcmp(avolInfo->adjacent0->getDomain2().c_str(), avolInfo->domainId) == 0) {
        ad0 = searchAvolInfoByDomain(avolInfoList, avolInfo->adjacent0->getDomain1().c_str());
      }
      //adjacent1
      if (strcmp(avolInfo->adjacent0->getDomain1().c_str(), avolInfo->domainId) == 0) {
        ad1 = searchAvolInfoByDomain(avolInfoList, avolInfo->adjacent1->getDomain2().c_str());
      } else if (strcmp(avolInfo->adjacent0->getDomain2().c_str(), avolInfo->domainId) == 0) {
        ad1 = searchAvolInfoByDomain(avolInfoList, avolInfo->adjacent1->getDomain1().c_str());
      }

      avolInfo->isDomain = new double[numOfVolIndexes];
      INIT_DOUBLE(avolInfo->isDomain, numOfVolIndexes);
      avolInfo->isBoundary = new int[numOfVolIndexes];
      INIT_INT(avolInfo->isBoundary, numOfVolIndexes);
      if (ad0 == NULL && ad1 == NULL) 
        continue;
      if (ad0 == NULL && ad1 != NULL) 
      {
        memcpy(avolInfo->isDomain, ad1->isDomain, sizeof(double)*numOfVolIndexes);
        continue;
      }
      if (ad0 != NULL && ad1 == NULL) 
      {
        memcpy(avolInfo->isDomain, ad0->isDomain, sizeof(double)*numOfVolIndexes);
        continue;
      }


      // 			if (ad0->isEdge == 0) {
      // 				ad0->isEdge = new bool[numOfVolIndexes];
      // 				for (j = 0; j < static_cast<unsigned int>(numOfVolIndexes); j++) {
      // 					ad0->isEdge[j] = false;
      // 				}
      // 			}
      // 			if (ad1->isEdge == 0) {
      // 				ad1->isEdge = new bool[numOfVolIndexes];
      // 				for (j = 0; j < static_cast<unsigned int>(numOfVolIndexes); j++) {
      // 					ad1->isEdge[j] = false;
      // 				}
      // 			}
      switch (dimension) {
      case 1:
        for (X = 0; X < Xindex; X++) {
          if (X % 2 != 0) {
            Xplus = Y * Xindex + (X + 1);
            Xminus = Y * Xindex + (X - 1);
            if (X != 0 && X != (Xindex - 1) && Y != 0 && Y != (Yindex - 1) && Z != 0 && Z != (Zindex - 1)) {
              if ((ad0->isDomain[Xplus] == 1 && ad1->isDomain[Xminus] == 1) ||
                (ad0->isDomain[Xminus] == 1 && ad1->isDomain[Xplus] == 1)) {
                  avolInfo->isDomain[X] = 1;
                  // 								ad0->isDomain[X] = 1;
                  // 								ad0->isEdge[X] = true;
                  // 								ad1->isDomain[X] = 1;
                  // 								ad1->isEdge[X] = true;
              }
            }
          }
        }
        break;
      case 2:
        for (Y = 0; Y < Yindex; Y++) {
          for (X = 0; X < Xindex; X++) {
            if ((Y * Xindex + X) % 2 != 0) {
              Xplus = Y * Xindex + (X + 1);
              Xminus = Y * Xindex + (X - 1);
              Yplus = (Y + 1) * Xindex + X;
              Yminus = (Y - 1) * Xindex + X;
              if (X != 0 && X != (Xindex - 1) && Y != 0 && Y != (Yindex - 1)) {
                if ((ad0->isDomain[Xplus] == 1 && ad1->isDomain[Xminus] == 1) ||
                  (ad0->isDomain[Xminus] == 1 && ad1->isDomain[Xplus] == 1) ||
                  (ad0->isDomain[Yplus] == 1 && ad1->isDomain[Yminus] == 1) ||
                  (ad0->isDomain[Yminus] == 1 && ad1->isDomain[Yplus] == 1)) {
                    avolInfo->isDomain[Y * Xindex + X] = 1;
                    // 									ad0->isDomain[Y * Xindex + X] = 1;
                    // 									ad0->isEdge[Y * Xindex + X] = true;
                    // 									ad1->isDomain[Y * Xindex + X] = 1;
                    // 									ad1->isEdge[Y * Xindex + X] = true;
                }
              } else if (X == 0 || X == Xindex - 1) {
                if ((ad0->isDomain[Yplus] == 1 && ad1->isDomain[Yminus] == 1) ||
                  (ad0->isDomain[Yminus] == 1 && ad1->isDomain[Yplus] == 1)) {
                    avolInfo->isDomain[Y * Xindex + X] = 1;
                    // 									ad0->isDomain[Y * Xindex + X] = 1;
                    // 									ad0->isEdge[Y * Xindex + X] = true;
                    // 									ad1->isDomain[Y * Xindex + X] = 1;
                    // 									ad1->isEdge[Y * Xindex + X] = true;
                }
              } else if (Y == 0 || Y == Yindex - 1) {
                if ((ad0->isDomain[Xplus] == 1 && ad1->isDomain[Xminus] == 1) ||
                  (ad0->isDomain[Xminus] == 1 && ad1->isDomain[Xplus] == 1)) {
                    avolInfo->isDomain[Y * Xindex + X] = 1;
                    // 									ad0->isDomain[Y * Xindex + X] = 1;
                    // 									ad0->isEdge[Y * Xindex + X] = true;
                    // 									ad1->isDomain[Y * Xindex + X] = 1;
                    // 									ad1->isEdge[Y * Xindex + X] = true;
                }
              }
            }
          }
        }
        //間の点を補完
        for (Y = 0; Y < Yindex; Y++) {
          for (X = 0; X < Xindex; X++) {
            //if (!(X % 2 == 0 && Y % 2 == 0)) {
            if ((X % 2 != 0 && Y % 2 != 0)) {
              Xplus = Y * Xindex + (X + 1);
              Xminus = Y * Xindex + (X - 1);
              Yplus = (Y + 1) * Xindex + X;
              Yminus = (Y - 1) * Xindex + X;
              if (X != 0 && X != Xindex - 1 && Y != 0 && Y != Yindex - 1){
                if ((avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Xminus] == 1) ||
                  (avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Yminus] == 1)) {
                    avolInfo->isDomain[Y * Xindex + X] = 2;
                }
                if ((avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Yplus] == 1) ||
                  (avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Yminus] == 1) ||
                  (avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Xminus] == 1) ||
                  (avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Xminus] == 1)) {
                    avolInfo->isDomain[Y * Xindex + X] = 2;
                }
              } else if (X == 0 || X == Xindex - 1) {
                if (avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Yminus] == 1) {
                  avolInfo->isDomain[Y * Xindex + X] = 2;
                }
              } else if (Y == 0 || Y == Yindex - 1) {
                if (avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Xminus] == 1) {
                  avolInfo->isDomain[Y * Xindex + X] = 2;
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
                Xplus = Z * Yindex * Xindex + Y * Xindex + (X + 1);
                Xminus = Z * Yindex * Xindex + Y * Xindex + (X - 1);
                Yplus = Z * Yindex * Xindex + (Y + 1) * Xindex + X;
                Yminus = Z * Yindex * Xindex + (Y - 1) * Xindex + X;
                Zplus = (Z + 1) * Yindex * Xindex + Y * Xindex + X;
                Zminus = (Z - 1) * Yindex * Xindex + Y * Xindex + X;
                if (X != 0 && X != (Xindex - 1) && Y != 0 && Y != (Yindex - 1) && Z != 0 && Z != (Zindex - 1)) {
                  if ((ad0->isDomain[Xplus] == 1 && ad1->isDomain[Xminus] == 1) ||
                    (ad0->isDomain[Xminus] == 1 && ad1->isDomain[Xplus] == 1) ||
                    (ad0->isDomain[Yplus] == 1 && ad1->isDomain[Yminus] == 1) ||
                    (ad0->isDomain[Yminus] == 1 && ad1->isDomain[Yplus] == 1) ||
                    (ad0->isDomain[Zplus] == 1 && ad1->isDomain[Zminus] == 1) ||
                    (ad0->isDomain[Zminus] == 1 && ad1->isDomain[Zplus] == 1)) {
                      avolInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      // 										ad0->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      // 										ad0->isEdge[Z * Yindex * Xindex + Y * Xindex + X] = true;
                      // 										ad1->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      // 										ad1->isEdge[Z * Yindex * Xindex + Y * Xindex + X] = true;
                  }
                } else if (X == 0 || X == Xindex - 1) {
                  if ((ad0->isDomain[Yplus] == 1 && ad1->isDomain[Yminus] == 1) ||
                    (ad0->isDomain[Yminus] == 1 && ad1->isDomain[Yplus] == 1) ||
                    (ad0->isDomain[Zplus] == 1 && ad1->isDomain[Zminus] == 1) ||
                    (ad0->isDomain[Zminus] == 1 && ad1->isDomain[Zplus] == 1)) {
                      avolInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      // 										ad0->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      // 										ad0->isEdge[Z * Yindex * Xindex + Y * Xindex + X] = true;
                      // 										ad1->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      // 										ad1->isEdge[Z * Yindex * Xindex + Y * Xindex + X] = true;
                  }
                } else if (Y == 0 || Y == Yindex - 1) {
                  if ((ad0->isDomain[Xplus] == 1 && ad1->isDomain[Xminus] == 1) ||
                    (ad0->isDomain[Xminus] == 1 && ad1->isDomain[Xplus] == 1) ||
                    (ad0->isDomain[Zplus] == 1 && ad1->isDomain[Zminus] == 1) ||
                    (ad0->isDomain[Zminus] == 1 && ad1->isDomain[Zplus] == 1)) {
                      avolInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      // 										ad0->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      // 										ad0->isEdge[Z * Yindex * Xindex + Y * Xindex + X] = true;
                      // 										ad1->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      // 										ad1->isEdge[Z * Yindex * Xindex + Y * Xindex + X] = true;
                  }
                } else if (Z == 0 || Z == Zindex - 1) {
                  if ((ad0->isDomain[Xplus] == 1 && ad1->isDomain[Xminus] == 1) ||
                    (ad0->isDomain[Xminus] == 1 && ad1->isDomain[Xplus] == 1) ||
                    (ad0->isDomain[Yplus] == 1 && ad1->isDomain[Yminus] == 1) ||
                    (ad0->isDomain[Yminus] == 1 && ad1->isDomain[Yplus] == 1)) {
                      avolInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      // 										ad0->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      // 										ad0->isEdge[Z * Yindex * Xindex + Y * Xindex + X] = true;
                      // 										ad1->isDomain[Z * Yindex * Xindex + Y * Xindex + X] = 1;
                      // 										ad1->isEdge[Z * Yindex * Xindex + Y * Xindex + X] = true;
                  }
                }
              }
            }
          }
        }
        //間の点を補完
        for (Z = 0; Z < Zindex; Z++) {
          for (Y = 0; Y < Yindex; Y++) {
            for (X = 0; X < Xindex; X++) {
              if ((X % 2 != 0 && Y % 2 != 0) || (Y % 2 != 0 && Z % 2 != 0) || (Z % 2 != 0 && X % 2 != 0)) {
                index = Z * Yindex * Xindex + Y * Xindex + X;
                Xplus = Z * Yindex * Xindex + Y * Xindex + (X + 1);
                Xminus = Z * Yindex * Xindex + Y * Xindex + (X - 1);
                Yplus = Z * Yindex * Xindex + (Y + 1) * Xindex + X;
                Yminus = Z * Yindex * Xindex + (Y - 1) * Xindex + X;
                Zplus = (Z + 1) * Yindex * Xindex + Y * Xindex + X;
                Zminus = (Z - 1) * Yindex * Xindex + Y * Xindex + X;
                if (X != 0 && X != Xindex - 1 && Y != 0 && Y != Yindex - 1 && Z != 0 && Z != Zindex - 1){
                  if ((avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Xminus] == 1) ||
                    (avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Yminus] == 1) ||
                    (avolInfo->isDomain[Zplus] == 1 && avolInfo->isDomain[Zminus] == 1)) {
                      avolInfo->isDomain[index] = 2;
                  }
                  if ((avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Yplus] == 1) ||
                    (avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Yminus] == 1) ||
                    (avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Zplus] == 1) ||
                    (avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Zminus] == 1) ||
                    (avolInfo->isDomain[Xminus] == 1 && avolInfo->isDomain[Yplus] == 1) ||
                    (avolInfo->isDomain[Xminus] == 1 && avolInfo->isDomain[Yminus] == 1) ||
                    (avolInfo->isDomain[Xminus] == 1 && avolInfo->isDomain[Zplus] == 1) ||
                    (avolInfo->isDomain[Xminus] == 1 && avolInfo->isDomain[Zminus] == 1) ||
                    (avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Zplus] == 1) ||
                    (avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Zminus] == 1) ||
                    (avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Zplus] == 1) ||
                    (avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Zminus] == 1)) {
                      avolInfo->isDomain[index] = 2;
                  }
                } else if (X == 0 || X == Xindex - 1) {
                  if ((avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Yminus] == 1) ||
                    (avolInfo->isDomain[Zplus] == 1 && avolInfo->isDomain[Zminus] == 1) ||
                    (avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Zplus] == 1) ||
                    (avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Zminus] == 1) ||
                    (avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Zplus] == 1) ||
                    (avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Zminus] == 1)) {
                      avolInfo->isDomain[index] = 2;
                  }
                } else if (Y == 0 || Y == Yindex - 1) {
                  if ((avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Xminus] == 1) ||
                    (avolInfo->isDomain[Zplus] == 1 && avolInfo->isDomain[Zminus] == 1) ||
                    (avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Zplus] == 1) ||
                    (avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Zminus] == 1) ||
                    (avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Zplus] == 1) ||
                    (avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Zminus] == 1)) {
                      avolInfo->isDomain[index] = 2;
                  }
                } else if (Z == 0 || Z == Zindex - 1) {
                  if ((avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Xminus] == 1) ||
                    (avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Yminus] == 1)) {
                      avolInfo->isDomain[index] = 2;
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
  /*
  //output geometries
  int *geo_edge = new int[numOfVolIndexes];
  for (i = 0; i < avolInfoList.size(); i++) {
  analyticVolInfo *avolInfo = avolInfoList[i];
  if (avolInfo->isVol == false) {//avol is membrane
  for (Z = 0; Z < Zindex; Z++) {
  for (Y = 0; Y < Yindex; Y++) {
  for (X = 0; X < Xindex; X++) {
  if (avolInfo->isDomain[Z * Yindex * Xindex + Y * Xindex + X] > 0) {
  geo_edge[Z * Yindex * Xindex + Y * Xindex + X] = 1;
  }
  }
  }
  }
  }
  }

  ofstream ofs;
  string geo_filename = "./result/txt/geo_sum.txt";
  ofs.open(geo_filename.c_str());
  switch(dimension) {
  case 1:
  for (X = 0; X < Xindex; X++) {
  if (geo_edge[X] == 1) {
  ofs << xInfo->value[X] << " " << 0 << " " << 1 << " " << 150 << " " << 150 << " " << 150 << " " << 255 << endl;
  } else if (geo_edge[index] == 0) {
  ofs << xInfo->value[X] << " " << 0 << " " << 0 << " " << 150 << " " << 150 << " " << 150 << " " << 0 << endl;
  }
  }
  break;
  case 2:
  for (Y = 0; Y < Yindex; Y++) {
  for (X = 0; X < Xindex; X++) {
  index = Y * Xindex + X;
  if (geo_edge[index] == 1) {
  ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << 1 << " " << 150 << " " << 150 << " " << 150 << " " << 255 << endl;
  //ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << 0 << " " << 150 << " " << 150 << " " << 150 << " " << 0 << endl;
  } else if (geo_edge[index] == 0) {
  ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << 0 << " " << 150 << " " << 150 << " " << 150 << " " << 0 << endl;
  }
  }
  ofs << endl;
  }
  break;
  case 3:
  for (Z = 0; Z < Zindex; Z++) {
  for (Y = 0; Y < Yindex; Y++) {
  for (X = 0; X < Xindex; X++) {
  index = Z * Yindex * Xindex + Y * Xindex + X;
  if (geo_edge[index] == 1) {
  ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << zInfo->value[index] << " " << 1 << " " << 150 << " " << 150 << " " << 150 << " " << 255 << endl;
  } else if (geo_edge[index] == 0) {
  ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << zInfo->value[index] << " " << 0 << " " << 150 << " " << 150 << " " << 150 << " " << 0 << endl;
  }
  }
  ofs << endl;
  }
  ofs << endl;
  }
  break;
  default:
  break;
  }
  ofs.close();
  delete[] geo_edge;

  //output membrane
  for (i = 0; i < avolInfoList.size(); i++) {
  analyticVolInfo *avolInfo = avolInfoList[i];
  if (avolInfo->isVol == false) {//avol is membrane
  string tmp_dtId = avolInfo->domainTypeId;
  string geo_filename = "";
  geo_filename = "./result/txt/" + tmp_dtId + ".txt";
  memList.push_back(avolInfo->domainTypeId);
  ofs.open(geo_filename.c_str());
  switch(dimension) {
  case 1:
  for (X = 0; X < Xindex; X++) {
  if (static_cast<int>(avolInfo->isDomain[X]) > 0) {
  ofs << xInfo->value[X] << " " << 0 << " " << 1 << " " << 150 << " " << 150 << " " << 150 << " " << 255 << endl;
  } else {
  ofs << xInfo->value[X] << " " << 0 << " " << 0 << " " << 150 << " " << 150 << " " << 150 << " " << 0 << endl;
  }
  }
  break;
  case 2:
  for (Y = 0; Y < Yindex; Y++) {
  for (X = 0; X < Xindex; X++) {
  index = Y * Xindex + X;
  if (static_cast<int>(avolInfo->isDomain[index]) > 0) {
  ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << 1 << " " << 150 << " " << 150 << " " << 150 << " " << 255 << endl;
  } else {
  ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << 0 << " " << 150 << " " << 150 << " " << 150 << " " << 0 << endl;
  }
  }
  ofs << endl;
  }
  break;
  case 3:
  for (Z = 0; Z < Zindex; Z++) {
  for (Y = 0; Y < Yindex; Y++) {
  for (X = 0; X < Xindex; X++) {
  index = Z * Yindex * Xindex + Y * Xindex + X;
  if (static_cast<int>(avolInfo->isDomain[index]) > 0) {
  ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << zInfo->value[index] << " " << 1 << " " << 150 << " " << 150 << " " << 150 << " " << 255 << endl;
  } else {
  ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << zInfo->value[index] << " " << 0 << " " << 150 << " " << 150 << " " << 150 << " " << 0 << endl;
  }
  }
  ofs << endl;
  }
  ofs << endl;
  }
  break;
  default:
  break;
  }
  ofs.close();
  }
  }
  */
  // 	vector<memArrayIndex*> memArrayIndexList = vector<memArrayIndex*>();
  // 	//info->value = new double[numOfVolIndexes];
  // 	vector<memArrayIndex*> mailZList = new vector<memArrayIndex*>[Zindex];
  // 	//membraneの座標の配列


  //set species' initial condition
  //boundary type(Xp, Xm, Yp, Ym, Zx, Zm)
  //parse dependence among species, compartments, parameters
  list<variableInfo*> notOrderedInfo;
  vector<variableInfo*> orderedInfo;
  for (i = 0; i < varInfoList.size(); i++) {
    variableInfo *info = varInfoList[i];
    ast = 0;
    if (model->getInitialAssignment(info->id) != 0) {//initial assignment
      if (info->value == 0) {//value is not set yet
        info->value = new double[numOfVolIndexes];
        INIT_DOUBLE(info->value, numOfVolIndexes);
        if (info->sp != 0 && (!info->sp->isSetConstant() || !info->sp->getConstant())) {
          info->delta = new double[4 * numOfVolIndexes];
          INIT_DOUBLE(info->delta, 4*numOfVolIndexes);          
        }
      }
      ast = const_cast<ASTNode*>((model->getInitialAssignment(info->id))->getMath());
    } else if (model->getRule(info->id) != 0 && model->getRule(info->id)->isAssignment()) {//assignment rule
      info->hasAssignmentRule = true;
      if (info->value == 0) {//value is not set yet
        info->value = new double[numOfVolIndexes];
        INIT_DOUBLE(info->value, numOfVolIndexes);
        if (info->sp != 0 && (!info->sp->isSetConstant() || !info->sp->getConstant())) {
          //the species is variable
          info->delta = new double[4 * numOfVolIndexes];
          INIT_DOUBLE(info->delta, 4*numOfVolIndexes);
        }
      }
      ast = const_cast<ASTNode*>(((AssignmentRule*)model->getRule(info->id))->getMath());
    }
    if (ast != 0) {
#ifdef PRINT_DEBUG
      cout << info->id << ": " << SBML_formulaToString(ast) << endl;
#endif
      rearrangeAST(ast);
      numOfASTNodes = 0;
      countAST(ast, numOfASTNodes);
      info->rpInfo = new reversePolishInfo();
      info->rpInfo->varList = new double*[numOfASTNodes];
      info->rpInfo->deltaList = 0;
      info->rpInfo->constList = new double*[numOfASTNodes];
      info->rpInfo->opfuncList = new int[numOfASTNodes];
      info->rpInfo->listNum = numOfASTNodes;
      info->isResolved = false;
      parseAST(ast, info->rpInfo->varList, info->rpInfo->deltaList, info->rpInfo->constList, info->rpInfo->opfuncList, varInfoList, info->rpInfo->listNum);
      parseDependence(ast, info->dependence, varInfoList);
      notOrderedInfo.push_back(info);
    }
  }
  //依存関係を考慮しつつ計算(未テスト)
  while (!notOrderedInfo.empty()) {
    list<variableInfo*>::reverse_iterator it = notOrderedInfo.rbegin();
    while (it != notOrderedInfo.rend()) {
      variableInfo* current = (*it);
      if (isResolvedAll(current->dependence)) {
        reversePolishInitial(
          current->rpInfo->varList, 
          current->rpInfo->constList, 
          current->rpInfo->opfuncList, 
          current->value, 
          current->rpInfo->listNum, 
          Xindex, Yindex, Zindex, 
          current->inVol);
        current->isResolved = true;
        //orderedInfo.push_back(*it);
        //notOrderedInfo.erase(--(it.base()));     
        notOrderedInfo.remove(current);

#if !(defined _MSC_VER)
        // removing already updates the iterator
        //it++;
#endif
      }      
    }
  }

  for (i = 0; i < numOfSpecies; i++) {
    Species *s = los->get(i);
    variableInfo *sInfo = searchInfoById(varInfoList, s->getId().c_str());
    if (sInfo != 0) {
      if ((sInfo->avi = searchAvolInfoByCompartment(avolInfoList, s->getCompartment().c_str()))) {
        double *isD = sInfo->avi->isDomain;
        int *isB = sInfo->avi->isBoundary;
        switch (dimension) {
        case 1:
          if (sInfo->avi->isVol == true) {
            for (X = 0; X < Xindex; X += 2) {
              if (static_cast<int>(sInfo->avi->isDomain[X]) == 0) {
                sInfo->value[X] = 0.0;
              }
              if (isD[X] == 1 && isB[X] == 1) {
                if (X == Xindex - 1) {//the edge of simulation area
                  sInfo->avi->bt[X] = Xp;
                } else if (X == 0) {
                  sInfo->avi->bt[index] = Xm;
                } else {//not the edge of simulation area
                  if (isD[X - 2] == 1 && isD[X + 2] == 0 && isB[X - 2] == 0) {//Xp
                    sInfo->avi->bt[X] = Xp;
                  } else if (isD[X + 2] == 1 && isD[X - 2] == 0 && isB[X + 2] == 0) {//Xm
                    sInfo->avi->bt[X] = Xm;
                  }
                }
              }
            }
          } else {//membrane
            for (X = 0; X < Xindex; X++) {
              Xplus = Y * Xindex + (X + 1);
              Xminus = Y * Xindex + (X - 1);
              if (static_cast<int>(sInfo->avi->isDomain[X]) == 0) {
                sInfo->value[X] = 0.0;
              }
              if (X % 2 != 0) {
                if (sInfo->avi->isDomain[X] == 2 && sInfo->avi->isDomain[Xplus] == 1 && sInfo->avi->isDomain[Xminus] == 1) {
                  sInfo->value[X] = (sInfo->value[Xplus] + sInfo->value[Xminus]) / 2.0;
                }
              }
            }
          }
          break;
        case 2:
          if (sInfo->avi->isVol == true) {//volume
            for (Y = 0; Y < Yindex; Y += 2) {
              for (X = 0; X < Xindex; X += 2) {
                index = Y * Xindex + X;
                Xplus = Y * Xindex + (X + 2);
                Xminus = Y * Xindex + (X - 2);
                Yplus = (Y + 2) * Xindex + X;
                Yminus = (Y - 2) * Xindex + X;
                if (static_cast<int>(sInfo->avi->isDomain[index]) == 0) {
                  sInfo->value[index] = 0.0;
                }
                if (isD[index] == 1 && isB[index] == 1) {
                  if (X == Xindex - 1) {//the edge of simulation area
                    sInfo->avi->bt[index] = Xp;
                  } else if (X == 0) {
                    sInfo->avi->bt[index] = Xm;
                  } else if (Y == Yindex - 1) {
                    sInfo->avi->bt[index] = Yp;
                  } else if (Y == 0) {
                    sInfo->avi->bt[index] = Ym;
                  } else {//not the edge of simulation area
                    if ((isD[Xminus] == 1 && isD[Xplus] == 0) &&
                      (isB[Xminus] == 0 ||
                      ((isD[Yplus] == 0 || isB[Yplus] == 1 ) && (isD[Yminus] == 0 || isB[Yminus] == 1)))) {//Xp
                        sInfo->avi->bt[index] = Xp;
                        //sInfo->value[index] = 2.0;
                    } else if ((isD[Xplus] == 1 && isD[Xminus] == 0) &&
                      (isB[Xplus] == 0 ||
                      ((isD[Yplus] == 0 || isB[Yplus] == 1) && (isD[Yminus] == 0 || isB[Yminus] == 1)))) {//Xm
                        sInfo->avi->bt[index] = Xm;
                        //sInfo->value[index] = 3.0;
                    } else if ((isD[Yminus] == 1 && isD[Yplus] == 0) &&
                      (isB[Yminus] == 0 ||
                      ((isD[Xplus] == 0 || isB[Xplus] == 1) && (isD[Xminus] == 0 || isB[Xminus] == 1)))) {//Yp
                        sInfo->avi->bt[index] = Yp;
                        //sInfo->value[index] = 4.0;
                    } else if ((isD[Yplus] == 1 && isD[Yminus] == 0) &&
                      (isB[Yplus] == 0 ||
                      ((isD[Xplus] == 0 || isB[Xplus] == 1) && (isD[Xminus] == 0 || isB[Xminus] == 1)))) {//Ym
                        sInfo->avi->bt[index] = Ym;
                        //sInfo->value[index] = 5.0;
                    } else {
                      //sInfo->value[index] = 1.0;
                    }
                  }
                }
              }
            }
          } else {//membrane(角の値を補完)
            for (Y = 0; Y < Yindex; Y++) {
              for (X = 0; X < Xindex; X++) {
                index = Y * Xindex + X;
                Xplus = Y * Xindex + (X + 1);
                Xminus = Y * Xindex + (X - 1);
                Yplus = (Y + 1) * Xindex + X;
                Yminus = (Y - 1) * Xindex + X;
                if (static_cast<int>(sInfo->avi->isDomain[index]) == 0) {
                  sInfo->value[index] = 0.0;
                }
                if (X % 2 != 0 && Y % 2 != 0) {
                  if (sInfo->avi->isDomain[index] == 2) {
                    if (sInfo->avi->isDomain[Xplus] == 1 && sInfo->avi->isDomain[Xminus] == 1) {
                      sInfo->value[index] = (sInfo->value[Xplus] + sInfo->value[Xminus]) / 2.0;
                    } else if (sInfo->avi->isDomain[Yplus] == 1 && sInfo->avi->isDomain[Yminus] == 1) {
                      sInfo->value[index] = (sInfo->value[Yplus] + sInfo->value[Yminus]) / 2.0;
                    } else if (sInfo->avi->isDomain[Xplus] == 1 && sInfo->avi->isDomain[Yplus] == 1) {
                      sInfo->value[index] = (sInfo->value[Xplus] + sInfo->value[Yplus]) / 2.0;
                    } else if (sInfo->avi->isDomain[Xplus] == 1 && sInfo->avi->isDomain[Yminus] == 1) {
                      sInfo->value[index] = (sInfo->value[Xplus] + sInfo->value[Yminus]) / 2.0;
                    } else if (sInfo->avi->isDomain[Xminus] == 1 && sInfo->avi->isDomain[Yplus] == 1) {
                      sInfo->value[index] = (sInfo->value[Xminus] + sInfo->value[Yplus]) / 2.0;
                    } else if (sInfo->avi->isDomain[Xminus] == 1 && sInfo->avi->isDomain[Yminus] == 1) {
                      sInfo->value[index] = (sInfo->value[Xminus] + sInfo->value[Yminus]) / 2.0;
                    }
                  }
                }
              }
            }
          }
          break;
        case 3:
          if (sInfo->avi->isVol == true) {//volume
            for (Z = 0; Z < Zindex; Z += 2) {
              for (Y = 0; Y < Yindex; Y += 2) {
                for (X = 0; X < Xindex; X += 2) {
                  index = Z * Yindex * Xindex + Y * Xindex + X;
                  Xplus = Z * Yindex * Xindex + Y * Xindex + (X + 2);
                  Xminus = Z * Yindex * Xindex + Y * Xindex + (X - 2);
                  Yplus = Z * Yindex * Xindex + (Y + 2) * Xindex + X;
                  Yminus = Z * Yindex * Xindex + (Y - 2) * Xindex + X;
                  Zplus = (Z + 2) * Yindex * Xindex + Y * Xindex + X;
                  Zminus = (Z - 2) * Yindex * Xindex + Y * Xindex + X;
                  if (static_cast<int>(sInfo->avi->isDomain[index]) == 0) {
                    sInfo->value[index] = 0.0;
                  }
                  if (isD[index] == 1 && isB[index] == 1) {
                    if (X == Xindex - 1) {//the edge of simulation area
                      sInfo->avi->bt[index] = Xp;
                    } else if (X == 0) {
                      sInfo->avi->bt[index] = Xm;
                    } else if (Y == Yindex - 1) {
                      sInfo->avi->bt[index] = Yp;
                    } else if (Y == 0) {
                      sInfo->avi->bt[index] = Ym;
                    } else if (Z == Zindex - 1) {
                      sInfo->avi->bt[index] = Zp;
                    } else if (Z == 0) {
                      sInfo->avi->bt[index] = Zm;
                    } else {//not the edge of simulation area
                      if ((isD[Xminus] == 1 && isD[Xplus] == 0) &&
                        (isB[Xminus] == 0 ||
                        ((isD[Yplus] == 0 || isB[Yplus] == 1) && (isD[Yminus] == 0 || isB[Yminus] == 1) &&
                        (isD[Zplus] == 0 || isB[Zplus] == 1)  && (isD[Zminus] == 0 || isB[Zminus] == 1)))) {//Xp
                          sInfo->avi->bt[index] = Xp;
                          //sInfo->value[index] = 2.0;
                      } else if ((isD[Xplus] == 1 && isD[Xminus] == 0) &&
                        (isB[Xplus] == 0 ||
                        ((isD[Yplus] == 0 || isB[Yplus] == 1) && (isD[Yminus] == 0 || isB[Yminus] == 1) &&
                        (isD[Zplus] == 0 || isB[Zplus] == 1) && (isD[Zminus] == 0 || isB[Zminus] == 1)))) {//Xm
                          sInfo->avi->bt[index] = Xm;
                          //sInfo->value[index] = 3.0;
                      } else if ((isD[Yminus] == 1 && isD[Yplus] == 0) &&
                        (isB[Yminus] == 0 ||
                        ((isD[Xplus] == 0 || isB[Xplus] == 1) && (isD[Xminus] == 0 || isB[Xminus] == 1) &&
                        (isD[Zplus] == 0 || isB[Zplus] == 1) && (isD[Zminus] == 0 || isB[Zminus] == 1)))) {//Yp
                          sInfo->avi->bt[index] = Yp;
                          //sInfo->value[index] = 4.0;
                      } else if ((isD[Yplus] == 1 && isD[Yminus] == 0) &&
                        (isB[Yplus] == 0 ||
                        ((isD[Xplus] == 0 || isB[Xplus] == 1) && (isD[Xminus] == 0 || isB[Xminus] == 1) &&
                        (isD[Zplus] == 0 || isB[Zplus] == 1) && (isD[Zminus] == 0 || isB[Zminus] == 1)))) {//Ym
                          sInfo->avi->bt[index] = Ym;
                          //sInfo->value[index] = 5.0;
                      } else if ((isD[Zminus] == 1 && isD[Zplus] == 0) &&
                        (isB[Zminus] == 0 ||
                        ((isD[Xplus] == 0 || isB[Xplus] == 1) && (isD[Xminus] == 0 || isB[Xminus] == 1) &&
                        (isD[Yplus] == 0 || isB[Yplus] == 1) && (isD[Yminus] == 0 || isB[Yminus] == 1)))) {//Zp
                          sInfo->avi->bt[index] = Zp;
                          //sInfo->value[index] = 6.0;
                      } else if ((isD[Zplus] == 1 && isD[Zminus] == 0) &&
                        (isB[Zplus] == 0 ||
                        ((isD[Xplus] == 0 || isB[Xplus] == 1) && (isD[Xminus] == 0 || isB[Xminus] == 1) &&
                        (isD[Yplus] == 0 || isB[Yplus] == 1) && (isD[Yminus] == 0 || isB[Yminus] == 1)))) {//Zm
                          sInfo->avi->bt[index] = Zm;
                          //sInfo->value[index] = 7.0;
                      }
                    }
                  }
                }
              }
            }
          } else {//membrane
            /*
            if ((X % 2 != 0 && Y % 2 != 0) || (Y % 2 != 0 && Z % 2 != 0) || (Z % 2 != 0 && X % 2 != 0)) {
            if (X != 0 && X != Xindex - 1 && Y != 0 && Y != Yindex - 1 && Z != 0 && Z != Zindex - 1){
            if (avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Xminus] == 1 ||
            avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Yminus] == 1 ||
            avolInfo->isDomain[Zplus] == 1 && avolInfo->isDomain[Zminus] == 1) {
            avolInfo->isDomain[index] = 2;
            }
            if (avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Yplus] == 1 ||
            avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Yminus] == 1 ||
            avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Zplus] == 1 ||
            avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Zminus] == 1 ||
            avolInfo->isDomain[Xminus] == 1 && avolInfo->isDomain[Yplus] == 1 ||
            avolInfo->isDomain[Xminus] == 1 && avolInfo->isDomain[Yminus] == 1 ||
            avolInfo->isDomain[Xminus] == 1 && avolInfo->isDomain[Zplus] == 1 ||
            avolInfo->isDomain[Xminus] == 1 && avolInfo->isDomain[Zminus] == 1 ||
            avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Zplus] == 1 ||
            avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Zminus] == 1 ||
            avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Zplus] == 1 ||
            avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Zminus] == 1) {
            avolInfo->isDomain[index] = 2;
            }
            } else if (X == 0 || X == Xindex - 1) {
            if (avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Yminus] == 1 ||
            avolInfo->isDomain[Zplus] == 1 && avolInfo->isDomain[Zminus] == 1 ||
            avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Zplus] == 1 ||
            avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Zminus] == 1 ||
            avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Zplus] == 1 ||
            avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Zminus] == 1) {
            avolInfo->isDomain[index] = 2;
            }
            } else if (Y == 0 || Y == Yindex - 1) {
            if (avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Xminus] == 1 ||
            avolInfo->isDomain[Zplus] == 1 && avolInfo->isDomain[Zminus] == 1 ||
            avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Zplus] == 1 ||
            avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Zminus] == 1 ||
            avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Zplus] == 1 ||
            avolInfo->isDomain[Yminus] == 1 && avolInfo->isDomain[Zminus] == 1) {
            avolInfo->isDomain[index] = 2;
            }
            } else if (Z == 0 || Z == Zindex - 1) {
            if (avolInfo->isDomain[Xplus] == 1 && avolInfo->isDomain[Xminus] == 1 ||
            avolInfo->isDomain[Yplus] == 1 && avolInfo->isDomain[Yminus] == 1) {
            avolInfo->isDomain[index] = 2;
            }
            }
            }
            */
            for (Z = 0; Z < Zindex; Z++) {
              for (Y = 0; Y < Yindex; Y++) {
                for (X = 0; X < Xindex; X++) {
                  index = Z * Yindex * Xindex + Y * Xindex + X;
                  Xplus = Z * Yindex * Xindex + Y * Xindex + (X + 1);
                  Xminus = Z * Yindex * Xindex + Y * Xindex + (X - 1);
                  Yplus = Z * Yindex * Xindex + (Y + 1) * Xindex + X;
                  Yminus = Z * Yindex * Xindex + (Y - 1) * Xindex + X;
                  Zplus = (Z + 1) * Yindex * Xindex + Y * Xindex + X;
                  Zminus = (Z - 1) * Yindex * Xindex + Y * Xindex + X;
                  if (static_cast<int>(sInfo->avi->isDomain[index]) == 0) {
                    sInfo->value[index] = 0.0;
                  }
                  if ((X % 2 != 0 && Y % 2 != 0) || (Y % 2 != 0 && Z % 2 != 0) || (Z % 2 != 0 && X % 2 != 0)) {
                    if (sInfo->avi->isDomain[index] == 2) {
                      if (sInfo->avi->isDomain[Xplus] == 1 && sInfo->avi->isDomain[Xminus] == 1) {
                        sInfo->value[index] = (sInfo->value[Xplus] + sInfo->value[Xminus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Yplus] == 1 && sInfo->avi->isDomain[Yminus] == 1) {
                        sInfo->value[index] = (sInfo->value[Yplus] + sInfo->value[Yminus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Zplus] == 1 && sInfo->avi->isDomain[Zminus] == 1) {
                        sInfo->value[index] = (sInfo->value[Zplus] + sInfo->value[Zminus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Xplus] == 1 && sInfo->avi->isDomain[Yplus] == 1) {
                        sInfo->value[index] = (sInfo->value[Xplus] + sInfo->value[Yplus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Xplus] == 1 && sInfo->avi->isDomain[Yminus] == 1) {
                        sInfo->value[index] = (sInfo->value[Xplus] + sInfo->value[Yminus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Xplus] == 1 && sInfo->avi->isDomain[Zplus] == 1) {
                        sInfo->value[index] = (sInfo->value[Xplus] + sInfo->value[Zplus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Xplus] == 1 && sInfo->avi->isDomain[Zminus] == 1) {
                        sInfo->value[index] = (sInfo->value[Xplus] + sInfo->value[Zminus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Xminus] == 1 && sInfo->avi->isDomain[Yplus] == 1) {
                        sInfo->value[index] = (sInfo->value[Xminus] + sInfo->value[Yplus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Xminus] == 1 && sInfo->avi->isDomain[Yminus] == 1) {
                        sInfo->value[index] = (sInfo->value[Xminus] + sInfo->value[Yminus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Xminus] == 1 && sInfo->avi->isDomain[Zplus] == 1) {
                        sInfo->value[index] = (sInfo->value[Xminus] + sInfo->value[Zplus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Xminus] == 1 && sInfo->avi->isDomain[Zminus] == 1) {
                        sInfo->value[index] = (sInfo->value[Xminus] + sInfo->value[Zminus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Yplus] == 1 && sInfo->avi->isDomain[Zplus] == 1) {
                        sInfo->value[index] = (sInfo->value[Yplus] + sInfo->value[Zplus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Yplus] == 1 && sInfo->avi->isDomain[Zminus] == 1) {
                        sInfo->value[index] = (sInfo->value[Yplus] + sInfo->value[Zminus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Yminus] == 1 && sInfo->avi->isDomain[Zplus] == 1) {
                        sInfo->value[index] = (sInfo->value[Yminus] + sInfo->value[Zplus]) / 2.0;
                      } else if (sInfo->avi->isDomain[Yminus] == 1 && sInfo->avi->isDomain[Zminus] == 1) {
                        sInfo->value[index] = (sInfo->value[Yminus] + sInfo->value[Zminus]) / 2.0;
                      }
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
  }

  //if compartment of reactants and products, the reaction is flux
  /*
  for (i = 0; i < numOfReactions; i++) {
  Reaction *r = model->getReaction(i);
  for (j = 0; j < r->getNumReactants(); j++) {
  for (k = 0; k < r->getNumProducts(); k++) {

  }
  }
  }
  */
  //reaction information
  for (i = 0; i < numOfReactions; i++) {
    Reaction *r = lor->get(i);
    const KineticLaw *kl = r->getKineticLaw();
    if (kl != 0) {
      for (j = 0; j < kl->getNumLocalParameters(); j++) {
        const LocalParameter *lp = kl->getLocalParameter(j);
        variableInfo *info = new variableInfo;
        InitializeVarInfo(info);
        info->id = lp->getId().c_str();
        info->value = new double[numOfVolIndexes];
        INIT_DOUBLE(info->value, numOfVolIndexes);
        for (Z = 0; Z < Zindex; Z++) {
          for (Y = 0; Y < Yindex; Y++) {
            for (X = 0; X < Xindex; X++) {
              info->value[Z * Yindex * Xindex + Y * Xindex + X] = lp->getValue();
            }
          }
        }
        varInfoList.push_back(info);
      }

      reactionInfo *rInfo = new reactionInfo;
      rInfo->id = r->getId().c_str();
      rInfo->value = new double[numOfVolIndexes];
      INIT_DOUBLE(rInfo->value, numOfVolIndexes);
      ast = const_cast<ASTNode*>(kl->getMath());
      int tmp = 0;
      countAST(ast, tmp);
#ifdef PRINT_DEBUG
      cout << "before reaction: " << SBML_formulaToString(ast) << endl;
      cerr << "num_of_nodes: " << tmp << endl;
#endif
      rearrangeAST(ast);
#ifdef PRINT_DEBUG
      cout << "after reaction: " << SBML_formulaToString(ast) << endl;
#endif
      numOfASTNodes = 0;
      countAST(ast, numOfASTNodes);
      //cerr << "num_of_nodes: " << numOfASTNodes << endl;
      rInfo->rpInfo = new reversePolishInfo();
      rInfo->rpInfo->varList = new double*[numOfASTNodes];
      rInfo->rpInfo->deltaList = new double*[numOfASTNodes];
      rInfo->rpInfo->constList = new double*[numOfASTNodes];
      rInfo->rpInfo->opfuncList = new int[numOfASTNodes];
      rInfo->rpInfo->listNum = numOfASTNodes;
      parseAST(ast, rInfo->rpInfo->varList, rInfo->rpInfo->deltaList, rInfo->rpInfo->constList, rInfo->rpInfo->opfuncList, varInfoList, numOfASTNodes);
      //parseAST(ast, rInfo->rpInfo->varList, rInfo->rpInfo->constList, rInfo->rpInfo->opfuncList, varInfoList, numOfASTNodes);
      rInfoList.push_back(rInfo);
    }
  }
  //rate rule information
  for (i = 0; i < numOfRules; i++) {
    if (model->getRule(i)->isRate()) {
      RateRule *rrule = (RateRule*)model->getRule(i);
      reactionInfo *rInfo = new reactionInfo;
      rInfo->id = rrule->getVariable().c_str();
      rInfo->value = new double[Xdiv * Ydiv * Zdiv];
      INIT_DOUBLE(rInfo->value, Xdiv * Ydiv * Zdiv);
      ast = const_cast<ASTNode*>(rrule->getMath());
      rearrangeAST(ast);
#ifdef PRINT_DEBUG
      cout << "rate rule: " << SBML_formulaToString(ast) << endl;
#endif
      numOfASTNodes = 0;
      countAST(ast, numOfASTNodes);
      rInfo->rpInfo = new reversePolishInfo();
      rInfo->rpInfo->varList = new double*[numOfASTNodes];
      rInfo->rpInfo->deltaList = new double*[numOfASTNodes];
      rInfo->rpInfo->constList = new double*[numOfASTNodes];
      rInfo->rpInfo->opfuncList = new int[numOfASTNodes];
      rInfo->rpInfo->listNum = numOfASTNodes;
      parseAST(ast, rInfo->rpInfo->varList, rInfo->rpInfo->deltaList, rInfo->rpInfo->constList, rInfo->rpInfo->opfuncList, varInfoList, numOfASTNodes);
      //parseAST(ast, rInfo->rpInfo->varList, rInfo->rpInfo->constList, rInfo->rpInfo->opfuncList, varInfoList, numOfASTNodes);
      rInfoList.push_back(rInfo);
    }
  }
}

SpatialSimulator::~SpatialSimulator()
{
  freeVarInfo(varInfoList, t_info->value);
  freeBcOfSpeciesInfo(bcOfSpeciesInfoList);
  freeAvolInfo(avolInfoList);
  freeRInfo(rInfoList);
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

double* SpatialSimulator::getGeometry(const std::string &compartmentId, int &length)
{
  analyticVolInfo* info = searchAvolInfoByCompartment(avolInfoList, compartmentId.c_str());
  length = numOfVolIndexes;
  return info->isDomain;
}

int SpatialSimulator::getGeometryLength() const
{
return numOfVolIndexes;
}

boundaryType* SpatialSimulator::getBoundaryType(const std::string &compartmentId, int &length)
{
  analyticVolInfo* info = searchAvolInfoByCompartment(avolInfoList, compartmentId.c_str());
  length = numOfVolIndexes;
  return info->bt;
}

int* SpatialSimulator::getBoundary(const std::string &compartmentId, int &length)
{
  analyticVolInfo* info = searchAvolInfoByCompartment(avolInfoList, compartmentId.c_str());
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
  struct stat st;

  if(stat("./result", &st) != 0) {    
    MKDIR("./result");
  }

  if(stat("./result/txt", &st) != 0) {    

    MKDIR("./result/txt");
  }
  if(stat("./result/png", &st) != 0) {
    MKDIR("./result/png");
  }


  int Xindex = 2 * Xdiv - 1, Yindex = 2 * Ydiv - 1, Zindex = 2 * Zdiv - 1;
#pragma region // print current step
  //int numOfVolIndexes = Xindex * Yindex * Zindex;
  unsigned int i, j; 
  int X = 0, Y = 0, Z = 0, index = 0;

  for (i = 0; i < numOfSpecies; i++) {
    Species *s = los->get(i);
    variableInfo *sInfo = searchInfoById(varInfoList, s->getId().c_str());
    if (sInfo != 0) {
      stringstream ss;
      ss << file_num;
      string dir_txt = "", dir_png = "";
      if (dimension <= 2 || (dimension == 3 && sInfo->inVol == true)) {
        dir_txt = "./result/txt/" + s->getId();
        dir_png = "./result/png/" + s->getId();
      } else if (dimension == 3 && sInfo->inVol == false) {
        if(stat("./result/txt/3Dmembrane", &st) != 0) {
          MKDIR("./result/txt/3Dmembrane");
        }
        if(stat("./result/png/3Dmembrane", &st) != 0) {
          MKDIR("./result/png/3Dmembrane");
        }
        dir_txt = "./result/txt/3Dmembrane/" + s->getId();
        dir_png = "./result/png/3Dmembrane/" + s->getId();
      }
      if(stat(dir_txt.c_str(), &st) != 0) {
        MKDIR(dir_txt.c_str());						
      }
      if(stat(dir_png.c_str(), &st) != 0) {
        MKDIR(dir_png.c_str());
      }
      string filename = "";

      if (dimension <= 2 || (dimension == 3 && sInfo->inVol == true)) {
        //filename =  "/Users/matsui/grad/SBMLSimulator/result/txt/" + s->getId() + "/" + ss.str() + ".txt";
        filename =  dir_txt + "/" + ss.str() + ".txt";
      } else if (dimension == 3 && sInfo->inVol == false) {
        //filename =  "/Users/matsui/grad/SBMLSimulator/result/txt/membrane/" + s->getId() + "/" + ss.str() + ".txt";
        filename =  dir_txt + "/" + ss.str() + ".txt";
      }
      ofstream ofs;
      ofs.open(filename.c_str());
      switch(dimension) {
      case 1:
        for (X = 0; X < Xindex; X += 2) {
          ofs << xInfo->value[X] << " " << 0 << " " << sInfo->value[X] << endl;
        }
        break;
      case 2:
        if (sInfo->avi->isVol == true) {
          for (Y = 0; Y < Yindex; Y += 2) {
            for (X = 0; X < Xindex; X += 2) {
              index = Y * Xindex + X;
              ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << sInfo->value[index] << endl;
            }
            ofs << endl;
          }
        } else if (sInfo->avi->isVol == false) {
          for (Y = 0; Y < Yindex; Y++) {
            for (X = 0; X < Xindex; X++) {
              index = Y * Xindex + X;
              ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << sInfo->value[index] << endl;
            }
            ofs << endl;
          }
        }
        break;
      case 3:
        if (sInfo->avi->isVol == true) {
          for (Z = 0; Z < Zindex; Z += 2) {
            for (Y = 0; Y < Yindex; Y += 2) {
              for (X = 0; X < Xindex; X += 2) {
                index = Z * Yindex * Xindex + Y * Xindex + X;
                ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << zInfo->value[index] << " " << sInfo->value[index] << endl;
              }
              if (Y != Yindex - 1) {
                ofs << endl;
              }
            }
            /*
            if (Z != Zindex - 1) {
            ofs << endl;
            }
            */
          }
        } else if (sInfo->avi->isVol == false) {
          bool flag = false;
          for (Z = 0; Z < Zindex; Z++) {
            for (Y = 0; Y < Yindex; Y++) {
              for (X = 0; X < Xindex; X++) {
                index = Z * Yindex * Xindex + Y * Xindex + X;
                if (sInfo->avi->isDomain[index] > 0) {
                  ofs << xInfo->value[index] << " " << yInfo->value[index] << " " << zInfo->value[index] << " " << sInfo->value[index] << endl;
                  flag = true;
                }
              }
              if (Y != Yindex - 1 && flag == true) {
                ofs << endl;
              }
              flag = false;
            }
          }
        }
        break;
      default:
        break;
      }
      ofs.close();
      //output png
      if (dimension != 3) {
        fprintf(gp, "set ticslevel 0\n");
        fprintf(gp, "set size square\n");
        fprintf(gp, "unset key\n");
        fprintf(gp, "set pm3d map corners2color c1\n");
        //fprintf(gp, "set pm3d map\n");
        fprintf(gp, "set mxtics %d\n", static_cast<int>((deltaX * (Xdiv - 1)) / 10));
        fprintf(gp, "set mxtics %d\n", static_cast<int>((deltaY * (Ydiv - 1)) / 10));
        //fprintf(gp, "set grid lw 0.5\n");
        fprintf(gp, "set cbrange[0.0:3.0]\n");
        fprintf(gp, "set palette defined (0 \"dark-blue\", 2 \"blue\", 4 \"green\", 8 \"yellow\", 10 \"red\")\n");
        fprintf(gp, "set xlabel \"x\"\n");
        fprintf(gp, "set ylabel \"y\"\n");
        fprintf(gp, "set title \"t=%lf\"\n", sim_time);
        fprintf(gp, "set terminal png truecolor\n");
#ifdef WIN32
        //fprintf(gp, "set output \"/dev/null\"\n");
#else
        fprintf(gp, "set output \"/dev/null\"\n");
#endif
        //fprintf(gp, "set output \"/Users/matsui/grad/SBMLSimulator/result/png/%s/%d.png\"\n", s->getId().c_str(), file_num);
        fprintf(gp, "splot \"%s\" with image failsafe\n", filename.c_str());
        if (sInfo->avi->isVol == true) {
          // 							fprintf(gp, "set output \"/Users/matsui/grad/SBMLSimulator/result/png/%s/%d.png\"\n", s->getId().c_str(), file_num);
          // 							fprintf(gp, "replot \"/Users/matsui/grad/SBMLSimulator/result/txt/geo_sum.txt\" with rgbalpha failsafe t\"\"\n");
          fprintf(gp, "set output \"./result/png/%s/%d.png\"\n", s->getId().c_str(), file_num);
          fprintf(gp, "replot \"./result/txt/geo_sum.txt\" with rgbalpha failsafe t\"\"\n");
        } else {
          for (j = 0; j < memList.size(); j++) {
            if (strcmp(sInfo->avi->domainTypeId, memList[j]) != 0) {
              if (j != memList.size() - 1 && memList.size() != 2) {
                fprintf(gp, "set output \"/dev/null\"\n");
                // 										fprintf(gp, "replot \"/Users/matsui/grad/SBMLSimulator/result/txt/%s.txt\" with rgbalpha failsafe t\"\"\n", memList[j]);
                fprintf(gp, "replot \"./result/txt/%s.txt\" with rgbalpha failsafe t\"\"\n", memList[j]);
              } else {
                // 										fprintf(gp, "set output \"/Users/matsui/grad/SBMLSimulator/result/png/%s/%d.png\"\n", s->getId().c_str(), file_num);
                // 										fprintf(gp, "replot \"/Users/matsui/grad/SBMLSimulator/result/txt/%s.txt\" with rgbalpha failsafe t\"\"\n", memList[j]);
                fprintf(gp, "set output \"./result/png/%s/%d.png\"\n", s->getId().c_str(), file_num);
                fprintf(gp, "replot \"./result/txt/%s.txt\" with rgbalpha failsafe t\"\"\n", memList[j]);
              }
            }
          }
        }
      } else if (dimension == 3 && sInfo->inVol == false) {
        fprintf(gp, "set ticslevel 0\n");
        fprintf(gp, "set size square\n");
        fprintf(gp, "unset key\n");
        fprintf(gp, "set cbrange[0.0:6.5]\n");
        fprintf(gp, "set palette defined (0 \"dark-blue\", 2 \"blue\", 4 \"green\", 8 \"yellow\", 10 \"red\")\n");
        fprintf(gp, "set xlabel \"x\"\n");
        fprintf(gp, "set ylabel \"y\"\n");
        fprintf(gp, "set zlabel \"z\"\n");
        fprintf(gp, "set title \"t=%lf\"\n", sim_time);
        fprintf(gp, "set terminal png truecolor\n");
        //fprintf(gp, "set output \"/dev/null\"\n");
        //fprintf(gp, "set output \"/Users/matsui/grad/SBMLSimulator/result/png/%s/%d.png\"\n", s->getId().c_str(), file_num);
        // 						fprintf(gp, "set output \"/Users/matsui/grad/SBMLSimulator/result/png/3Dmembrane/%s/%d.png\"\n", s->getId().c_str(), file_num);
        fprintf(gp, "set output \"./result/png/3Dmembrane/%s/%d.png\"\n", s->getId().c_str(), file_num);
        //fprintf(gp, "splot \"%s\" with image failsafe\n", filename.c_str());
        fprintf(gp, "splot \"%s\"with points palette\n", filename.c_str());
      }
    }
  }
  file_num++;


#pragma endregion
}

void SpatialSimulator::performStep(double t, double dt)
{
  unsigned int /* i,*/ j , m;
  int X = 0, Y = 0, Z = 0, index = 0;

  int Xplus = 0, Xminus = 0, Yplus = 0, Yminus = 0, Zplus = 0, Zminus = 0;
#pragma region // integration step 
  // 		stringstream ss;
  // 		ss << t;
  // 		string filename =  "./speed/" + ss.str() + ".txt";
  // 		ofstream ofs;
  // 		ofs.open(filename.c_str());
  //calculation
  //runge-kutta
  for (m = 0; m < 4; m++) {
    //reaction
    int i;
#pragma omp parallel for      
    for (i = 0; i < (int)numOfReactions; i++) {
      Reaction *r = model->getReaction(i);
      //reactants
      for (j = 0; j < r->getNumReactants(); j++) {
        SpeciesReference *sr = r->getReactant(j);
        variableInfo *sInfo = searchInfoById(varInfoList, sr->getSpecies().c_str());
        Species *s = model->getSpecies(sr->getSpecies());
        if (!s->isSetConstant() || !s->getConstant()) {
          reversePolishRK(rInfoList[i], sInfo, Xindex, Yindex, Zindex, dt, m, reactants, sInfo->inVol, sInfo->dci);
          if (m == 0) {
            //ofs << "test" << endl;
            //cout << rInfoList[i]->rpInfo->listNum << " " << ((sim_end - sim_start) / static_cast<double>(CLOCKS_PER_SEC)) << endl;
            //ofs << rInfoList[i]->rpInfo->listNum << " " << ((ast_end - ast_start) / static_cast<double>(CLOCKS_PER_SEC)) << endl;
          }
        }
      }//the end of reactants
      //products
      for (j = 0; j < r->getNumProducts(); j++) {
        SpeciesReference *sr = r->getProduct(j);
        variableInfo *sInfo = searchInfoById(varInfoList, sr->getSpecies().c_str());
        Species *s = model->getSpecies(sr->getSpecies());
        if (!s->isSetConstant() || !s->getConstant()) {
          reversePolishRK(rInfoList[i], sInfo, Xindex, Yindex, Zindex, dt, m, products, sInfo->inVol, sInfo->dci);
        }
      }//end of products
    }//end of reaction
    //rate rule
    for (i = 0; i < (int) numOfRules; i++) {
      if (model->getRule(i)->isRate()) {
        RateRule *rrule = (RateRule*)model->getRule(i);
        variableInfo *sInfo = searchInfoById(varInfoList, rrule->getVariable().c_str());
        Species *s = model->getSpecies(rrule->getVariable());
        if ((!s->isSetConstant() || !s->getConstant()) && sInfo->hasAssignmentRule == false) {
          reversePolishRK(rInfoList[i + numOfReactions], sInfo, Xindex, Yindex, Zindex, dt, m, products, sInfo->inVol, sInfo->dci);
        }
      }
    }
    //diffusion, advection
    for (i = 0; i < (int)numOfSpecies; i++) {
      variableInfo *sInfo = searchInfoById(varInfoList, los->get(i)->getId().c_str());
      double rk[4] = {0, 0.5, 0.5, 1.0};
      //diffusion
      if (sInfo->dci != 0) {
        for (Z = 0; Z < Zindex; Z += 2) {
          for (Y = 0; Y < Yindex; Y += 2) {
            for (X = 0; X < Xindex; X += 2) {
              index = Z * Yindex * Xindex + Y * Xindex + X;
              Xplus = Z * Yindex * Xindex + Y * Xindex + (X + 2);
              Xminus = Z * Yindex * Xindex + Y * Xindex + (X - 2);
              Yplus = Z * Yindex * Xindex + (Y + 2) * Xindex + X;
              Yminus = Z * Yindex * Xindex + (Y - 2) * Xindex + X;
              Zplus = (Z + 2) * Yindex * Xindex + Y * Xindex + X;
              Zminus = (Z - 2) * Yindex * Xindex + Y * Xindex + X;
              if (sInfo->avi->isDomain[index] == 1 && sInfo->avi->isBoundary[index] == 0) {
                double* val = sInfo->value;
                double *d = sInfo->delta;
                if (m == 0) {
                  if (sInfo->dci[0] != 0) {//x-diffusion
                    sInfo->delta[m * numOfVolIndexes + index] += sInfo->dci[0]->value[index]
                    * (val[Xplus] - 2 * val[index] + val[Xminus]) / pow(deltaX, 2);
                  }
                  if (sInfo->dci[1] != 0) {//y-diffusion
                    sInfo->delta[m * numOfVolIndexes + index] += sInfo->dci[1]->value[index]
                    * (val[Yplus] - 2 * val[index] + val[Yminus]) / pow(deltaY, 2);
                  }
                  if (sInfo->dci[2] != 0) {//z-diffusion
                    sInfo->delta[m * numOfVolIndexes + index] += sInfo->dci[2]->value[index]
                    * (val[Zplus] - 2 * val[index] + val[Zminus]) / pow(deltaZ, 2);
                  }
                } else {
                  if (sInfo->dci[0] != 0) {//x-diffusion
                    sInfo->delta[m * numOfVolIndexes + index] += sInfo->dci[0]->value[index]
                    * ((val[Xplus] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + Xplus])
                      - 2 * (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])
                      + (val[Xminus] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + Xminus])) / pow(deltaX, 2);
                  }
                  if (sInfo->dci[1] != 0) {//y-diffusion
                    sInfo->delta[m * numOfVolIndexes + index] += sInfo->dci[1]->value[index]
                    * ((val[Yplus] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + Yplus])
                      - 2 * (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])
                      + (val[Yminus] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + Yminus])) / pow(deltaY, 2);
                  }
                  if (sInfo->dci[2] != 0) {//z-diffusion
                    sInfo->delta[m * numOfVolIndexes + index] += sInfo->dci[2]->value[index]
                    * ((val[Zplus] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + Zplus])
                      - 2 * (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])
                      + (val[Zminus] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + Zminus])) / pow(deltaZ, 2);
                  }
                }
              }
            }
          }
        }
        //boundary condition
        //if no flux, nuemann(0)
        for (Z = 0; Z < Zindex; Z += 2) {
          for (Y = 0; Y < Yindex; Y += 2) {
            for (X = 0; X < Xindex; X += 2) {
              index = Z * Yindex * Xindex + Y * Xindex + X;
              if (sInfo->avi->isDomain[index] == 1 && sInfo->avi->isBoundary[index] == 1) {
                switch(sInfo->avi->bt[index]) {
                case Xp:
                  sInfo->value[index] = sInfo->value[Z * Yindex * Xindex + Y * Xindex + (X - 2)]
                  + dt * sInfo->delta[m * numOfVolIndexes + (Z * Yindex * Xindex + Y * Xindex + (X - 2))];
                  break;
                case Xm:
                  sInfo->value[index] = sInfo->value[Z * Yindex * Xindex + Y * Xindex + (X + 2)]
                  + dt * sInfo->delta[m * numOfVolIndexes + (Z * Yindex * Xindex + Y * Xindex + (X + 2))];
                  break;
                case Yp:
                  sInfo->value[index] = sInfo->value[Z * Yindex * Xindex + (Y - 2) * Xindex + X]
                  + dt * sInfo->delta[m * numOfVolIndexes + (Z * Yindex * Xindex + (Y - 2) * Xindex + X)];
                  break;
                case Ym:
                  sInfo->value[index] = sInfo->value[Z * Yindex * Xindex + (Y + 2) * Xindex + X]
                  + dt * sInfo->delta[m * numOfVolIndexes + (Z * Yindex * Xindex + (Y + 2) * Xindex + X)];
                  break;
                case Zp:
                  sInfo->value[index] = sInfo->value[(Z - 2) * Yindex * Xindex + Y * Xindex + X]
                  + dt * sInfo->delta[m * numOfVolIndexes + ((Z - 2) * Yindex * Xindex + Y * Xindex + X)];
                  break;
                case Zm:
                  sInfo->value[index] = sInfo->value[(Z + 2) * Yindex * Xindex + Y * Xindex + X]
                  + dt * sInfo->delta[m * numOfVolIndexes + ((Z + 2) * Yindex * Xindex + Y * Xindex + X)];
                  break;
                default:
                  break;
                }
              }
            }
          }
        }
      }//end of diffusion

      //advection
      if (sInfo->aci != 0) {
        double *fx_cip = 0, *fy_cip = 0, *fz_cip = 0, *tmp = 0;
        double* val = sInfo->value;
        double *d = sInfo->delta;
        if (t == 0 && m == 0) {
          if (sInfo->dcip == 0) {
            sInfo->dcip = new derivativeCIP;
            sInfo->dcip->fx = 0;
            sInfo->dcip->fx_next = 0;
            sInfo->dcip->fy = 0;
            sInfo->dcip->fy_next = 0;
            sInfo->dcip->fz = 0;
            sInfo->dcip->fz_next = 0;
          }
          if (dimension >= 1 && sInfo->dcip->fx == 0) {
            sInfo->dcip->fx = new double[numOfVolIndexes];
            sInfo->dcip->fx_next = new double[numOfVolIndexes];
          }
          if (dimension >= 2 && sInfo->dcip->fy == 0) {
            sInfo->dcip->fy = new double[numOfVolIndexes];
            sInfo->dcip->fy_next = new double[numOfVolIndexes];
          }
          if (dimension >= 3 && sInfo->dcip->fz == 0) {
            sInfo->dcip->fz = new double[numOfVolIndexes];
            sInfo->dcip->fz_next = new double[numOfVolIndexes];
          }
        }

        //境界以外の微分係数
        if (!(t == 0 && m == 0)) {
          // 						for (Z = 0; Z < Zdiv; Z++) {
          // 							for (Y = 0; Y < Ydiv; Y++) {
          // 								for (X = 0; X < Xdiv; X++) {
          // 									index = Z * Ydiv * Xdiv + Y * Xdiv + X;
          // 									if (sInfo->dcip->fx != 0) {
          // 										sInfo->dcip->fx[index] = sInfo->dcip->fx_next[index];
          // 									}
          // 									if (sInfo->dcip->fy != 0) {
          // 										sInfo->dcip->fy[index] = sInfo->dcip->fy_next[index];
          // 									}
          // 									if (sInfo->dcip->fz != 0) {
          // 										sInfo->dcip->fz[index] = sInfo->dcip->fz_next[index];
          // 									}
          // 								}
          // 							}
          // 						}
          //fx
          tmp = sInfo->dcip->fx;
          sInfo->dcip->fx = sInfo->dcip->fx_next;
          sInfo->dcip->fx_next = tmp;
          //fy
          tmp = sInfo->dcip->fy;
          sInfo->dcip->fy = sInfo->dcip->fy_next;
          sInfo->dcip->fy_next = tmp;
          //fz
          tmp = sInfo->dcip->fz;
          sInfo->dcip->fz = sInfo->dcip->fz_next;
          sInfo->dcip->fz_next = tmp;

          fx_cip = sInfo->dcip->fx;
          fy_cip = sInfo->dcip->fy;
          fz_cip = sInfo->dcip->fz;
        } else {//t == 0 && m == 0
          fx_cip = sInfo->dcip->fx;
          fy_cip = sInfo->dcip->fy;
          fz_cip = sInfo->dcip->fz;
          for (Z = 0; Z < Zindex; Z += 2) {
            for (Y = 0; Y < Yindex; Y += 2) {
              for (X = 0; X < Xindex; X += 2) {
                index = Z * Yindex * Xindex + Y * Xindex + X;
                if (sInfo->avi->isDomain[index] == 1 && sInfo->avi->isBoundary[index] == 0) {
                  Xplus = Z * Yindex * Xindex + Y * Xindex + (X + 2);
                  Xminus = Z * Yindex * Xindex + Y * Xindex + (X - 2);
                  Yplus = Z * Yindex * Xindex + (Y + 2) * Xindex + X;
                  Yminus = Z * Yindex * Xindex + (Y - 2) * Xindex + X;
                  Zplus = (Z + 2) * Yindex * Xindex + Y * Xindex + X;
                  Zminus = (Z - 2) * Yindex * Xindex + Y * Xindex + X;
                  if (dimension >= 1) {
                    fx_cip[index] = (val[Xplus] - val[Xminus]) / (2.0 * deltaX);
                  }
                  if (dimension >= 2) {
                    fy_cip[index] = (val[Yplus] - val[Yminus]) / (2.0 * deltaY);
                  }
                  if (dimension >= 3) {
                    fy_cip[index] = (val[Zplus] - val[Zminus]) / (2.0 * deltaZ);
                  }
                }
              }
            }
          }
        }

        //境界での微分係数
        for (Z = 0; Z < Zindex; Z += 2) {
          for (Y = 0; Y < Yindex; Y += 2) {
            for (X = 0; X < Xindex; X += 2) {
              index = Z * Yindex * Xindex + Y * Xindex + X;
              if (sInfo->avi->isDomain[index] == 1 && sInfo->avi->isBoundary[index] == 1) {
                switch(sInfo->avi->bt[index]) {
                case Xp:
                  fx_cip[index] = fx_cip[Z * Yindex * Xindex + Y * Xindex + (X - 2)];
                  break;
                case Xm:
                  fx_cip[index] = fx_cip[Z * Yindex * Xindex + Y * Xindex + (X + 2)];
                  break;
                case Yp:
                  fy_cip[index] = fy_cip[Z * Yindex * Xindex + (Y - 2) * Xindex + X];
                  break;
                case Ym:
                  fy_cip[index] = fy_cip[Z * Yindex * Xindex + (Y + 2) * Xindex + X];
                  break;
                case Zp:
                  fz_cip[index] = fz_cip[(Z - 2) * Yindex * Xindex + Y * Xindex + X];
                  break;
                case Zm:
                  fz_cip[index] = fz_cip[(Z + 2) * Yindex * Xindex + Y * Xindex + X];
                  break;
                default:
                  break;
                }
              }
            }
          }
        }
        //cip
        int index_Xm = 0, index_Ym = 0, index_XmYm = 0;
        double ux = 0.0, uy = 0.0, Xsign = 0.0, Ysign = 0.0, XX = 0.0, YY = 0.0;
        double A = 0.0, B = 0.0, C = 0.0, D = 0.0, E = 0.0, F = 0.0, G = 0.0;

        switch (dimension) {
        case 1:
          for (X = 0; X < Xindex; X += 2) {
            if (sInfo->avi->isDomain[X] == 1 && sInfo->avi->isBoundary[X] == 0) {
              if (sInfo->aci[0] != 0) {
                ux = sInfo->aci[0]->value[index];
                Xsign = copysign(1.0, ux);
              }
            }
            index_Xm = X - 2 * static_cast<int>(Xsign);
            if (m == 0) {
              A = (Xsign * 2.0 * (val[index_Xm] - val[index]) + (fx_cip[X] + fx_cip[index_Xm]) * deltaX ) / pow(deltaX, 3);
              B = (3.0 * (val[index_Xm] - val[X]) + Xsign * (2.0 * fx_cip[X] + fx_cip[index_Xm]) * deltaX) / pow(deltaX, 2);
            } else {
              A = (Xsign * 2.0 * ((val[index_Xm] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index_Xm]) - (val[X] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + X])) + (fx_cip[X] + fx_cip[index_Xm]) * deltaX ) / pow(deltaX, 3);
              B = (3.0 * ((val[index_Xm] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index_Xm]) - (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])) + Xsign * (2.0 * fx_cip[X] + fx_cip[index_Xm]) * deltaX) / pow(deltaX, 2);
            }
            XX = -ux * dt;
            d[m * numOfVolIndexes + X] = A * pow(XX, 3) + B * pow(XX, 2) + fx_cip[X] * XX;
            sInfo->dcip->fx_next[X] = 3.0 * A * pow(XX, 2) + 2.0 * B * XX + fx_cip[X];
          }
          break;
        case 2:
          for (Y = 0; Y < Yindex; Y += 2) {
            for (X = 0; X < Xindex; X += 2) {
              index = Y * Xindex + X;
              if (sInfo->avi->isDomain[index] == 1 && sInfo->avi->isBoundary[index] == 0) {
                if (sInfo->aci[0] != 0) {
                  ux = sInfo->aci[0]->value[index];
                  Xsign = copysign(1.0, ux);
                }
                if (sInfo->aci[1] != 0) {
                  uy = sInfo->aci[1]->value[index];
                  Ysign = copysign(1.0, uy);
                }

                index_Xm = Y * Xindex + (X - 2 * static_cast<int>(Xsign));
                index_Ym = ((Y - 2 * static_cast<int>(Ysign))) * Xindex + X;
                index_XmYm = ((Y - 2 * static_cast<int>(Ysign))) * Xindex + (X - static_cast<int>(Xsign));
                if (m == 0) {
                  A = ((fx_cip[index] + fx_cip[index_Xm]) * deltaX - Xsign * 2.0 * (val[index] - val[index_Xm])) / pow(deltaX, 3);
                  B = ((fy_cip[index] + fy_cip[index_Ym]) * deltaY - Ysign * 2.0 * (val[index] - val[index_Ym])) / pow(deltaY, 3);
                  C = Ysign * (-1.0 * (val[index] - val[index_Ym] - val[index_Xm] + val[index_XmYm])
                    - Xsign * (fx_cip[index_Ym] - fx_cip[index]) * deltaX) / (pow(deltaX, 2) * deltaY);
                  D = Xsign * (-1.0 * (val[index] - val[index_Ym] - val[index_Xm] + val[index_XmYm])
                    - Ysign * (fy_cip[index_Xm] - fy_cip[index]) * deltaY) / (deltaX * pow(deltaY, 2));
                  E = (3.0 * (val[index_Xm] - val[index]) + Xsign * (fx_cip[index_Xm] + 2.0 * fx_cip[index]) * deltaX) / pow(deltaX, 2);
                  F = (3.0 * (val[index_Ym] - val[index]) + Ysign * (fy_cip[index_Ym] + 2.0 * fy_cip[index]) * deltaY) / pow(deltaY, 2);
                  G = Xsign * (fy_cip[index] - fy_cip[index_Xm] + C * pow(deltaX, 2)) / deltaX;
                } else {
                  A = ((fx_cip[index] + fx_cip[index_Xm]) * deltaX - Xsign * 2.0 * ((val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index]) - (val[index_Xm] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index_Xm]))) / pow(deltaX, 3);
                  B = ((fy_cip[index] + fy_cip[index_Ym]) * deltaY - Ysign * 2.0 * ((val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index]) - (val[index_Ym] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index_Ym]))) / pow(deltaY, 3);
                  C = Ysign * (-1.0 * ((val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index]) - (val[index_Ym] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index_Ym]) - (val[index_Xm] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index_Xm]) + (val[index_XmYm] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index_XmYm])) - Xsign * (fx_cip[index_Ym] - fx_cip[index]) * deltaX) / (pow(deltaX, 2) * deltaY);
                  D = Xsign * (-1.0 * ((val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index]) - (val[index_Ym] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index_Ym]) - (val[index_Xm] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index_Xm]) + (val[index_XmYm] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index_XmYm])) - Ysign * (fy_cip[index_Xm] - fy_cip[index]) * deltaY) / (deltaX * pow(deltaY, 2));
                  E = (3.0 * ((val[index_Xm] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index_Xm]) - (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])) + Xsign * (fx_cip[index_Xm] + 2.0 * fx_cip[index]) * deltaX) / pow(deltaX, 2);
                  F = (3.0 * ((val[index_Ym] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index_Ym]) - (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])) + Ysign * (fy_cip[index_Ym] + 2.0 * fy_cip[index]) * deltaY)
                    / pow(deltaY, 2);
                  G = Xsign * (fy_cip[index] - fy_cip[index_Xm] + C * pow(deltaX, 2)) / deltaX;
                }
                XX = -ux * dt;
                YY = -uy * dt;
                //cout << XX << endl;
                d[m * numOfVolIndexes + index] = A * pow(XX, 3) + B * pow(YY, 3) + C * pow(XX, 2) * YY + D * XX * pow(YY, 2) + E * pow(XX, 2) + F * pow(YY, 2) + G * XX * YY + fx_cip[index] * XX + fy_cip[index] * YY;
                sInfo->dcip->fx_next[index] = 3.0 * A * pow(XX, 2) + 2.0 * C * XX * YY + D * pow(YY, 2) + 2.0 * E * XX + G * YY + fx_cip[index];
                sInfo->dcip->fy_next[index] = 3.0 * B * pow(YY, 2) + C * pow(XX, 2) + 2.0 * D * XX * YY + 2.0 * F * YY + G * XX + fy_cip[index];
              }
            }
          }
          break;
        case 3:
          break;
        default:
          break;
        }
        //boundary condition
      }//end of advection
      //membrane flux
      //under construction

    }
  }//end of runge-kutta
  //ofs.close();

#pragma endregion

}

void SpatialSimulator::updateValues(double dt)
{

#pragma region     // update values
  //update values
#pragma omp parallel for 
  for (int i = 0; i < (int)numOfSpecies; i++) 
  {
    Species *s = los->get(i);
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
  unsigned int j;
  int X = 0, Y = 0, Z = 0, index = 0;

  int Xplus = 0, Xminus = 0, Yplus = 0, Yminus = 0, Zplus = 0, Zminus = 0;
  if (sInfo == NULL) return;

#pragma region // update species
  for (Z = 0; Z < Zindex; Z++) 
  {
    for (Y = 0; Y < Yindex; Y++) 
    {
      for (X = 0; X < Xindex; X++) 
      {
        if ((sInfo->inVol == true && X % 2 == 0 && Y % 2 == 0 && Z % 2 == 0) ||
          (sInfo->inVol == false && (Z * Yindex * Xindex + Y * Xindex + X) % 2 != 0)) 
        {
          index = Z * Yindex * Xindex + Y * Xindex + X;
          //update values for the next time
          if (static_cast<int>(sInfo->avi->isDomain[index] == 1) && sInfo->avi->isBoundary[index] == 0) 
          {
            sInfo->value[index] += dt * (sInfo->delta[index] + 2.0 * sInfo->delta[numOfVolIndexes + index] + 2.0 * sInfo->delta[2 * numOfVolIndexes + index] + sInfo->delta[3 * numOfVolIndexes + index]) / 6.0;
            //reset delta
            for (j = 0; j < 4; j++) 
            {
              sInfo->delta[j * numOfVolIndexes + index] = 0.0;
            }
          }
        }
      }
    }
  }
#pragma endregion
#pragma region // update boundary in domain
  //boundary condition
  //in domain
  //if no flux, nuemann(0)
  for (Z = 0; Z < Zindex; Z += 2) {
    for (Y = 0; Y < Yindex; Y += 2) {
      for (X = 0; X < Xindex; X += 2) {
        index = Z * Yindex * Xindex + Y * Xindex + X;
        if (sInfo->avi->isDomain[index] == 1 && sInfo->avi->isBoundary[index] == 1) {
          switch(sInfo->avi->bt[index]) {
          case Xp:
            sInfo->value[index] = sInfo->value[Z * Yindex * Xindex + Y * Xindex + (X - 2)];
            break;
          case Xm:
            sInfo->value[index] = sInfo->value[Z * Yindex * Xindex + Y * Xindex + (X + 2)];
            break;
          case Yp:
            sInfo->value[index] = sInfo->value[Z * Yindex * Xindex + (Y - 2) * Xindex + X];
            break;
          case Ym:
            sInfo->value[index] = sInfo->value[Z * Yindex * Xindex + (Y + 2) * Xindex + X];
            // 									if (X == 88 && Y == 24 && strcmp(sInfo->id, "u") == 0) {
            // 										cout << sInfo->value[index] << " " << sInfo->value[Z * Yindex * Xindex + (Y + 2) * Xindex + X] << " " << sInfo->value[Z * Yindex * Xindex + (Y + 4) * Xindex + X];
            // 									}
            break;
          case Zp:
            sInfo->value[index] = sInfo->value[(Z - 2) * Yindex * Xindex + Y * Xindex + X];
            break;
          case Zm:
            sInfo->value[index] = sInfo->value[(Z + 2) * Yindex * Xindex + Y * Xindex + X];
            break;
          default:
            break;
          }
        }
      }
    }
  }
#pragma endregion
#pragma region // update boundary at the edge
  //boundary
  //at the edge of simulation area
  //parameter
  //Xp, Xm
  bcOfSpeciesInfo *bcos = searchBCOSInfoBySpeId(bcOfSpeciesInfoList, sInfo->id);


  bool ypFlux = bcos->bcYpType == "Flux";
  bool ypValue = bcos->bcYpType == "Value";
  bool ymFlux = bcos->bcYmType == "Flux";
  bool ymValue = bcos->bcYmType == "Value";

#pragma region // one dimension
  if (dimension >= 1) 
  {
    for (Z = 0; Z < Zindex; Z += 2) 
    {
      for (Y = 0; Y < Yindex; Y += 2) 
      {
        Xplus = Z * Yindex * Xindex + Y * Xindex + Xindex - 1;
        Xminus = Z * Yindex * Xindex + Y * Xindex;

        if (Xplus >= numOfVolIndexes  || Xminus >=numOfVolIndexes)
          continue;

        if (static_cast<int>(sInfo->avi->isDomain[Xplus]) == 1) 
        {
          if (bcos->bcXpType == "Flux") 
          {
            sInfo->value[Xplus] = bcos->bcXp->value[Xplus] * deltaX + sInfo->value[Xplus - 2];
          } 
          else if (bcos->bcXpType == "Value") 
          {
            sInfo->value[Xplus] = bcos->bcXp->value[Xplus];
          }
        }
        if (static_cast<int>(sInfo->avi->isDomain[Xminus]) == 1) 
        {
          if (bcos->bcXmType == "Flux") 
          {
            sInfo->value[Xminus] = bcos->bcXm->value[Xminus] * deltaX + sInfo->value[Xminus + 2];
          } 
          else if (bcos->bcXmType == "Value") 
          {
            sInfo->value[Xminus] = bcos->bcXm->value[Xminus];
          }
        }
      }
    }
  }
#pragma endregion
#pragma region // two dimensions
  if (dimension >= 2) 
  {
    //Yp, Ym

    for (Z = 0; Z < Zindex; Z += 2) 
    {
      for (X = 0; X < Xindex; X += 2) 
      {
        Yplus = Z * Yindex * Xindex + (Yindex - 1) * Xindex + X;
        Yminus = Z * Yindex * Xindex + X;

        if (Yplus >= numOfVolIndexes  || Yminus >=numOfVolIndexes)
          continue;

        bool isInYPlus = static_cast<int>(sInfo->avi->isDomain[Yplus]) == 1;
        bool isInYMinus = static_cast<int>(sInfo->avi->isDomain[Yminus]) == 1 && sInfo->avi->isBoundary[Yminus] == 1;


        if (isInYPlus) 
        {
          if (ypFlux) 
          {
            sInfo->value[Yplus] = bcos->bcYp->value[Yplus] * deltaY + sInfo->value[Yplus - 2];
          } 
          else if (ypValue) 
          {
            sInfo->value[Yplus] = bcos->bcYp->value[Yplus];
          }
        }
        if (isInYMinus) 
        {
          if (ymFlux) 
          {
            sInfo->value[Yminus] = bcos->bcYm->value[Yminus] * deltaY + sInfo->value[Yminus + 2];
          } 
          else if (ymValue) 
          {
            sInfo->value[Yminus] = bcos->bcYm->value[Yminus];
          }
        }
      }
    }

  }
#pragma endregion
#pragma region // three dimensions
  //Zp, Zm
  if (dimension == 3) 
  {
    for (Y = 0; Y < Yindex; Y += 2) 
    {
      for (X = 0; X < Xindex; X += 2) 
      {
        Zplus = (Zindex - 1) * Yindex * Xindex + Y * Xindex  + X;
        Zminus = Y * Xindex + Xindex+ X;

        if (Zplus >= numOfVolIndexes  || Zminus >=numOfVolIndexes)
          continue;

        if (static_cast<int>(sInfo->avi->isDomain[Zplus]) == 1) 
        {
          if (bcos->bcZpType == "Flux") 
          {
            sInfo->value[Zplus] = bcos->bcZp->value[Zplus] * deltaZ + sInfo->value[Zplus - 2];
          } 
          else if (bcos->bcZpType == "Value") 
          {
            sInfo->value[Zplus] = bcos->bcZp->value[Zplus];
          }
        }
        if (static_cast<int>(sInfo->avi->isDomain[Zminus]) == 1) 
        {
          if (bcos->bcZmType == "Flux") 
          {
            sInfo->value[Zminus] = bcos->bcZm->value[Zminus] * deltaZ + sInfo->value[Zminus + 2];
          } 
          else if (bcos->bcZmType == "Value") 
          {
            sInfo->value[Zminus] = bcos->bcZm->value[Zminus];
          }
        }
      }
    }
  }
#pragma endregion
#pragma endregion

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
        if (static_cast<int>(info->avi->isDomain[index]) == 0) {
          info->value[index] = 0.0;
        }
      }
    }
  }
}

void SpatialSimulator::updateAssignmentRules()
{
  //assignment rule
  for (unsigned int i = 0; i < varInfoList.size(); i++) {
    variableInfo *info = varInfoList[i];
    if (info->hasAssignmentRule) {
      reversePolishInitial(info->rpInfo->varList, info->rpInfo->constList, info->rpInfo->opfuncList, info->value, info->rpInfo->listNum, Xindex, Yindex, Zindex, info->inVol);
      deleteValuesOutsideDomain(info);      
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
