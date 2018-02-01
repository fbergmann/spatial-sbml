#include "sbml/SBMLTypes.h"
#include "sbml/extension/SBMLExtensionRegistry.h"
#include "sbml/packages/spatial/common/SpatialExtensionTypes.h"
#include "sbml/packages/spatial/extension/SpatialModelPlugin.h"
#include "sbml/packages/spatial/extension/SpatialExtension.h"
#include <vector>
#include <cmath>
#include "mystruct.h"
#include "searchFunction.h"

#ifndef M_E
#define M_E        2.71828182845904523536
#endif
#ifndef  M_PI
#define M_PI       3.14159265358979323846
#endif

void parseAST(ASTNode *ast, reversePolishInfo *rpInfo, vector<variableInfo*> &varInfoList, int index_max, vector<double*> &freeConstList)
{
	static int index = 0;
	unsigned int i;
	for (i = 0; i < ast->getNumChildren(); i++) {
		parseAST(ast->getChild(i), rpInfo, varInfoList, index_max, freeConstList);
	}
	if (ast->isFunction() || ast->isOperator() || ast->isRelational() || ast->isLogical()) {
		//ast is function, operator, relational or logical
		rpInfo->varList[index] = 0;
		rpInfo->constList[index] = 0;
		rpInfo->opfuncList[index] = ast->getType();
		if (rpInfo->deltaList != 0) rpInfo->deltaList[index] = 0;
	} else if (ast->isReal()) {//ast is real number
		rpInfo->varList[index] = 0;
		rpInfo->constList[index] = new double(ast->getReal());
		rpInfo->opfuncList[index] = 0;
		if (rpInfo->deltaList != 0) rpInfo->deltaList[index] = 0;
		freeConstList.push_back(rpInfo->constList[index]);
	} else if (ast->isInteger()) {//ast is integer
		rpInfo->varList[index] = 0;
		rpInfo->constList[index] = new double((double)ast->getInteger());
		rpInfo->opfuncList[index] = 0;
		if (rpInfo->deltaList != 0) rpInfo->deltaList[index] = 0;
		freeConstList.push_back(rpInfo->constList[index]);
	} else if (ast->isConstant()) {//ast is constant
		int type = (int)ast->getType();
		if (type == AST_CONSTANT_E) {
			rpInfo->varList[index] = 0;
			rpInfo->constList[index] = new double(M_E);
			rpInfo->opfuncList[index] = 0;
			if (rpInfo->deltaList != 0) rpInfo->deltaList[index] = 0;
			freeConstList.push_back(rpInfo->constList[index]);
		} else if (type == AST_CONSTANT_PI) {
			rpInfo->varList[index] = 0;
			rpInfo->constList[index] = new double(M_PI);
			rpInfo->opfuncList[index] = 0;
			if (rpInfo->deltaList != 0) rpInfo->deltaList[index] = 0;
			freeConstList.push_back(rpInfo->constList[index]);
		}  else if (type == AST_CONSTANT_FALSE) {
			rpInfo->varList[index] = 0;
			rpInfo->constList[index] = new double(0.0);
			rpInfo->opfuncList[index] = 0;
			if (rpInfo->deltaList != 0) rpInfo->deltaList[index] = 0;
			freeConstList.push_back(rpInfo->constList[index]);
		}  else if (type == AST_CONSTANT_TRUE) {
			rpInfo->varList[index] = 0;
			rpInfo->constList[index] = new double(1.0);
			rpInfo->opfuncList[index] = 0;
			if (rpInfo->deltaList != 0) rpInfo->deltaList[index] = 0;
			freeConstList.push_back(rpInfo->constList[index]);
		}
	}
	else if (ast->isName()) {
		int type = (int)ast->getType();
		if (type == AST_NAME) {
			variableInfo *info = searchInfoById(varInfoList, ast->getName());
			if (info != 0) {
				if (info->isUniform) {
					//cout << info->id << " " << *(info->value) << endl;;
					rpInfo->varList[index] = 0;
					rpInfo->constList[index] = info->value;
					rpInfo->opfuncList[index] = 0;
					if (rpInfo->deltaList != 0) rpInfo->deltaList[index] = 0;
				} else {
					rpInfo->varList[index] = info->value;
					rpInfo->constList[index] = 0;
					rpInfo->opfuncList[index] = 0;
					if (rpInfo->deltaList != 0) rpInfo->deltaList[index] = info->delta;
				}
			}
		} else if (type == AST_NAME_AVOGADRO) {
			rpInfo->varList[index] = 0;
			rpInfo->constList[index] = new double(6.0221367e+23);
			rpInfo->opfuncList[index] = 0;
			if (rpInfo->deltaList != 0) rpInfo->deltaList[index] = 0;
			freeConstList.push_back(rpInfo->constList[index]);
		} else if (type == AST_NAME_TIME) {
			variableInfo *info = searchInfoById(varInfoList, "t");
			if (info != 0) {
				rpInfo->varList[index] = 0;
				rpInfo->constList[index] = info->value;
				rpInfo->opfuncList[index] = 0;
				if (rpInfo->deltaList != 0) rpInfo->deltaList[index] = 0;
			}
		}
	}
	index++;
	if (index == index_max) index = 0;
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
  if (ast == NULL) return;
//  static int _count = 0;
//  ++_count;
//  std::string before = SBML_formulaToString(ast);
	unsigned int i, j;
	int type = (int)ast->getType();
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
	} else if (type == AST_TIMES) {//0 * x or x * 0 to 0
		if ((ast->getLeftChild()->isReal() && fabs(ast->getLeftChild()->getReal()) == 0) ||
			(ast->getLeftChild()->isInteger() && ast->getLeftChild()->getInteger() == 0) ||
			(ast->getRightChild()->isReal() && fabs(ast->getRightChild()->getReal()) == 0) ||
			(ast->getRightChild()->isInteger() && ast->getRightChild()->getInteger() == 0)) {
			ast->removeChild(0);
			ast->removeChild(0);
			ast->setType(AST_REAL);
			ast->setValue(0.0);
		}
	}

#if LIBSBML_VERSION >= 50903
  if (type == AST_TIMES || type == AST_PLUS)
  {
    ast->reduceToBinary();
  }
#endif 

	for (i = 0; i < ast->getNumChildren(); i++) {
		rearrangeAST(ast->getChild(i));
	}
  //std::string after = SBML_formulaToString(ast);
  //if (after != before)
  //{
  //  cout << "before: " << before << endl;
  //  cout << "after: " << after << endl << endl;
  //}
}

void parseDependence(const ASTNode *ast, vector<variableInfo*> &dependence, vector<variableInfo*> &varInfoList)
{
	unsigned int i;
	
  for (i = 0; i < ast->getNumChildren(); ++i) 
  {
		parseDependence(ast->getChild(i), dependence, varInfoList);
	}
	
  if (!ast->isName()) return;

  if (searchInfoById(dependence, ast->getName()) != 0) 
    return;

  variableInfo* temp = searchInfoById(varInfoList, ast->getName());
  if (temp == NULL)
    return;

  dependence.push_back(temp);
}
