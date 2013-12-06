#include <cmath>
#include "sbml/SBMLTypes.h"
#include "sbml/extension/SBMLExtensionRegistry.h"
#include "sbml/packages/req/common/RequiredElementsExtensionTypes.h"
#include "sbml/packages/spatial/common/SpatialExtensionTypes.h"
#include "sbml/packages/spatial/extension/SpatialModelPlugin.h"
#include "sbml/packages/spatial/extension/SpatialExtension.h"
#include <vector>
#include "mystruct.h"
#include "searchFunction.h"

#define stackMax 50

void reversePolishInitial(vector<int> &indexList, reversePolishInfo *rpInfo, double *value, int numOfASTNodes, int Xindex, int Yindex, int Zindex, bool isAllArea)
{
	int X, Y, Z, i, j, it_end = 0;
	int st_index = 0, index = 0;
	double rpStack[stackMax] = {0};
	if (!isAllArea) it_end = (int)indexList.size();
	else it_end = Xindex * Yindex * Zindex;
	for (j = 0; j < it_end; j++) {
		if (!isAllArea) index = indexList[j];
		else index = j;
		Z = index / (Xindex * Yindex);
		Y = (index - Z * Xindex * Yindex) / Xindex;
		X = index - Z * Xindex * Yindex - Y * Xindex;
		st_index = 0;
		for (i = 0; i < numOfASTNodes; i++) {
			if (rpInfo->varList[i] != 0) {//set variable into the stack
				rpStack[st_index] = rpInfo->varList[i][index];
				st_index++;
			} else if (rpInfo->constList[i] != 0) {//set const into the stack
				rpStack[st_index] = *(rpInfo->constList[i]);
				st_index++;
			} else {//operation
				st_index--;
				switch (rpInfo->opfuncList[i]) {
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
					rpStack[st_index] = acosh(rpStack[st_index]);
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
		value[index] = rpStack[--st_index];
	}
}

void reversePolishRK(reactionInfo *rInfo, GeometryInfo *geoInfo, int Xindex, int Yindex, int Zindex, double dt, int m, int numOfReactants, bool isReaction)
{
	int X, Y, Z, i, j;
	int k;
	int st_index = 0, index = 0, numOfVolIndexes = Xindex * Yindex * Zindex;
	double rpStack[stackMax] = {0};
	double rk[4] = {0, 0.5, 0.5, 1.0};
	double **variable = rInfo->rpInfo->varList;
	double **constant = rInfo->rpInfo->constList;
	double **d = rInfo->rpInfo->deltaList;
	int *operation = rInfo->rpInfo->opfuncList;
	int numOfASTNodes = rInfo->rpInfo->listNum;
	for (j = 0; j < (int)geoInfo->domainIndex.size(); j++) {
		index = geoInfo->domainIndex[j];
		Z = index / (Xindex * Yindex);
		Y = (index - Z * Xindex * Yindex) / Xindex;
		X = index - Z * Xindex * Yindex - Y * Xindex;
		st_index = 0;
		for (i = 0; i < numOfASTNodes; i++) {
			if (variable[i] != 0) {//set variable into the stack
				if (d != 0 && d[i] != 0) {
					if (m == 0) rpStack[st_index] = variable[i][index];
					else rpStack[st_index] = variable[i][index] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + index];
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
					rpStack[st_index] = acosh(rpStack[st_index]);
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
				//rp(rpStack, st_index, operation[i]);
			}
		}
		st_index--;
		if (isReaction) {//Reaction
			for (k = 0; k < numOfReactants; k++) {//reactants
				if (rInfo->isVariable[k]) rInfo->spRefList[k]->delta[m * numOfVolIndexes + index] -= rInfo->srStoichiometry[k] * rpStack[st_index];
			}
			for (k = numOfReactants; k < (int)rInfo->spRefList.size(); k++) {//products
				if (rInfo->isVariable[k]) rInfo->spRefList[k]->delta[m * numOfVolIndexes + index] += rInfo->srStoichiometry[k] * rpStack[st_index];
			}
		} else if (rInfo->isVariable[0]){//RateRule
			rInfo->spRefList[0]->delta[m * numOfVolIndexes + index] += rInfo->srStoichiometry[0] * rpStack[st_index];
		}
	}
}

void calcDiffusion(variableInfo *sInfo, double deltaX, double deltaY, double deltaZ, int Xindex, int Yindex, int Zindex, int m, double dt)
{
	int X = 0, Y = 0, Z = 0, index = 0;
	unsigned int j;
	int Xplus2 = 0, Xminus2 = 0, Yplus2 = 0, Yminus2 = 0, Zplus2 = 0, Zminus2 = 0;
	int numOfVolIndexes = Xindex * Yindex * Zindex;
	int dcIndex = 0;
	double* val = sInfo->value;
	double* d = sInfo->delta;
	double rk[4] = {0, 0.5, 0.5, 1.0};
	GeometryInfo *geoInfo = sInfo->geoi;
	//flux
	//2d
	//J = -D * dval / deltaX
	//x direction: d = (-J * deltaY) / (deltaY * deltaX) = -J / deltaX
	//3d
	//x direction: d = (-J * deltaY * deltaZ) / (deltaY * deltaZ * deltaX) = -J / deltaX

	// double Dx = deltaX, Dy = deltaY, Dz = deltaZ;
	for (j = 0; j < geoInfo->domainIndex.size(); j++) {
		index = geoInfo->domainIndex[j];
		Z = index / (Xindex * Yindex);
		Y = (index - Z * Xindex * Yindex) / Xindex;
		X = index - Z * Xindex * Yindex - Y * Xindex;
		Xplus2 = Z * Yindex * Xindex + Y * Xindex + (X + 2);
		Xminus2 = Z * Yindex * Xindex + Y * Xindex + (X - 2);
		Yplus2 = Z * Yindex * Xindex + (Y + 2) * Xindex + X;
		Yminus2 = Z * Yindex * Xindex + (Y - 2) * Xindex + X;
		Zplus2 = (Z + 2) * Yindex * Xindex + Y * Xindex + X;
		Zminus2 = (Z - 2) * Yindex * Xindex + Y * Xindex + X;
		if (sInfo->geoi->isDomain[index] == 1) {
			if (m == 0) {
				if (sInfo->diffCInfo[0] != 0) {//x-diffusion
					if (sInfo->diffCInfo[0]->isUniform == false) dcIndex = index;
					if (sInfo->geoi->bType[index].isBofXp == false) {
						sInfo->delta[m * numOfVolIndexes + index] += sInfo->diffCInfo[0]->value[dcIndex] * (val[Xplus2] - val[index]) / pow(deltaX, 2);
					}
					if (sInfo->geoi->bType[index].isBofXm == false) {
						sInfo->delta[m * numOfVolIndexes + index] += sInfo->diffCInfo[0]->value[dcIndex] * (val[Xminus2] - val[index]) / pow(deltaX, 2);
					}
				}
				if (sInfo->diffCInfo[1] != 0) {//y-diffusion
					if (sInfo->diffCInfo[0]->isUniform == false) dcIndex = index;
					if (sInfo->geoi->bType[index].isBofYp == false) {
						sInfo->delta[m * numOfVolIndexes + index] += sInfo->diffCInfo[1]->value[dcIndex] * (val[Yplus2] - val[index]) / pow(deltaY, 2);
					}
					if (sInfo->geoi->bType[index].isBofYm == false) {
						sInfo->delta[m * numOfVolIndexes + index] += sInfo->diffCInfo[1]->value[dcIndex] * (val[Yminus2] - val[index]) / pow(deltaY, 2);
					}
				}
				if (sInfo->diffCInfo[2] != 0) {//z-diffusion
					if (sInfo->diffCInfo[0]->isUniform == false) dcIndex = index;
					if (sInfo->geoi->bType[index].isBofZp == false) {
						sInfo->delta[m * numOfVolIndexes + index] += sInfo->diffCInfo[2]->value[dcIndex] * (val[Zplus2] - val[index]) / pow (deltaZ, 2);
					}
					if (sInfo->geoi->bType[index].isBofZm == false) {
						sInfo->delta[m * numOfVolIndexes + index] += sInfo->diffCInfo[2]->value[dcIndex] * (val[Zminus2] - val[index]) / pow (deltaZ, 2);
					}
				}
			} else {
				if (sInfo->diffCInfo[0] != 0) {//x-diffusion
					if (sInfo->diffCInfo[0]->isUniform == false) dcIndex = index;
					if (sInfo->geoi->bType[index].isBofXp == false) {
						sInfo->delta[m * numOfVolIndexes + index]
							+= sInfo->diffCInfo[0]->value[dcIndex] *
							((val[Xplus2] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + Xplus2])
							 - (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])) / pow(deltaX, 2);
					}
					if (sInfo->geoi->bType[index].isBofXm == false) {
						sInfo->delta[m * numOfVolIndexes + index]
							+= sInfo->diffCInfo[0]->value[dcIndex] *
							((val[Xminus2] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + Xminus2])
							 - (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])) / pow(deltaX, 2);
					}
				}
				if (sInfo->diffCInfo[1] != 0) {//y-diffusion
					if (sInfo->diffCInfo[0]->isUniform == false) dcIndex = index;
					if (sInfo->geoi->bType[index].isBofYp == false) {
						sInfo->delta[m * numOfVolIndexes + index]
							+= sInfo->diffCInfo[1]->value[dcIndex] *
							((val[Yplus2] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + Yplus2])
							 - (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])) / pow(deltaY, 2);
					}
					if (sInfo->geoi->bType[index].isBofYm == false) {
						sInfo->delta[m * numOfVolIndexes + index]
							+= sInfo->diffCInfo[1]->value[dcIndex] *
							((val[Yminus2] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + Yminus2])
							 - (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])) / pow(deltaY, 2);
					}
				}
				if (sInfo->diffCInfo[2] != 0) {//z-diffusion
					if (sInfo->diffCInfo[0]->isUniform == false) dcIndex = index;
					if (sInfo->geoi->bType[index].isBofZp == false) {
						sInfo->delta[m * numOfVolIndexes + index]
							+= sInfo->diffCInfo[2]->value[dcIndex] *
							((val[Zplus2] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + Zplus2])
							 - (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])) / pow(deltaZ, 2);
					}
					if (sInfo->geoi->bType[index].isBofZm == false) {
						sInfo->delta[m * numOfVolIndexes + index]
							+= sInfo->diffCInfo[2]->value[dcIndex] *
							((val[Zminus2] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + Zminus2])
							 - (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])) / pow(deltaZ, 2);
					}
				}
			}
		}
	}
}

void cipCSLR(variableInfo *sInfo, double deltaX, double deltaY, double deltaZ, double dt, int Xindex, int Yindex, int Zindex, int dimension)
{
	double ux = 0.0, uy = 0.0, uz = 0.0, Dx = 0.0, Dy = 0.0, Dz = 0.0;
	double a_i = 0.0, b_i = 0.0, beta_i = 0.0;//parameter of CIP rational function
	double a_i2 = 0.0, b_i2 = 0.0, beta_i2 = 0.0;
	double g_in = 0.0, g_out = 0.0, xi_x = 0.0, xi_y = 0.0, xi_z = 0.0;
	int X = 0, Y = 0, Z = 0, index = 0;
	int Xplus1 = 0, Xminus1 = 0, Yplus1 = 0, Yminus1 = 0, Zplus1 = 0, Zminus1 = 0;
	int Xplus2 = 0, Xminus2 = 0, Yplus2 = 0, Yminus2 = 0, Zplus2 = 0, Zminus2 = 0;
	int Xplus3 = 0, Xminus3 = 0, Yplus3 = 0, Yminus3 = 0, Zplus3 = 0, Zminus3 = 0;
	unsigned int i;
	double *val = sInfo->value;
	double *val_delta =  new double[Xindex * Yindex * Zindex];
	fill_n(val_delta, Xindex * Yindex * Zindex, 0);
	boundaryType type;
	for (int j = 0; j < Xindex * Yindex * Zindex; j++) val_delta[j] = 0.0;
	if (sInfo->adCInfo[0] != 0) {//x-direction
		for (i = 0; i < sInfo->geoi->domainIndex.size(); i++) {
			index = sInfo->geoi->domainIndex[i];
			g_in = 0.0;
			g_out = 0.0;
			Z = index / (Xindex * Yindex);
			Y = (index - Z * Xindex * Yindex) / Xindex;
			X = index - Z * Xindex * Yindex - Y * Xindex;
			Xplus1 = Z * Xindex * Yindex + Y * Xindex + (X + 1);
			Xplus2 = Z * Xindex * Yindex + Y * Xindex + (X + 2);
			Xplus3 = Z * Xindex * Yindex + Y * Xindex + (X + 3);
			Xminus1 = Z * Xindex * Yindex + Y * Xindex + (X - 1);
			Xminus2 = Z * Xindex * Yindex + Y * Xindex + (X - 2);
			Xminus3 = Z * Xindex * Yindex + Y * Xindex + (X - 3);
			type = sInfo->geoi->bType[index];
			//x-direction
			if (!(X + 4 < Xindex && sInfo->geoi->isDomain[Z * Xindex * Yindex + Y * Xindex + (X + 4)] == 1)) {
				if (X + 2 < Xindex && sInfo->geoi->isDomain[Xplus2] == 1) Xplus3 = Xplus2;
				else Xplus3 = index;
			}
			if (!(X - 4 >= 0 && sInfo->geoi->isDomain[Z * Xindex * Yindex + Y * Xindex + (X - 4)] == 1)) {
				if (X - 2 >= 0 && sInfo->geoi->isDomain[Xminus2] == 1) Xminus3 = Xminus2;
				else Xminus3 = index;
			}
			if (!(X + 2 < Xindex && sInfo->geoi->isDomain[Xplus2] == 1)) {
				Xplus1 = index;
				Xplus2 = index;
			}
			if (!(X - 2 >= 0 && sInfo->geoi->isDomain[Xminus2] == 1)) {
				Xminus1 = index;
				Xminus2 = index;
			}
			if (sInfo->adCInfo[0]->isUniform == true) ux = sInfo->adCInfo[0]->value[0];
			else ux = (sInfo->adCInfo[0]->value[index] + sInfo->adCInfo[0]->value[Xplus2]) / 2.0;
			Dx = -copysign(1.0, ux) * deltaX;
			xi_x = ux * dt;//Taylor
			if (ux >= 0 && !type.isBofXp) {//flux_out
				a_i = val[Xplus1];
				beta_i = ((fabs(val[Xplus1] - val[index]) + DBL_EPSILON) / (fabs(val[index] - val[Xminus1]) + DBL_EPSILON) - 1.0) / Dx;
				b_i = ((1.0 + beta_i * Dx) * val[index] - val[Xplus1]) / Dx;
				//next time value
				val_delta[Xplus1] = (a_i + 2 * b_i * (-xi_x) + beta_i * b_i * pow(-xi_x, 2)) / pow((1.0 + beta_i * (-xi_x)), 2) - val[Xplus1];
				//correction due to velocity divergence
				if (!sInfo->adCInfo[0]->isUniform) {
					val_delta[Xplus1] += (dt / (2 * deltaX)) * (val[Xplus1] + val_delta[Xplus1]) * (sInfo->adCInfo[0]->value[Xplus3] - sInfo->adCInfo[0]->value[Xminus1]);
				}
				g_out = (-a_i * xi_x + b_i * pow(xi_x, 2)) / (-1.0 + beta_i * xi_x);
			} else if (ux < 0 && !type.isBofXp) {//flux_in
				a_i = val[Xplus1];
				beta_i = ((fabs(val[Xplus1] - val[Xplus2]) + DBL_EPSILON) / (fabs(val[Xplus2] - val[Xplus3]) + DBL_EPSILON) - 1.0) / Dx;
				b_i = ((1.0 + beta_i * Dx) * val[Xplus2] - val[Xplus1]) / Dx;
				//next time value
				val_delta[Xplus1] = (a_i + 2 * b_i * (-xi_x) + beta_i * b_i * pow(-xi_x, 2)) / pow((1.0 + beta_i * (-xi_x)), 2) - val[Xplus1];
				//correction due to velocity divergence
				if (!sInfo->adCInfo[0]->isUniform) {
					val_delta[Xplus1] += (dt / (2 * deltaX)) * (val[Xplus1] + val_delta[Xplus1]) * (sInfo->adCInfo[0]->value[Xplus3] - sInfo->adCInfo[0]->value[Xminus1]);
				}
				g_in = -(-a_i * xi_x + b_i * pow(xi_x, 2)) / (-1.0 + beta_i * xi_x);
			}
			if (ux >= 0 && !type.isBofXm) {//flux_in
				a_i2 = val[Xminus1];
				beta_i2 = ((fabs(val[Xminus1] - val[Xminus2]) + DBL_EPSILON) / (fabs(val[Xminus2] - val[Xminus3]) + DBL_EPSILON) - 1.0) / Dx;
				b_i2 = ((1.0 + beta_i2 * Dx) * val[Xminus2] - val[Xminus1]) / Dx;
				g_in = (-a_i2 * xi_x + b_i2 * pow(xi_x, 2)) / (-1.0 + beta_i2 * xi_x);
			} else if (ux < 0 && !type.isBofXm) {//flux_out
				a_i2 = val[Xminus1];
				beta_i2 = ((fabs(val[Xminus1] - val[index]) + DBL_EPSILON) / (fabs(val[index] - val[Xplus1]) + DBL_EPSILON) - 1.0) / Dx;
				b_i2 = ((1.0 + beta_i2 * Dx) * val[index] - val[Xminus1]) / Dx;
				g_out = -(-a_i2 * xi_x + b_i2 * pow(xi_x, 2)) / (-1.0 + beta_i2 * xi_x);
			}
			val_delta[index] = (g_in - g_out) / deltaX;
		}
	}
	//update
	for (Z = 0; Z < Zindex; Z++) {
		for (Y = 0; Y < Yindex; Y++) {
			for (X = 0; X < Xindex; X++) {
				index = Z * Xindex * Yindex + Y * Xindex + X;
				val[index] += val_delta[index];
				if (X % 2 == 0 && Y % 2 == 0 && Z % 2 == 0 && sInfo->geoi->isDomain[index] == 1) {
					if (dimension >= 2 && !sInfo->geoi->bType[index].isBofYp) {//update val(x, y+1/2, z)
						Yplus1 = Z * Xindex * Yindex + (Y + 1) * Xindex + X;
						Yplus2 = Z * Xindex * Yindex + (Y + 2) * Xindex + X;
						val[Yplus1] += (val_delta[Yplus2] + val_delta[index]) / 2.0;
					}
					if (dimension == 3 && !sInfo->geoi->bType[index].isBofZp) {//update val(x, y, z+1/2)
						Zplus1 = (Z + 1) * Xindex * Yindex + Y * Xindex + X;
						Zplus2 = (Z + 2) * Xindex * Yindex + Y * Xindex + X;
						val[Zplus1] += (val_delta[Zplus2] + val_delta[index]) / 2.0;
					}
				}
			}
		}
	}
	for (int j = 0; j < Xindex * Yindex * Zindex; j++) val_delta[j] = 0.0;
	if (dimension >= 2) {
		if (sInfo->adCInfo[1] != 0) {//y-direction
			for (i = 0; i < sInfo->geoi->domainIndex.size(); i++) {
				index = sInfo->geoi->domainIndex[i];
				g_in = 0.0;
				g_out = 0.0;
				Z = index / (Xindex * Yindex);
				Y = (index - Z * Xindex * Yindex) / Xindex;
				X = index - Z * Xindex * Yindex - Y * Xindex;
				Yplus1 = Z * Xindex * Yindex + (Y + 1) * Xindex + X;
				Yplus2 = Z * Xindex * Yindex + (Y + 2) * Xindex + X;
				Yplus3 = Z * Xindex * Yindex + (Y + 3) * Xindex + X;
				Yminus1 = Z * Xindex * Yindex + (Y - 1) * Xindex + X;
				Yminus2 = Z * Xindex * Yindex + (Y - 2) * Xindex + X;
				Yminus3 = Z * Xindex * Yindex + (Y - 3) * Xindex + X;
				type = sInfo->geoi->bType[index];
				if (!(Y + 4 < Yindex && sInfo->geoi->isDomain[Z * Xindex * Yindex + (Y + 4) * Xindex + X] == 1)) {
					if (Y + 2 < Yindex && sInfo->geoi->isDomain[Yplus2] == 1) Yplus3 = Yplus2;
					else Yplus3 = index;
				}
				if (!(Y - 4 >= 0 && sInfo->geoi->isDomain[Z * Xindex * Yindex + (Y - 4) * Xindex + X] == 1)) {
					if (Y - 2 >= 0 && sInfo->geoi->isDomain[Yminus2] == 1) Yminus3 = Yminus2;
					else Yminus3 = index;
				}

				if (!(Y + 2 < Yindex && sInfo->geoi->isDomain[Yplus2] == 1)) {
					Yplus1 = index;
					Yplus2 = index;
				}
				if (!(Y - 2 >= 0 && sInfo->geoi->isDomain[Yminus2] == 1)) {
					Yminus1 = index;
					Yminus2 = index;
				}
				if (sInfo->adCInfo[1]->isUniform == true) uy = sInfo->adCInfo[1]->value[0];
				else uy = (sInfo->adCInfo[1]->value[index] + sInfo->adCInfo[1]->value[Yplus2]) / 2.0;
				Dy = -copysign(1.0, uy) * deltaY;
				xi_y = uy * dt;//Taylor
				if (uy >= 0 && !type.isBofYp) {//flux_out
					a_i = val[Yplus1];
					beta_i = ((fabs(val[Yplus1] - val[index]) + DBL_EPSILON) / (fabs(val[index] - val[Yminus1]) + DBL_EPSILON) - 1.0) / Dy;
					b_i = ((1.0 + beta_i * Dy) * val[index] - val[Yplus1]) / Dy;
					//next time value
					val_delta[Yplus1] = (a_i + 2 * b_i * (-xi_y) + beta_i * b_i * pow(-xi_y, 2)) / pow((1.0 + beta_i * (-xi_y)), 2) - val[Yplus1];
					//correction due to velocity divergence
					if (!sInfo->adCInfo[1]->isUniform) {
						val_delta[Yplus1] += dt * (val[Yplus1] + val_delta[Yplus1]) * (sInfo->adCInfo[1]->value[Yplus3] - sInfo->adCInfo[1]->value[Yminus1]) / (2 * deltaY);
					}
					g_out = (-a_i * xi_y + b_i * pow(xi_y, 2)) / (-1.0 + beta_i * xi_y);
				} else if (uy < 0 && !type.isBofYp) {//flux_in
					a_i = val[Yplus1];
					beta_i = ((fabs(val[Yplus1] - val[Yplus2]) + DBL_EPSILON) / (fabs(val[Yplus2] - val[Yplus3]) + DBL_EPSILON) - 1.0) / Dy;
					b_i = ((1.0 + beta_i * Dy) * val[Yplus2] - val[Yplus1]) / Dy;
					//next time value
					val_delta[Yplus1] = (a_i + 2 * b_i * (-xi_y) + beta_i * b_i * pow(-xi_y, 2)) / pow((1.0 + beta_i * (-xi_y)), 2) - val[Yplus1];
					//correction due to velocity divergence
					if (!sInfo->adCInfo[1]->isUniform) {
						val_delta[Yplus1] += dt * (val[Yplus1] + val_delta[Yplus1]) * (sInfo->adCInfo[1]->value[Yplus3] - sInfo->adCInfo[1]->value[Yminus1]) / (2 * deltaY);
					}
					g_in = -(-a_i * xi_y + b_i * pow(xi_y, 2)) / (-1.0 + beta_i * xi_y);
				}
				if (uy >= 0 && !type.isBofYm) {//flux_in
					a_i2 = val[Yminus1];
					beta_i2 = ((fabs(val[Yminus1] - val[Yminus2]) + DBL_EPSILON) / (fabs(val[Yminus2] - val[Yminus3]) + DBL_EPSILON) - 1.0) / Dy;
					b_i2 = ((1.0 + beta_i2 * Dy) * val[Yminus2] - val[Yminus1]) / Dy;
					g_in = (-a_i2 * xi_y + b_i2 * pow(xi_y, 2)) / (-1.0 + beta_i2 * xi_y);
				} else if (uy < 0 && !type.isBofYm) {//flux_out
					a_i2 = val[Yminus1];
					beta_i2 = ((fabs(val[Yminus1] - val[index]) + DBL_EPSILON) / (fabs(val[index] - val[Yplus1]) + DBL_EPSILON) - 1.0) / Dy;
					b_i2 = ((1.0 + beta_i2 * Dy) * val[index] - val[Yminus1]) / Dy;
					g_out = -(-a_i2 * xi_y + b_i2 * pow(xi_y, 2)) / (-1.0 + beta_i2 * xi_y);
				}
				val_delta[index] = (g_in - g_out) / deltaY;
			}
		}
		//update
		for (Z = 0; Z < Zindex; Z++) {
			for (Y = 0; Y < Yindex; Y++) {
				for (X = 0; X < Xindex; X++) {
					index = Z * Xindex * Yindex + Y * Xindex + X;
					val[index] += val_delta[index];
					if (X % 2 == 0 && Y % 2 == 0 && Z % 2 == 0 && sInfo->geoi->isDomain[index] == 1) {
						if (!sInfo->geoi->bType[index].isBofXp) {//update val(x+1/2, y ,z)
							Xplus1 = Z * Xindex * Yindex + Y * Xindex + (X + 1);
							Xplus2 = Z * Xindex * Yindex + Y * Xindex + (X + 2);
							val[Xplus1] += (val_delta[Xplus2] + val_delta[index]) / 2.0;
						}
						if (dimension == 3 && !sInfo->geoi->bType[index].isBofZp) {//update val(x, y, z+1/2)
							Zplus1 = (Z + 1) * Xindex * Yindex + Y * Xindex + X;
							Zplus2 = (Z + 2) * Xindex * Yindex + Y * Xindex + X;
							val[Zplus1] += (val_delta[Zplus2] + val_delta[index]) / 2.0;
						}
					}
				}
			}
		}
		for (int j = 0; j < Xindex * Yindex * Zindex; j++) val_delta[j] = 0.0;
		if (dimension == 3) {
			if (sInfo->adCInfo[2] != 0) {//z-direction
				for (i = 0; i < sInfo->geoi->domainIndex.size(); i++) {
					index = sInfo->geoi->domainIndex[i];
					g_in = 0.0;
					g_out = 0.0;
					Z = index / (Xindex * Yindex);
					Y = (index - Z * Xindex * Yindex) / Xindex;
					X = index - Z * Xindex * Yindex - Y * Xindex;
					Zplus1 = (Z + 1) * Xindex * Yindex + Y * Xindex + X;
					Zplus2 = (Z + 2) * Xindex * Yindex + Y * Xindex + X;
					Zplus3 = (Z + 3)* Xindex * Yindex + Y * Xindex + X;
					Zminus1 = (Z - 1) * Xindex * Yindex + Y * Xindex + X;
					Zminus2 = (Z - 2) * Xindex * Yindex + Y * Xindex + X;
					Zminus3 = (Z - 3) * Xindex * Yindex + Y * Xindex + X;
					type = sInfo->geoi->bType[index];
					if (!(Z + 4 < Zindex && sInfo->geoi->isDomain[(Z + 4) * Xindex * Yindex + Y * Xindex + X] == 1)) {
						if (Z + 2 < Zindex && sInfo->geoi->isDomain[Zplus2] == 1) Zplus3 = Zplus2;
						else Zplus3 = index;
					}
					if (!(Z - 4 >= 0 && sInfo->geoi->isDomain[(Z - 4) * Xindex * Yindex + Y * Xindex + X] == 1)) {
						if (Z - 2 >= 0 && sInfo->geoi->isDomain[Zminus2] == 1) Zminus3 = Zminus2;
						else Zminus3 = index;
					}

					if (!(Z + 2 < Zindex && sInfo->geoi->isDomain[Zplus2] == 1)) {
						Zplus1 = index;
						Zplus2 = index;
					}
					if (!(Z - 2 >= 0 && sInfo->geoi->isDomain[Zminus2] == 1)) {
						Zminus1 = index;
						Zminus2 = index;
					}
					if (sInfo->adCInfo[2]->isUniform == true) uz = sInfo->adCInfo[2]->value[0];
					else uz = (sInfo->adCInfo[2]->value[index] + sInfo->adCInfo[2]->value[Zplus2]) / 2.0;
					Dz = -copysign(1.0, uz) * deltaZ;
					xi_z = uz * dt;
					if (uz >= 0 && !type.isBofZp) {//flux_out
						a_i = val[Zplus1];
						beta_i = ((fabs(val[Zplus1] - val[index]) + DBL_EPSILON) / (fabs(val[index] - val[Zminus1]) + DBL_EPSILON) - 1.0) / Dz;
						b_i = ((1.0 + beta_i * Dz) * val[index] - val[Zplus1]) / Dz;
						//next time value
						val_delta[Zplus1] = (a_i + 2 * b_i * (-xi_z) + beta_i * b_i * pow(-xi_z, 2)) / pow((1.0 + beta_i * (-xi_z)), 2) - val[Zplus1];
						//correction due to velocity divergence
						if (!sInfo->adCInfo[2]->isUniform) {
							val_delta[Zplus1] += dt * (val[Zplus1] + val_delta[Zplus1]) * (sInfo->adCInfo[2]->value[Zplus3] - sInfo->adCInfo[2]->value[Zminus1]) / (2 * deltaZ);
						}
						g_out = (-a_i * xi_z + b_i * pow(xi_z, 2)) / (-1.0 + beta_i * xi_z);
					} else if (uz < 0 && !type.isBofZp) {//flux_in
						a_i = val[Zplus1];
						beta_i = ((fabs(val[Zplus1] - val[Zplus2]) + DBL_EPSILON) / (fabs(val[Zplus2] - val[Zplus3]) + DBL_EPSILON) - 1.0) / Dz;
						b_i = ((1.0 + beta_i * Dz) * val[Zplus2] - val[Zplus1]) / Dz;
						//next time value
						val_delta[Zplus1] = (a_i + 2 * b_i * (-xi_z) + beta_i * b_i * pow(-xi_z, 2)) / pow((1.0 + beta_i * (-xi_z)), 2) - val[Zplus1];
						//correction due to velocity divergence
						if (!sInfo->adCInfo[2]->isUniform) {
							val_delta[Zplus1] += dt * (val[Zplus1] + val_delta[Zplus1]) * (sInfo->adCInfo[2]->value[Zplus3] - sInfo->adCInfo[2]->value[Zminus1]) / (2 * deltaZ);
						}
						g_in = -(-a_i * xi_z + b_i * pow(xi_z, 2)) / (-1.0 + beta_i * xi_z);
					}
					if (uz >= 0 && !type.isBofZm) {//flux_in
						a_i2 = val[Zminus1];
						beta_i2 = ((fabs(val[Zminus1] - val[Zminus2]) + DBL_EPSILON) / (fabs(val[Zminus2] - val[Zminus3]) + DBL_EPSILON) - 1.0) / Dz;
						b_i2 = ((1.0 + beta_i2 * Dz) * val[Zminus2] - val[Zminus1]) / Dz;
						g_in = (-a_i2 * xi_z + b_i2 * pow(xi_z, 2)) / (-1.0 + beta_i2 * xi_z);
					} else if (uz < 0 && !type.isBofZm) {//flux_out
						a_i2 = val[Zminus1];
						beta_i2 = ((fabs(val[Zminus1] - val[index]) + DBL_EPSILON) / (fabs(val[index] - val[Zplus1]) + DBL_EPSILON) - 1.0) / Dz;
						b_i2 = ((1.0 + beta_i2 * Dz) * val[index] - val[Zminus1]) / Dz;
						g_out = -(-a_i2 * xi_z + b_i2 * pow(xi_z, 2)) / (-1.0 + beta_i2 * xi_z);
					}
					val_delta[index] = (g_in - g_out) / deltaZ;
				}
			}
			//update
			for (Z = 0; Z < Zindex; Z++) {
				for (Y = 0; Y < Yindex; Y++) {
					for (X = 0; X < Xindex; X++) {
						index = Z * Xindex * Yindex + Y * Xindex + X;
						val[index] += val_delta[index];
						if (X % 2 == 0 && Y % 2 == 0 && Z % 2 == 0 && sInfo->geoi->isDomain[index] == 1) {
							if (!sInfo->geoi->bType[index].isBofXp) {//update val(x+1/2, y ,z)
								Xplus1 = Z * Xindex * Yindex + Y * Xindex + (X + 1);
								Xplus2 = Z * Xindex * Yindex + Y * Xindex + (X + 2);
								val[Xplus1] += (val_delta[Xplus2] + val_delta[index]) / 2.0;
							}
							if (dimension >= 2 && !sInfo->geoi->bType[index].isBofYp) {//update val(x, y+1/2, z)
								Yplus1 = Z * Xindex * Yindex + (Y + 1) * Xindex + X;
								Yplus2 = Z * Xindex * Yindex + (Y + 2) * Xindex + X;
								val[Yplus1] += (val_delta[Yplus2] + val_delta[index]) / 2.0;
							}
						}
					}
				}
			}
		}
	}
	delete[] val_delta;
}

void calcBoundary(variableInfo *sInfo, double deltaX, double deltaY, double deltaZ, int Xindex, int Yindex, int Zindex, int m, int dimension)
{
	int Xp = 0, Xm = 0, Yp = 0, Ym = 0, Zp = 0, Zm = 0, X = 0, Y = 0, Z = 0;
	int divIndexXp = 0, divIndexXm = 0, divIndexYp = 0,divIndexYm = 0, divIndexZp = 0, divIndexZm = 0;
	int numOfVolIndexes = Xindex * Yindex * Zindex;
	// 	int Xdiv = (Xindex + 1) / 2, Ydiv = (Yindex + 1) / 2, Zdiv = (Zindex + 1) / 2;
	BoundaryCondition *maxSideBC = 0, *minSideBC = 0;
	//boundary flux
	//2d
	//x direction: d = (-J * deltaY) / (deltaY * (deltaX / 2.0)) = -2.0 * J / deltaX
	//3d
	//x direction: d = (-J * deltaY * deltaZ) / (deltaY * deltaZ * (deltaX / 2.0)) = -2.0 * J / deltaX
	if (dimension >= 1) {
		maxSideBC = static_cast<SpatialParameterPlugin*>(sInfo->boundaryInfo[Xmax]->para->getPlugin("spatial"))->getBoundaryCondition();
		minSideBC = static_cast<SpatialParameterPlugin*>(sInfo->boundaryInfo[Xmin]->para->getPlugin("spatial"))->getBoundaryCondition();
		//Xp, Xm
		for (Z = 0; Z < Zindex; Z += 2) {
			for (Y = 0; Y < Yindex; Y += 2) {
				Xp = Z * Yindex * Xindex + Y * Xindex + (Xindex - 1);
				Xm = Z * Yindex * Xindex + Y * Xindex;
				if (!sInfo->boundaryInfo[Xmax]->isUniform) divIndexXp = Xp;
				if (!sInfo->boundaryInfo[Xmin]->isUniform) divIndexXm = Xm;
				if (sInfo->geoi->isDomain[Xp] == 1) {//Xp
					if (maxSideBC->getType() == "Flux") sInfo->delta[m * numOfVolIndexes + Xp] += 2.0 * (-sInfo->boundaryInfo[Xmax]->value[divIndexXp]) / deltaX;
					else if (maxSideBC->getType() == "Value") sInfo->value[Xp] = sInfo->boundaryInfo[Xmax]->value[divIndexXp];
				}
				if (sInfo->geoi->isDomain[Xm] == 1) {//Xm
					if (minSideBC->getType() == "Flux") sInfo->delta[m * numOfVolIndexes + Xm] += -2.0 * (-sInfo->boundaryInfo[Xmin]->value[divIndexXm]) / deltaX;
					else if (minSideBC->getType() == "Value") sInfo->value[Xm] = sInfo->boundaryInfo[Xmin]->value[divIndexXm];
				}
			}
		}
	}
	//Yp, Ym
	if (dimension >= 2) {
		maxSideBC = static_cast<SpatialParameterPlugin*>(sInfo->boundaryInfo[Ymax]->para->getPlugin("spatial"))->getBoundaryCondition();
		minSideBC = static_cast<SpatialParameterPlugin*>(sInfo->boundaryInfo[Ymin]->para->getPlugin("spatial"))->getBoundaryCondition();
		for (Z = 0; Z < Zindex; Z += 2) {
			for (X = 0; X < Xindex; X += 2) {
				Yp = Z * Yindex * Xindex + (Yindex - 1) * Xindex + X;
				Ym = Z * Yindex * Xindex + X;
				if (!sInfo->boundaryInfo[Ymax]->isUniform) divIndexYp = Yp;
				if (!sInfo->boundaryInfo[Ymin]->isUniform) divIndexYm = Ym;
				if (sInfo->geoi->isDomain[Yp] == 1) {//Yp
					if (maxSideBC->getType() == "Flux")	sInfo->delta[m * numOfVolIndexes + Yp] += 2.0 * (-sInfo->boundaryInfo[Ymax]->value[divIndexYp]) / deltaY;
					else if (maxSideBC->getType() == "Value") sInfo->value[Yp] = sInfo->boundaryInfo[Ymax]->value[divIndexYp];
				}
				if (sInfo->geoi->isDomain[Ym] == 1) {//Ym
					if (minSideBC->getType() == "Flux") sInfo->delta[m * numOfVolIndexes + Ym] += -2.0 * (-sInfo->boundaryInfo[Ymin]->value[divIndexYm]) / deltaY;
					else if (minSideBC->getType() == "Value") sInfo->value[Ym] = sInfo->boundaryInfo[Ymin]->value[divIndexYm];
				}
			}
		}
	}
	//Zp, Zm
	if (dimension >= 3) {
		maxSideBC = static_cast<SpatialParameterPlugin*>(sInfo->boundaryInfo[Zmax]->para->getPlugin("spatial"))->getBoundaryCondition();
		minSideBC = static_cast<SpatialParameterPlugin*>(sInfo->boundaryInfo[Zmin]->para->getPlugin("spatial"))->getBoundaryCondition();
		for (Y = 0; Y < Yindex; Y += 2) {
			for (X = 0; X < Xindex; X += 2) {
				Zp = (Zindex - 1) * Yindex * Xindex + Y * Xindex + X;
				Zm = Y * Xindex + X;
				if (!sInfo->boundaryInfo[Zmax]->isUniform) divIndexZp = Zp;
				if (!sInfo->boundaryInfo[Zmin]->isUniform) divIndexZm = Zm;
				if (sInfo->geoi->isDomain[Zp] == 1) {//Zp
					if (maxSideBC->getType() == "Flux") sInfo->delta[m * numOfVolIndexes + Zp] += 2.0 * (-sInfo->boundaryInfo[Zmax]->value[divIndexZp]) / deltaZ;
					else if (maxSideBC->getType() == "Value") sInfo->value[Zp] = sInfo->boundaryInfo[Zmax]->value[divIndexZp];
				}
				if (sInfo->geoi->isDomain[Zm] == 1) {//Zm
					if (minSideBC->getType() == "Flux") sInfo->delta[m * numOfVolIndexes + Zm] += -2.0 * (-sInfo->boundaryInfo[Zmin]->value[divIndexZm]) / deltaZ;
					else if (minSideBC->getType() == "Value") sInfo->value[Zm] = sInfo->boundaryInfo[Zmin]->value[divIndexZm];
				}
			}
		}
	}
}

void calcMemTransport(reactionInfo *rInfo, GeometryInfo *geoInfo, normalUnitVector *nuVec, int Xindex, int Yindex, int Zindex, double dt, int m, double deltaX, double deltaY, double deltaZ, int dimension, int numOfReactants)
{
	int X, Y, Z, i, j, k;
	int st_index = 0, index = 0, numOfVolIndexes = Xindex * Yindex * Zindex;
	int Xplus1 = 0, Xminus1 = 0, Yplus1 = 0, Yminus1 = 0, Zplus1 = 0, Zminus1 = 0;
	int Xplus3 = 0, Xminus3 = 0, Yplus3 = 0, Yminus3 = 0, Zplus3 = 0, Zminus3 = 0;
	double rpStack[stackMax] = {0};
	double rk[4] = {0, 0.5, 0.5, 1.0};
	double **variable = rInfo->rpInfo->varList;
	double **constant = rInfo->rpInfo->constList;
	double **d = rInfo->rpInfo->deltaList;
	int *operation = rInfo->rpInfo->opfuncList;
	int numOfASTNodes = rInfo->rpInfo->listNum;
	variableInfo *symbolInfo = 0;
	for (k = 0; k < (int)geoInfo->domainIndex.size(); k++) {
		index = geoInfo->domainIndex[k];
		Z = index / (Xindex * Yindex);
		Y = (index - Z * Xindex * Yindex) / Xindex;
		X = index - Z * Xindex * Yindex - Y * Xindex;
		if ((dimension == 3 && (X * Y * Z == 0 || X == Xindex - 1 || Y == Yindex - 1 || Z == Zindex - 1))
			|| (dimension == 2 && (X * Y == 0 || X == Xindex - 1 || Y == Yindex - 1))) {
			continue;
		}
		st_index = 0;
		Xplus3 = Z * Yindex * Xindex + Y * Xindex + (X + 3);
		Xminus3 = Z * Yindex * Xindex + Y * Xindex + (X - 3);
		Yplus3 = Z * Yindex * Xindex + (Y + 3) * Xindex + X;
		Yminus3 = Z * Yindex * Xindex + (Y - 3) * Xindex + X;
		Zplus3 = (Z + 3) * Yindex * Xindex + Y * Xindex + X;
		Zminus3 = (Z - 3) * Yindex * Xindex + Y * Xindex + X;
		Xplus1 = Z * Yindex * Xindex + Y * Xindex + (X + 1);
		Xminus1 = Z * Yindex * Xindex + Y * Xindex + (X - 1);
		Yplus1 = Z * Yindex * Xindex + (Y + 1) * Xindex + X;
		Yminus1 = Z * Yindex * Xindex + (Y - 1) * Xindex + X;
		Zplus1 = (Z + 1) * Yindex * Xindex + Y * Xindex + X;
		Zminus1 = (Z - 1) * Yindex * Xindex + Y * Xindex + X;

		if (geoInfo->isDomain[index] == 1) {//not pseudo membrane
			for (i = 0; i < numOfASTNodes; i++) {
				if (variable[i] != 0) {//set variable into the stack
					if (d != 0 && d[i] != 0) {
						for (j = 0; j < static_cast<int>(rInfo->spRefList.size()); j++) {
							//search compartment of equation's symbol
							if (variable[i] == rInfo->spRefList[j]->value) {
								symbolInfo = rInfo->spRefList[j];
								break;
							}
						}
						/*
						  a volume symbol's value at membrane is calculated with linear approximation
						  value = value(boundary) + (value(boundary) - value(boundary_next)) / 2
						  = 1.5 * value(boundary) - 0.5 * value(boundary_next)
						*/
						if (symbolInfo->geoi->isVol) {//symbol is in volume
							//x transport
							if (static_cast<int>(symbolInfo->geoi->isDomain[Xplus1]) == 1) {//right of membrane
								if (Xplus3 < numOfVolIndexes && symbolInfo->geoi->isDomain[Xplus3] == 1) {
									if (m == 0) rpStack[st_index] = 1.5 * variable[i][Xplus1] - 0.5 * variable[i][Xplus3];
									else {
										rpStack[st_index] = 1.5 * (variable[i][Xplus1] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Xplus1])
											- 0.5 * (variable[i][Xplus3] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Xplus3]);
									}
								} else {
									if (m == 0) rpStack[st_index] = variable[i][Xplus1];
									else rpStack[st_index] = variable[i][Xplus1] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Xplus1];
								}
							} else if (symbolInfo->geoi->isDomain[Xminus1] == 1) {//left of membrane
								if (Xminus3 < numOfVolIndexes && symbolInfo->geoi->isDomain[Xminus3] == 1) {
									if (m == 0) rpStack[st_index] = 1.5 * variable[i][Xminus1] - 0.5 * variable[i][Xminus3];
									else {
										rpStack[st_index] = 1.5 * (variable[i][Xminus1] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Xminus1])
											- 0.5 * (variable[i][Xminus3] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Xminus3]);
									}
								} else {
									if (m == 0) rpStack[st_index] = variable[i][Xminus1];
									else rpStack[st_index] = variable[i][Xminus1] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Xminus1];
								}
							}
							//y transport
							if (static_cast<int>(symbolInfo->geoi->isDomain[Yplus1]) == 1) {//upper of membrane
								if (Yplus3 < numOfVolIndexes && symbolInfo->geoi->isDomain[Yplus3] == 1) {
									if (m == 0) rpStack[st_index] = 1.5 * variable[i][Yplus1] - 0.5 * variable[i][Yplus3];
									else {
										rpStack[st_index] = 1.5 * (variable[i][Yplus1] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Yplus1])
											- 0.5 * (variable[i][Yplus3] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Yplus3]);
									}
								} else {
									if (m == 0) rpStack[st_index] = variable[i][Yplus1];
									else rpStack[st_index] = variable[i][Yplus1] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Yplus1];
								}
							} else if (static_cast<int>(symbolInfo->geoi->isDomain[Yminus1]) == 1) {//downer of membrane
								if (Yminus3 < numOfVolIndexes && symbolInfo->geoi->isDomain[Yminus3] == 1) {
									if (m == 0) {
										rpStack[st_index] = 1.5 * variable[i][Yminus1] - 0.5 * variable[i][Yminus3];
									} else {
										rpStack[st_index] = 1.5 * (variable[i][Yminus1] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Yminus1])
											- 0.5 * (variable[i][Yminus3] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Yminus3]);
									}
								} else {
									if (m == 0) rpStack[st_index] = variable[i][Yminus1];
									else rpStack[st_index] = variable[i][Yminus1] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Yminus1];
								}
							}
							//z transport
							if (dimension == 3) {
								if (static_cast<int>(symbolInfo->geoi->isDomain[Zplus1]) == 1) {//higher of membrane
									if (Zplus3 < numOfVolIndexes && symbolInfo->geoi->isDomain[Zplus3] == 1) {
										if (m == 0) rpStack[st_index] = 1.5 * variable[i][Zplus1] - 0.5 * variable[i][Zplus3];
										else {
											rpStack[st_index] = 1.5 * (variable[i][Zplus1] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Zplus1])
												- 0.5 * (variable[i][Zplus3] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Zplus3]);
										}
									} else {
										if (m == 0) rpStack[st_index] = variable[i][Zplus1];
										else rpStack[st_index] = variable[i][Zplus1] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Zplus1];
									}
								} else if (static_cast<int>(symbolInfo->geoi->isDomain[Zminus1]) == 1) {//lowner of membrane
									if (Zminus3 < numOfVolIndexes && symbolInfo->geoi->isDomain[Zminus3] == 1) {
										if (m == 0) rpStack[st_index] = 1.5 * variable[i][Zminus1] - 0.5 * variable[i][Zminus3];
										else {
											rpStack[st_index] = 1.5 * (variable[i][Zminus1] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Zminus1])
												- 0.5 * (variable[i][Zminus3] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Zminus3]);
										}
									} else {
										if (m == 0) rpStack[st_index] = variable[i][Zminus1];
										else rpStack[st_index] = variable[i][Zminus1] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + Zminus1];
									}
								}
							}
						} else {//symbol is in membrane
							if (m == 0) rpStack[st_index] = variable[i][index];
							else rpStack[st_index] = variable[i][index] + rk[m] * dt * d[i][(m - 1) * numOfVolIndexes + index];
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
						rpStack[st_index] = acosh(rpStack[st_index]);
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
			//area = dx * dy
			//du / dt = -(Jx * dy) / area = -Jx / dx, du / dt = -Jy / dy
			st_index--;
			for (j = 0; j < numOfReactants; j++) {//reactants
				if (rInfo->isVariable[j]) {
					if (geoInfo->bType[index].isBofXp || geoInfo->bType[index].isBofXm) {//x transport or x binding
						//transport
						if (rInfo->spRefList[j]->geoi->isDomain[Xplus1] == 1) {//right of membrane
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + Xplus1] -= fabs(nuVec[index].nx) * rInfo->srStoichiometry[j] * rpStack[st_index] / deltaX;
						}
						if (rInfo->spRefList[j]->geoi->isDomain[Xminus1] == 1) {//left of membrane
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + Xminus1] -= fabs(nuVec[index].nx) * rInfo->srStoichiometry[j] * rpStack[st_index] / deltaX;
						}
						//binding
						if (rInfo->spRefList[j]->geoi->isDomain[index] == 1) {
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + index] -= fabs(nuVec[index].nx) * rInfo->srStoichiometry[j] * rpStack[st_index] / (deltaX / 2.0);
						}
					} else if (geoInfo->bType[index].isBofYp || geoInfo->bType[index].isBofYm) {//y transport or y binding
						//transport
						if (rInfo->spRefList[j]->geoi->isDomain[Yplus1] == 1) {//upper of membrane
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + Yplus1] -= fabs(nuVec[index].ny) * rInfo->srStoichiometry[j] * rpStack[st_index] / deltaY;
						}
						if (rInfo->spRefList[j]->geoi->isDomain[Yminus1] == 1) {//downer of membrane
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + Yminus1] -= fabs(nuVec[index].ny) * rInfo->srStoichiometry[j] * rpStack[st_index] / deltaY;
						}
						//binding
						if (rInfo->spRefList[j]->geoi->isDomain[index] == 1) {
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + index] -= fabs(nuVec[index].ny) * rInfo->srStoichiometry[j] * rpStack[st_index] / (deltaY / 2.0);
						}
					} else if (dimension == 3 && (geoInfo->bType[index].isBofZp || geoInfo->bType[index].isBofZm)) {//z transport
						//transport
						if (rInfo->spRefList[j]->geoi->isDomain[Zplus1] == 1) {//higher of membrane
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + Zplus1] -= fabs(nuVec[index].nz) * rInfo->srStoichiometry[j] * rpStack[st_index] / deltaZ;
						}
						if (rInfo->spRefList[j]->geoi->isDomain[Zminus1] == 1) {//lower of membrane
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + Zminus1] -= fabs(nuVec[index].nz) * rInfo->srStoichiometry[j] * rpStack[st_index] / deltaZ;
						}
						//binding
						if (rInfo->spRefList[j]->geoi->isDomain[index] == 1) {
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + index] -= fabs(nuVec[index].nz) * rInfo->srStoichiometry[j] * rpStack[st_index] / (deltaZ / 2.0);
						}
					}
				}
			}
			for (j = numOfReactants; j < (int)rInfo->spRefList.size(); j++) {//products
				if (rInfo->isVariable[j]) {
					if (geoInfo->bType[index].isBofXp || geoInfo->bType[index].isBofXm) {//x transport or x binding
						//transport
						if (rInfo->spRefList[j]->geoi->isDomain[Xplus1] == 1) {//right of membrane
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + Xplus1] += fabs(nuVec[index].nx) * rInfo->srStoichiometry[j] * rpStack[st_index] / deltaX;
						}
						if (rInfo->spRefList[j]->geoi->isDomain[Xminus1] == 1) {//left of membrane
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + Xminus1] += fabs(nuVec[index].nx) * rInfo->srStoichiometry[j] * rpStack[st_index] / deltaX;
						}
						//binding
						if (rInfo->spRefList[j]->geoi->isDomain[index] == 1) {
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + index] += fabs(nuVec[index].nx) * rInfo->srStoichiometry[j] * rpStack[st_index] / (deltaX / 2.0);
						}
					} else if (geoInfo->bType[index].isBofYp || geoInfo->bType[index].isBofYm) {//y transport or y binding
						//transport
						if (rInfo->spRefList[j]->geoi->isDomain[Yplus1] == 1) {//upper of membrane
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + Yplus1] += fabs(nuVec[index].ny) * rInfo->srStoichiometry[j] * rpStack[st_index] / deltaY;
						}
						if (rInfo->spRefList[j]->geoi->isDomain[Yminus1] == 1) {//downer of membrane
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + Yminus1] += fabs(nuVec[index].ny) * rInfo->srStoichiometry[j] * rpStack[st_index] / deltaY;
						}
						//binding
						if (rInfo->spRefList[j]->geoi->isDomain[index] == 1) {
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + index] += fabs(nuVec[index].ny) * rInfo->srStoichiometry[j] * rpStack[st_index] / (deltaY / 2.0);
						}
					} else if (dimension == 3 && (geoInfo->bType[index].isBofZp || geoInfo->bType[index].isBofZm)) {//z transport
						//transport
						if (rInfo->spRefList[j]->geoi->isDomain[Zplus1] == 1) {//higher of membrane
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + Zplus1] += fabs(nuVec[index].nz) * rInfo->srStoichiometry[j] * rpStack[st_index] / deltaZ;
						}
						if (rInfo->spRefList[j]->geoi->isDomain[Zminus1] == 1) {//lower of membrane
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + Zminus1] += fabs(nuVec[index].nz) * rInfo->srStoichiometry[j] * rpStack[st_index] / deltaZ;
						}
						//binding
						if (rInfo->spRefList[j]->geoi->isDomain[index] == 1) {
							rInfo->spRefList[j]->delta[m * numOfVolIndexes + index] += fabs(nuVec[index].nz) * rInfo->srStoichiometry[j] * rpStack[st_index] / (deltaZ / 2.0);
						}
					}
				}
			}
		}
	}
}

void calcMemDiffusion(variableInfo *sInfo, voronoiInfo *vorI, int Xindex, int Yindex, int Zindex, int m, double dt, int dimension)
{
	int X = 0, Y = 0, Z = 0, index = 0;
	unsigned int i, j;
	int Xplus2 = 0, Xminus2 = 0, Yplus2 = 0, Yminus2 = 0, Zplus2 = 0, Zminus2 = 0;
	int numOfVolIndexes = Xindex * Yindex * Zindex;
	int dcIndex = 0;
	double* val = sInfo->value;
	double* d = sInfo->delta;
	double rk[4] = {0, 0.5, 0.5, 1.0};
	GeometryInfo *geoInfo = sInfo->geoi;
	double area = 0.0;
	//flux
	//2d
	//J = -D * dval / deltaX
	//area = \Sum (d_ij * s_ij) / 4
	//x direction: delta = ((-J / d_ij) * s_ij)/ area
	//3d
	//x direction: delta = (-J * deltaY * deltaZ) / (deltaY * deltaZ * deltaX) = -J / deltaX

	// double Dx = deltaX, Dy = deltaY, Dz = deltaZ;
	for (i = 0; i < geoInfo->domainIndex.size(); i++) {
		index = geoInfo->domainIndex[i];
		Z = index / (Xindex * Yindex);
		Y = (index - Z * Xindex * Yindex) / Xindex;
		X = index - Z * Xindex * Yindex - Y * Xindex;
		Xplus2 = Z * Yindex * Xindex + Y * Xindex + (X + 2);
		Xminus2 = Z * Yindex * Xindex + Y * Xindex + (X - 2);
		Yplus2 = Z * Yindex * Xindex + (Y + 2) * Xindex + X;
		Yminus2 = Z * Yindex * Xindex + (Y - 2) * Xindex + X;
		Zplus2 = (Z + 2) * Yindex * Xindex + Y * Xindex + X;
		Zminus2 = (Z - 2) * Yindex * Xindex + Y * Xindex + X;

		if (sInfo->diffCInfo[0] != 0) {
			if (sInfo->diffCInfo[0]->isUniform == false) dcIndex = index;
			area = 0.0;
			for (j = 0; j < 2; j++) {
				area += vorI[index].diXY[j] * vorI[index].siXY[j];
				area += vorI[index].diYZ[j] * vorI[index].siYZ[j];
				area += vorI[index].diXZ[j] * vorI[index].siXZ[j];
			}
			if (dimension == 2) area /= 2.0;
			else if (dimension == 3) area /= 4.0;
			if (m == 0) {
				//xy plane
				if ((geoInfo->bType[index].isBofXp && geoInfo->bType[index].isBofXm) || (geoInfo->bType[index].isBofYp && geoInfo->bType[index].isBofYm)) {
					for (j = 0; j < 2; j++) {
						sInfo->delta[m * numOfVolIndexes + index] +=
							((sInfo->diffCInfo[0]->value[dcIndex] * (val[vorI[index].adjacentIndexXY[j]] - val[index]) * vorI[index].siXY[j]) / vorI[index].diXY[j]) / area;
					}
				}
				//yz plane (only 3D)
				if (dimension == 3 && (geoInfo->bType[index].isBofYp && geoInfo->bType[index].isBofYm) || (geoInfo->bType[index].isBofZp && geoInfo->bType[index].isBofZm)) {
					for (j = 0; j < 2; j++) {
						sInfo->delta[m * numOfVolIndexes + index] +=
							((sInfo->diffCInfo[0]->value[dcIndex] * (val[vorI[index].adjacentIndexYZ[j]] - val[index]) * vorI[index].siYZ[j]) / vorI[index].diYZ[j]) / area;
					}
				}
				//xz plane (only 3D)
				if (dimension == 3 && (geoInfo->bType[index].isBofXp && geoInfo->bType[index].isBofXm) || (geoInfo->bType[index].isBofZp && geoInfo->bType[index].isBofZm)) {
					for (j = 0; j < 2; j++) {
						sInfo->delta[m * numOfVolIndexes + index] +=
							((sInfo->diffCInfo[0]->value[dcIndex] * (val[vorI[index].adjacentIndexXZ[j]] - val[index]) * vorI[index].siXZ[j]) / vorI[index].diXZ[j]) / area;
					}
				}
			} else {
				//xy plane
				if ((geoInfo->bType[index].isBofXp && geoInfo->bType[index].isBofXm) || (geoInfo->bType[index].isBofYp && geoInfo->bType[index].isBofYm)) {
					for (j = 0; j < 2; j++) {
						sInfo->delta[m * numOfVolIndexes + index]
							+= sInfo->diffCInfo[0]->value[dcIndex] *
							((((val[vorI[index].adjacentIndexXY[j]] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + vorI[index].adjacentIndexXY[j]])
							   - (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])) * vorI[index].siXY[j]) / vorI[index].diXY[j]) / area;
					}
				}
				//yz plane (only 3D)
				if (dimension == 3 && (geoInfo->bType[index].isBofYp && geoInfo->bType[index].isBofYm) || (geoInfo->bType[index].isBofZp && geoInfo->bType[index].isBofZm)) {
					for (j = 0; j < 2; j++) {
						sInfo->delta[m * numOfVolIndexes + index]
							+= sInfo->diffCInfo[0]->value[dcIndex] *
							((((val[vorI[index].adjacentIndexYZ[j]] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + vorI[index].adjacentIndexYZ[j]])
							   - (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])) * vorI[index].siYZ[j]) / vorI[index].diYZ[j]) / area;
					}
				}
				//xz plane (only 3D)
				if (dimension == 3 && (geoInfo->bType[index].isBofXp && geoInfo->bType[index].isBofXm) || (geoInfo->bType[index].isBofZp && geoInfo->bType[index].isBofZm)) {
					for (j = 0; j < 2; j++) {
						sInfo->delta[m * numOfVolIndexes + index]
							+= sInfo->diffCInfo[0]->value[dcIndex] *
							((((val[vorI[index].adjacentIndexXZ[j]] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + vorI[index].adjacentIndexXZ[j]])
							   - (val[index] + rk[m] * dt * d[(m - 1) * numOfVolIndexes + index])) * vorI[index].siXZ[j]) / vorI[index].diXZ[j]) / area;
					}
				}
			}
		}
	}
}
