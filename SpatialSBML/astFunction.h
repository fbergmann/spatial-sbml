#ifndef ASTFUNCTION_H_
#define ASTFUNCTION_H_

void rearrangeAST(ASTNode *ast);
void countAST(ASTNode *ast, int &numOfASTNode);
void parseAST(ASTNode *ast, reversePolishInfo *rpInfo, vector<variableInfo*> &varInfoList, int index_max, vector<double*> &freeConstList);
void parseDependence(const ASTNode *ast, vector<variableInfo*> &dependence, vector<variableInfo*> &varInfoList);

#endif
