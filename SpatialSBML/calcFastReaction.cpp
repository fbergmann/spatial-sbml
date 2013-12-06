void calcFastReaction(reactionInfo *rInfo, variableInfo *sInfo, int Xindex, int Yindex, int Zindex, double dt, int m, materialType mType, double stoichimetry)
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
	GeometryInfo *geoInfo = sInfo->geoi;
	int increment = 0;
	ASTNode *ast, *cloneAst, *reCoefficient, *proCoefficient;
	if (sInfo->inVol == true) {
		increment = 2;
	} else {
		increment = 1;
	}

	if (rInfo->reaction->getReversible()) {//equilibrium
		const KineticLaw *kl = r->getKineticLaw();
		ast = const_cast<ASTNode*>(kl->getMath());
		cloneAst = ast->deepCopy();
		/*reactantの係数を調べる
		 *1.astのreactantの値を1に置き換える
		 *2.reactantの最上流の*を見つける
		 *3.見つけた*より下を計算
		 */
		cout << "before:  " << SBML_formulaToString(ast) << endl;

		cout << "after:  " << SBML_formulaToString(ast) << endl;

	} else {//all reactants change to products
		for (Z = 0; Z < Zindex; Z += increment) {
			for (Y = 0; Y < Yindex; Y += increment) {
				for (X = 0; X < Xindex; X += increment) {
					if (sInfo->inVol == true) {//normal reaction
						index = Z * Yindex * Xindex + Y * Xindex + X;
					} else if (sInfo->inVol == false && (Z * Yindex * Xindex + Y * Xindex + X) % 2 != 0) {//membrane trans
					}
				}
			}
		}
	}
}

//astのreactantの値を1に置き換える
void replaceOne(ASTNode *ast)
{
	int i;
	for (i = 0; i < ast->getNumChildren(); i++) {
		replaceOne(ast->getChild(i), rpInfo, varInfoList, index_max, freeConstList);
	}
	if (ast->isName() && ast->getType() == AST_NAME) {
		for (i = 0; i < r->getNumReactants(); i++) {
			if (strcmp(r->getReactants(j)->getSpecies(), ast->getName() == 0)
				variableInfo *info = searchInfoById(varInfoList, ast->getName());
				if (info != 0) {
					if (info->isUniform) {
					} else {
						ast->setType(AST_REAL);
						ast->setValue(1.0);
					}
				}
			}
		}
	}
}
