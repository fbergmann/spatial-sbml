#include "sbml/SBMLTypes.h"
#include "sbml/extension/SBMLExtensionRegistry.h"
#include "sbml/packages/spatial/common/SpatialExtensionTypes.h"
#include "sbml/packages/spatial/extension/SpatialModelPlugin.h"
#include "sbml/packages/spatial/extension/SpatialExtension.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include "mystruct.h"
#include "searchFunction.h"

void setBoundaryType(Model *model, vector<variableInfo*> &varInfoList, vector<GeometryInfo*> &geoInfoList, int Xindex, int Yindex, int Zindex, int dimension)
{
	ListOfSpecies *los = model->getListOfSpecies();
	unsigned int i;
	int X, Y, Z, index;
	unsigned int numOfSpecies = static_cast<unsigned int>(model->getNumSpecies());
	int Xplus1 = 0, Xminus1 = 0, Yplus1 = 0, Yminus1 = 0, Zplus1 = 0, Zminus1 = 0;
	int Xplus2 = 0, Xminus2 = 0, Yplus2 = 0, Yminus2 = 0, Zplus2 = 0, Zminus2 = 0;
	for (i = 0; i < numOfSpecies; i++) {
		Species *s = los->get(i);
		variableInfo *sInfo = searchInfoById(varInfoList, s->getId().c_str());
    if (sInfo != 0 &&  !sInfo->isUniform && searchAvolInfoByCompartment(geoInfoList, s->getCompartment().c_str()) != 0) {
    	sInfo->geoi = searchAvolInfoByCompartment(geoInfoList, s->getCompartment().c_str());
			int *isD = sInfo->geoi->isDomain;
			int *isB = sInfo->geoi->isBoundary;

      if (isD == NULL || isB == NULL)
        continue;

			switch (dimension) {
			case 1:
				if (sInfo->geoi->isVol == true) {
					for (X = 0; X < Xindex; X += 2) {
						if (isD[X] == 0) sInfo->value[X] = 0.0;
						if (isD[X] == 1 && isB[X] == 1) {
							if (X == Xindex - 1 || X == 0) {
								if (X == 0) sInfo->geoi->bType[X].isBofXm = true;
								if (X == Xindex - 1) sInfo->geoi->bType[X].isBofXp = true;
							} else {//not the edge of simulation area
								if (isD[X - 2] == 0) sInfo->geoi->bType[X].isBofXm = true;
								if (isD[X + 2] == 0) sInfo->geoi->bType[X].isBofXp = true;
							}
						}
					}
				} else {//membrane
					for (X = 0; X < Xindex; X++) {
						Xplus1 = X + 1;
						Xminus1 = X - 1;
						if (isD[X] == 0) sInfo->value[X] = 0.0;
						if (X % 2 != 0) {
							if (isD[X] == 2 && isD[Xplus1] == 1 && isD[Xminus1] == 1) {
								sInfo->value[X] = min(sInfo->value[Xplus1], sInfo->value[Xminus1]);
							}
						}
					}
				}
				break;
			case 2:
				if (sInfo->geoi->isVol == true) {//volume
					for (Y = 0; Y < Yindex; Y += 2) {
						for (X = 0; X < Xindex; X += 2) {
							index = Y * Xindex + X;
							Xplus2 = Y * Xindex + (X + 2);
							Xminus2 = Y * Xindex + (X - 2);
							Yplus2 = (Y + 2) * Xindex + X;
							Yminus2 = (Y - 2) * Xindex + X;
							//Xplus1 = Y * Xindex + (X + 1);
							//Xminus1 = Y * Xindex + (X - 1);
							//Yplus1 = (Y + 1) * Xindex + X;
							//Yminus1 = (Y - 1) * Xindex + X;
							if (isD[index] == 0) sInfo->value[index] = 0.0;
							if (isD[index] == 1 && isB[index] == 1) {
								if (X == Xindex - 1 || X == 0 || Y == Yindex - 1 || Y == 0) {
									if (X == 0) sInfo->geoi->bType[index].isBofXm = true;
									if (X == Xindex - 1) sInfo->geoi->bType[index].isBofXp = true;
									if (Y == 0) sInfo->geoi->bType[index].isBofYm = true;
									if (Y == Yindex - 1) sInfo->geoi->bType[index].isBofYp = true;
								} else {//not the edge of simulation area
									if (isD[Xminus2] == 0) sInfo->geoi->bType[index].isBofXm = true;
									if (isD[Xplus2] == 0) sInfo->geoi->bType[index].isBofXp = true;
									if (isD[Yminus2] == 0) sInfo->geoi->bType[index].isBofYm = true;
									if (isD[Yplus2] == 0) sInfo->geoi->bType[index].isBofYp = true;
								}
							}
						}
					}
				} else {//membrane(角の値を補完)
					for (Y = 0; Y < Yindex; Y++) {
						for (X = 0; X < Xindex; X++) {
							index = Y * Xindex + X;
							Xplus1 = Y * Xindex + (X + 1);
							Xminus1 = Y * Xindex + (X - 1);
							Yplus1 = (Y + 1) * Xindex + X;
							Yminus1 = (Y - 1) * Xindex + X;
							if (isD[index] == 0) sInfo->value[index] = 0.0;
							if (X % 2 != 0 && Y % 2 != 0) {
								//values at pseudo membrane
								if (isD[index] == 2) {
									if (isD[Xplus1] == 1 && isD[Xminus1] == 1) {
										sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Xminus1]);
									} else if (isD[Yplus1] == 1 && isD[Yminus1] == 1) {
										sInfo->value[index] = min(sInfo->value[Yplus1], sInfo->value[Yminus1]);
									} else if (isD[Xplus1] == 1 && isD[Yplus1] == 1) {
										sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Yplus1]);
									} else if (isD[Xplus1] == 1 && isD[Yminus1] == 1) {
										sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Yminus1]);
									} else if (isD[Xminus1] == 1 && isD[Yplus1] == 1) {
										sInfo->value[index] = min(sInfo->value[Xminus1], sInfo->value[Yplus1]);
									} else if (isD[Xminus1] == 1 && isD[Yminus1] == 1) {
										sInfo->value[index] = min(sInfo->value[Xminus1], sInfo->value[Yminus1]);
									}
								}
							}
						}
					}
				}
				break;
			case 3:
				if (sInfo->geoi->isVol == true) {//volume
					for (Z = 0; Z < Zindex; Z += 2) {
						for (Y = 0; Y < Yindex; Y += 2) {
							for (X = 0; X < Xindex; X += 2) {
								index = Z * Yindex * Xindex + Y * Xindex + X;
								Xplus2 = Z * Yindex * Xindex + Y * Xindex + (X + 2);
								Xminus2 = Z * Yindex * Xindex + Y * Xindex + (X - 2);
								Yplus2 = Z * Yindex * Xindex + (Y + 2) * Xindex + X;
								Yminus2 = Z * Yindex * Xindex + (Y - 2) * Xindex + X;
								Zplus2 = (Z + 2) * Yindex * Xindex + Y * Xindex + X;
								Zminus2 = (Z - 2) * Yindex * Xindex + Y * Xindex + X;
								//Xplus1 = Z * Yindex * Xindex + Y * Xindex + (X + 1);
								//Xminus1 = Z * Yindex * Xindex + Y * Xindex + (X - 1);
								//Yplus1 = Z * Yindex * Xindex + (Y + 1) * Xindex + X;
								//Yminus1 = Z * Yindex * Xindex + (Y - 1) * Xindex + X;
								//Zplus1 = (Z + 1) * Yindex * Xindex + Y * Xindex + X;
								//Zminus1 = (Z - 1) * Yindex * Xindex + Y * Xindex + X;
								if (isD[index] == 0) sInfo->value[index] = 0.0;
								if (isD[index] == 1 && isB[index] == 1) {
									if (X == Xindex - 1 || X == 0 || Y == Yindex - 1 || Y == 0 || Z == Zindex - 1 || Z == 0) {
										if (X == 0) sInfo->geoi->bType[index].isBofXm = true;
										if (X == Xindex - 1) sInfo->geoi->bType[index].isBofXp = true;
										if (Y == 0) sInfo->geoi->bType[index].isBofYm = true;
										if (Y == Yindex - 1) sInfo->geoi->bType[index].isBofYp = true;
										if (Z == 0) sInfo->geoi->bType[index].isBofZm = true;
										if (Z == Zindex - 1) sInfo->geoi->bType[index].isBofZp = true;
									} else {//not the edge of simulation area
										if (isD[Xminus2] == 0) sInfo->geoi->bType[index].isBofXm = true;
										if (isD[Xplus2] == 0) sInfo->geoi->bType[index].isBofXp = true;
										if (isD[Yminus2] == 0) sInfo->geoi->bType[index].isBofYm = true;
										if (isD[Yplus2] == 0) sInfo->geoi->bType[index].isBofYp = true;
										if (isD[Zminus2] == 0) sInfo->geoi->bType[index].isBofZm = true;
										if (isD[Zplus2] == 0) sInfo->geoi->bType[index].isBofZp = true;
									}
								}
							}
						}
					}
				} else {//membrane
					for (Z = 0; Z < Zindex; Z++) {
						for (Y = 0; Y < Yindex; Y++) {
							for (X = 0; X < Xindex; X++) {
								index = Z * Yindex * Xindex + Y * Xindex + X;
								Xplus1 = Z * Yindex * Xindex + Y * Xindex + (X + 1);
								Xminus1 = Z * Yindex * Xindex + Y * Xindex + (X - 1);
								Yplus1 = Z * Yindex * Xindex + (Y + 1) * Xindex + X;
								Yminus1 = Z * Yindex * Xindex + (Y - 1) * Xindex + X;
								Zplus1 = (Z + 1) * Yindex * Xindex + Y * Xindex + X;
								Zminus1 = (Z - 1) * Yindex * Xindex + Y * Xindex + X;
								if (isD[index] == 0) sInfo->value[index] = 0.0;
								if ((X % 2 != 0 && Y % 2 != 0) || (Y % 2 != 0 && Z % 2 != 0) || (Z % 2 != 0 && X % 2 != 0)) {
									//values at pseudo membrane
									if (isD[index] == 2) {
										if (isD[Xplus1] == 1 && isD[Xminus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Xminus1]);
										} else if (isD[Yplus1] == 1 && isD[Yminus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Yplus1], sInfo->value[Yminus1]);
										} else if (isD[Zplus1] == 1 && isD[Zminus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Zplus1], sInfo->value[Zminus1]);
										} else if (isD[Xplus1] == 1 && isD[Yplus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Yplus1]);
										} else if (isD[Xplus1] == 1 && isD[Yminus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Yminus1]);
										} else if (isD[Xplus1] == 1 && isD[Zplus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Zplus1]);
										} else if (isD[Xplus1] == 1 && isD[Zminus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Xplus1], sInfo->value[Zminus1]);
										} else if (isD[Xminus1] == 1 && isD[Yplus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Xminus1], sInfo->value[Yplus1]);
										} else if (isD[Xminus1] == 1 && isD[Yminus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Xminus1], sInfo->value[Yminus1]);
										} else if (isD[Xminus1] == 1 && isD[Zplus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Xminus1], sInfo->value[Zplus1]);
										} else if (isD[Xminus1] == 1 && isD[Zminus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Xminus1], sInfo->value[Zminus1]);
										} else if (isD[Yplus1] == 1 && isD[Zplus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Yplus1], sInfo->value[Zplus1]);
										} else if (isD[Yplus1] == 1 && isD[Zminus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Yplus1], sInfo->value[Zminus1]);
										} else if (isD[Yminus1] == 1 && isD[Zplus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Yminus1], sInfo->value[Zplus1]);
										} else if (isD[Yminus1] == 1 && isD[Zminus1] == 1) {
											sInfo->value[index] = min(sInfo->value[Yminus1], sInfo->value[Zminus1]);
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
