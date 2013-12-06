//============================================================================
// Name        : SBMLSimulator.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <sbml/xml/XMLError.h>
#include <cstdlib>
#include <getopt.h>
//#include <unistd.h>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/stat.h>
#include "sbml/SBMLTypes.h"
#include "sbml/extension/SBMLExtensionRegistry.h"
#include "sbml/packages/req/common/RequiredElementsExtensionTypes.h"
#include "sbml/packages/spatial/common/SpatialExtensionTypes.h"
#include "sbml/packages/spatial/extension/SpatialModelPlugin.h"
#include "sbml/packages/spatial/extension/SpatialExtension.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "mystruct.h"
#include "initializeFunction.h"
#include "freeFunction.h"
#include "searchFunction.h"
#include "astFunction.h"
#include "calcPDE.h"
#include "setInfoFunction.h"
#include "boundaryFunction.h"
#include "outputFunction.h"
#include "checkStability.h"
#include <zlib.h>

using namespace std;

void spatialSimulator(SBMLDocument *doc, int argc, char *argv[]);
bool isResolvedAll(vector<variableInfo*> &dependence);
void printErrorMessage()
{
	cout << "illegal option" << endl;
	cout << "how to use" << endl;
	cout << "-x #(int): the number of points at x coordinate (ex. -x 100)" << endl;
	cout << "-y #(int): the number of points at y coordinate (ex. -y 100)" << endl;
	cout << "-z #(int): the number of points at z coordinate (ex. -z 100)" << endl;
	cout << "-t #(double or int): simulation time (ex. -t 10)" << endl;
	cout << "-d #(double or int): delta t (ex. -d 0.01)" << endl;
	cout << "-o #(int): output results every # step(ex. -o 10)" << endl;
	cout << "-c #(double or int): max of color bar range # (ex. -c 1)" << endl;
	exit(1);
}

int main(int argc, char *argv[])
{
	clock_t start = clock();
	if (argc == 1) printErrorMessage();

	SBMLDocument *doc = readSBML(argv[argc - 1]);
	//if (doc->getErrorLog()->getNumFailsWithSeverity(LIBSBML_SEV_ERROR) > 0) {
	if (doc->getErrorLog()->contains(XMLFileUnreadable) || doc->getErrorLog()->contains(BadlyFormedXML)
		|| doc->getErrorLog()->contains(MissingXMLEncoding) || doc->getErrorLog()->contains(BadXMLDecl)) {
		doc->printErrors();
		exit(1);
	}

	struct stat st;
	if(stat("./result", &st) != 0) system("mkdir ./result");

	//cout << doc->getNamespaces() << endl;
	//XMLNamespaces *x = doc->getNamespaces();
	if (doc->getModel()->getPlugin("spatial") != 0 && doc->getPkgRequired("spatial") && doc->getPkgRequired("req")) {//PDE
		spatialSimulator(doc, argc, argv);
	} else {//ODE
	}
	SBMLDocument_free(doc);
	clock_t end = clock();
	//cout << endl << "time: " << ((end - start) / static_cast<double>(CLOCKS_PER_SEC)) << " sec" << endl;
	cerr << "time: " << ((end - start) / static_cast<double>(CLOCKS_PER_SEC)) << endl;
    return 0;
}

// int main(int argc, const char* argv[]){
// 	int ret = main2(argc, argv);
// 	char buf[256];
// 	printf("\n");
// 	printf("\n");
// 	sprintf(buf, "leaks %d", getpid());
// 	system(buf);
// 	return ret;
// }

void spatialSimulator(SBMLDocument *doc, int argc, char *argv[])
{
	struct stat st;
	unsigned int i, j, k, m;
	int X = 0, Y = 0, Z = 0, index = 0, divIndex = 0, t = 0, count = 0, file_num = 0, percent = 0, out_step = 1;
	int Xdiv = 101, Ydiv = 101, Zdiv = 101;//num of reaction diffsion mesh
	double end_time = 1.0, dt = 0.01, out_time = 1.0;
	double *sim_time = new double(0.0);
	double deltaX = 0.0, deltaY = 0.0, deltaZ = 0.0;
	double Xsize = 0.0, Ysize = 0.0, Zsize = 0.0;
	int numOfASTNodes = 0;
	//int numOfIndexes = 0;
	char *xaxis = 0, *yaxis = 0, *zaxis = 0;

	//sbml core
	Model *model = doc->getModel();
	ASTNode *ast = 0;
	Species *s;
	SpeciesReference *sr;
	XMLNamespaces *xns = doc->getNamespaces();
	string spatialPrefix = xns->getPrefix("http://www.sbml.org/sbml/level3/version1/spatial/version1");
	string reqPrefix = xns->getPrefix("http://www.sbml.org/sbml/level3/version1/requiredElements/version1");
	ListOfSpecies *los = model->getListOfSpecies();
	ListOfCompartments *loc = model->getListOfCompartments();
	ListOfRules *lorules = model->getListOfRules();
	SpatialCompartmentPlugin *cPlugin = 0;
	//sbml spatial package
	SpatialModelPlugin *spPlugin = static_cast<SpatialModelPlugin*>(model->getPlugin(spatialPrefix));
	Geometry *geometry = spPlugin->getGeometry();

	//size of list
	unsigned int numOfSpecies = static_cast<unsigned int>(model->getNumSpecies());
	unsigned int numOfReactions = static_cast<unsigned int>(model->getNumReactions());
	unsigned int numOfCompartments = static_cast<unsigned int>(model->getNumCompartments());
	unsigned int numOfParameters = static_cast<unsigned int>(model->getNumParameters());
	unsigned int numOfRules = static_cast<unsigned int>(model->getNumRules());

	vector<variableInfo*> varInfoList = vector<variableInfo*>();
	varInfoList.reserve(numOfCompartments + numOfSpecies + numOfParameters);
	vector<GeometryInfo*> geoInfoList = vector<GeometryInfo*>();
	geoInfoList.reserve(geometry->getNumGeometryDefinitions());
	vector<reactionInfo*> rInfoList = vector<reactionInfo*>();
	rInfoList.reserve(numOfReactions + numOfRules);
	vector<reactionInfo*> fast_rInfoList = vector<reactionInfo*>();
	fast_rInfoList.reserve(numOfReactions + numOfRules);
	vector<const char*> memList = vector<const char*>();
	memList.reserve(numOfCompartments);
	unsigned int dimension = geometry->getNumCoordinateComponents();

	vector<double*> freeConstList;
	vector<string> spIdList = vector<string>();
	spIdList.reserve(numOfSpecies);
	for (i = 0; i < numOfSpecies; i++) {
		spIdList.push_back(model->getSpecies(i)->getId());
	}

	int Xplus1 = 0, Xminus1 = 0, Yplus1 = 0, Yminus1 = 0, Zplus1 = 0, Zminus1 = 0;

	//option
	double range_max = 1.0;
	int opt_result = 0;
	extern char	*optarg;
	extern int optind;
	while ((opt_result = getopt(argc - 1, argv, "x:y:z:t:d:o:c:")) != -1) {
		switch(opt_result) {
		case 'x':
			for (i = 0; i < string(optarg).size(); i++) {
				if (!isdigit(optarg[i])) printErrorMessage();
			}
			Xdiv = atoi(optarg) + 1;
			break;
		case 'y':
			for (i = 0; i < string(optarg).size(); i++) {
				if (!isdigit(optarg[i])) printErrorMessage();
			}
			Ydiv = atoi(optarg) + 1;
			break;
		case 'z':
			for (i = 0; i < string(optarg).size(); i++) {
				if (!isdigit(optarg[i])) printErrorMessage();
			}
			Zdiv = atoi(optarg) + 1;
			break;
		case 't':
			for (i = 0; i < string(optarg).size(); i++) {
				if (!isdigit(optarg[i]) && optarg[i] != '.') printErrorMessage();
			}
			end_time = atof(optarg);
			break;
		case 'd':
			for (i = 0; i < string(optarg).size(); i++) {
				if (!isdigit(optarg[i]) && optarg[i] != '.') printErrorMessage();
			}
			dt = atof(optarg);
		case 'o':
			for (i = 0; i < string(optarg).size(); i++) {
				if (!isdigit(optarg[i]) && optarg[i] != '.') printErrorMessage();
			}
			out_step = atoi(optarg);
			break;
		case 'c':
			for (i = 0; i < string(optarg).size(); i++) {
				if (!isdigit(optarg[i]) && optarg[i] != '.') printErrorMessage();
			}
		    range_max = atof(optarg);
			break;
		default:
			printErrorMessage();
			break;
		}
	}

	//div
	if (dimension <= 1) {
		Ydiv = 1;
		Zdiv = 1;
	}
	if (dimension <= 2) {
		Zdiv = 1;
	}

	//filename
	string fname(argv[optind]);
	fname = fname.substr((int)fname.find_last_of("/") + 1, (int)fname.find_last_of(".") - (int)fname.find_last_of("/") - 1);
	if (stat(string("result/" + fname).c_str(), &st) != 0) system(string("mkdir result/" + fname).c_str());

	bool isImageBased = false;
	for (i = 0; i < geometry->getNumGeometryDefinitions(); i++) {
		if (geometry->getGeometryDefinition(i)->isSampledFieldGeometry()) {
			//SampleFieldGeometry
			SampledFieldGeometry *sfGeo	= static_cast<SampledFieldGeometry*>(geometry->getGeometryDefinition(i));
			SampledField *samField = sfGeo->getSampledField();
			cout << "image size:" << endl << "width * height * depth = " << samField->getNumSamples1() << " * " << samField->getNumSamples2() << " * " << samField->getNumSamples3() << endl << endl;
			isImageBased = true;
			Xdiv = samField->getNumSamples1();
			Ydiv = samField->getNumSamples2();
			Zdiv = samField->getNumSamples3();
		}
	}
	out_time = (double)out_time * dt;
	cout << "x mesh num: " << ((isImageBased)? Xdiv: Xdiv - 1) << endl;
	if (dimension >= 2) cout << "y mesh num: " << ((isImageBased)? Ydiv: Ydiv - 1) << endl;
	if (dimension == 3) cout << "z mesh num: " << ((isImageBased)? Zdiv: Zdiv - 1) << endl;
	cout << "simulation time = " << end_time << endl;
	cout << "dt = " << dt << endl;
	cout << "output results every " << out_step << " step" << endl;
	cout << "color bar range max: " << range_max << endl << endl;

	int Xindex = 2 * Xdiv - 1, Yindex = 2 * Ydiv - 1, Zindex = 2 * Zdiv - 1;//num of mesh
	int numOfVolIndexes = Xindex * Yindex * Zindex;
	int indexMax = Zindex * Yindex * Xindex;
	int indexMin = -1;

	//int numOfIndexes = Xdiv * Ydiv * Zdiv;
	//unit
	unsigned int volDimension = 0, memDimension = 0;
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
	t_info->value = sim_time;
	t_info->isResolved = true;
	t_info->isUniform = true;

	//volume index
	vector<int> volumeIndexList;
	for (Z = 0; Z < Zindex; Z += 2) {
		for (Y = 0; Y < Yindex; Y += 2) {
			for (X = 0; X < Xindex; X += 2) {
				volumeIndexList.push_back(Z * Yindex * Xindex + Y * Xindex + X);
			}
		}
	}

	cout << "defining geometry... " << endl;
	//geometryDefinition
	double *tmp_isDomain = new double[numOfVolIndexes];
	for (i = 0; i < geometry->getNumGeometryDefinitions(); i++) {
		if (geometry->getGeometryDefinition(i)->isAnalyticGeometry()) {
			//AnalyticVolumes
			AnalyticGeometry *analyticGeo = static_cast<AnalyticGeometry*>(geometry->getGeometryDefinition(i));
			//gather information of compartment, domainType, analyticVolume
			for (j = 0; j < numOfCompartments; j++) {
				Compartment *c = loc->get(j);
				cPlugin = static_cast<SpatialCompartmentPlugin*>(c->getPlugin(spatialPrefix));
				if (cPlugin != 0) {
					if (AnalyticVolume *analyticVol = analyticGeo->getAnalyticVolume(cPlugin->getCompartmentMapping()->getDomainType())) {
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
						//cout << "before: " << SBML_formulaToString(ast) << endl;
						rearrangeAST(ast);
						//cout << "after: " << SBML_formulaToString(ast) << endl;
						numOfASTNodes = 0;
						countAST(ast, numOfASTNodes);
						geoInfo->rpInfo = new reversePolishInfo();
						geoInfo->rpInfo->varList = new double*[numOfASTNodes];
						fill_n(geoInfo->rpInfo->varList, numOfASTNodes, reinterpret_cast<double*>(0));
						geoInfo->rpInfo->constList = new double*[numOfASTNodes];
						fill_n(geoInfo->rpInfo->constList, numOfASTNodes, reinterpret_cast<double*>(0));
						geoInfo->rpInfo->opfuncList = new int[numOfASTNodes];
						fill_n(geoInfo->rpInfo->opfuncList, numOfASTNodes, 0);
						geoInfo->rpInfo->listNum = numOfASTNodes;
						geoInfo->isDomain = new int[numOfVolIndexes];
						fill_n(geoInfo->isDomain, numOfVolIndexes, 0);
						geoInfo->isBoundary = new int[numOfVolIndexes];
						fill_n(geoInfo->isBoundary, numOfVolIndexes, 0);
						geoInfo->adjacent0 = 0;
						geoInfo->adjacent1 = 0;
						parseAST(ast, geoInfo->rpInfo, varInfoList, numOfASTNodes, freeConstList);
						//judge if the coordinate point is inside the analytic volume
						fill_n(tmp_isDomain, numOfVolIndexes, 0);
						reversePolishInitial(volumeIndexList, geoInfo->rpInfo, tmp_isDomain, numOfASTNodes, Xindex, Yindex, Zindex, false);
						//cout << geoInfo->domainTypeId << endl;
						for (k = 0; k < (unsigned int)numOfVolIndexes; k++) {
							index = k;
							geoInfo->isDomain[k] = (int)tmp_isDomain[k];
							Z = index / (Xindex * Yindex);
							Y = (index - Z * Xindex * Yindex) / Xindex;
							X = index - Z * Xindex * Yindex - Y * Xindex;
							if (geoInfo->isDomain[k] == 1 && string(geoInfo->domainTypeId) == "nucleus") {
								//cerr << X << ", " << Y << ", " << Z << ", " << geoInfo->isDomain[k] << endl;
							} else {
								//cout << geoInfo->domainTypeId << endl;
							}
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
				if (c->getSpatialDimensions() == volDimension) {
					cPlugin = static_cast<SpatialCompartmentPlugin*>(c->getPlugin(spatialPrefix));
					if (cPlugin != 0) {
						SampledField *samField = sfGeo->getSampledField();
						SampledVolume *samVol = 0;
						for (k = 0; k < sfGeo->getNumSampledVolumes(); k++) {
							if (sfGeo->getSampledVolume(k)->getDomainType() == cPlugin->getCompartmentMapping()->getDomainType()) {
								samVol = sfGeo->getSampledVolume(k);
							}
						}
						const ImageData *idata = samField->getImageData();
						uLong uncomprLen = 0;
						switch (dimension) {
						case 1:
							uncomprLen = samField->getNumSamples1();
							break;
						case 2:
							uncomprLen = samField->getNumSamples1() * samField->getNumSamples2();
							break;
						case 3:
							uncomprLen = samField->getNumSamples1() * samField->getNumSamples2() * samField->getNumSamples3();
						default:
							break;
						}
						uLong comprLen = numOfVolIndexes / 4;
						Byte *compr = (Byte*)calloc(sizeof(Byte), comprLen);
						int *compr_int = (int*)calloc(sizeof(int), comprLen);
						Byte *uncompr = (Byte*)calloc(sizeof(Byte), uncomprLen);
						idata->getSamples(compr_int);
						for (k = 0; k < comprLen; k++) {
							compr[k] = (Byte)(compr_int[k]);
						}

						/*
						  #define Z_OK            0
						  171 #define Z_STREAM_END    1
						  172 #define Z_NEED_DICT     2
						  173 #define Z_ERRNO        (-1)
						  174 #define Z_STREAM_ERROR (-2)
						  175 #define Z_DATA_ERROR   (-3)
						  176 #define Z_MEM_ERROR    (-4)
						  177 #define Z_BUF_ERROR    (-5)
						  178 #define Z_VERSION_ERROR (-6)
						*/

						int err = uncompress(uncompr, &uncomprLen, (const Byte*)compr, comprLen);
						if (err != Z_OK) {
							cout << "err with uncompress" << endl;
							cout << err << endl;
							exit(1);
						}

						//ofstream sam_ofs;
						//string sam_fname = "/Users/matsui/Documents/SBMLSimulator/build/Debug/sample.csv";
						//sam_ofs.open(sam_fname.c_str());
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
						geoInfo->isDomain = new int[numOfVolIndexes];
						fill_n(geoInfo->isDomain, numOfVolIndexes, 0);
						geoInfo->isBoundary = new int[numOfVolIndexes];
						fill_n(geoInfo->isBoundary, numOfVolIndexes, 0);
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
									if (static_cast<unsigned int>(uncompr[Z / 2 * samField->getNumSamples2() * samField->getNumSamples1() + (samField->getNumSamples2() - Y / 2 - 1) * samField->getNumSamples1() + X / 2]) == static_cast<unsigned int>(samVol->getSampledValue())) {
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
						free(compr);
						free(compr_int);
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
				for (k = 0; k < analyticGeo->getNumAnalyticVolumes(); k++) {
					AnalyticVolume *analyticVolIn = analyticGeo->getAnalyticVolume(k);
					if (analyticVolEx->getOrdinal() < analyticVolIn->getOrdinal()) {//ex's ordinal is smaller than in's ordinal
						GeometryInfo *geoInfoIn = searchAvolInfoByDomainType(geoInfoList, analyticVolIn->getDomainType().c_str());
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
		SpatialCompartmentPlugin *cPlugin = static_cast<SpatialCompartmentPlugin*>(loc->get(i)->getPlugin(spatialPrefix));
		if (cPlugin == 0) continue;
		DomainType *dType = geometry->getDomainType(cPlugin->getCompartmentMapping()->getDomainType());
		GeometryInfo *geoInfo = searchAvolInfoByCompartment(geoInfoList, loc->get(i)->getId().c_str());
		if (geoInfo == 0) {
			geoInfo = new GeometryInfo;
			InitializeAVolInfo(geoInfo);
			geoInfo->compartmentId = loc->get(i)->getId().c_str();
			geoInfo->domainTypeId = cPlugin->getCompartmentMapping()->getDomainType().c_str();
			//geoInfo->bt = new boundaryType[numOfVolIndexes];
			if (dType->getSpatialDimensions() == volDimension) {
				geoInfo->isVol = true;
			} else if (dType->getSpatialDimensions() == memDimension) {
				geoInfo->isVol = false;
			}
			geoInfoList.push_back(geoInfo);
		}
		for (j = 0; j < geometry->getNumDomains(); j++) {//Domain
			Domain *domain = geometry->getDomain(j);
			if (domain->getDomainType() == geoInfo->domainTypeId) {
				geoInfo->domainIdList.push_back(domain->getSpatialId().c_str());
			}
		}
	}

	for (i = 0; i < geometry->getNumAdjacentDomains(); i++) {//AdjacentDomain
		AdjacentDomains *adDomain = geometry->getAdjacentDomains(i);
		Domain *ad1 = geometry->getDomain(adDomain->getDomain1());
		Domain *ad2 = geometry->getDomain(adDomain->getDomain2());
		GeometryInfo *geoi1 = searchAvolInfoByDomainType(geoInfoList, ad1->getDomainType().c_str());
		GeometryInfo *geoi2 = searchAvolInfoByDomainType(geoInfoList, ad2->getDomainType().c_str());
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
			fill_n(geoInfo->isDomain, numOfVolIndexes, 0);
			geoInfo->isBoundary = new int[numOfVolIndexes];
			fill_n(geoInfo->isBoundary, numOfVolIndexes, 0);
			geoInfo->bType = new boundaryType[numOfVolIndexes];
			for (j = 0; j < (unsigned int)numOfVolIndexes; j++) {
				geoInfo->bType[j].isBofXp = false;
				geoInfo->bType[j].isBofXm = false;
				geoInfo->bType[j].isBofYp = false;
				geoInfo->bType[j].isBofYm = false;
				geoInfo->bType[j].isBofZp = false;
				geoInfo->bType[j].isBofZm = false;
			}
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
	cout << "finished" << endl << endl;
	//make directrise to output result (txt and img)
	if(stat(string("./result/" + fname + "/txt").c_str(), &st) != 0) {
		system(string("mkdir ./result/" + fname + "/txt").c_str());
	}
	if(stat(string("./result/" + fname + "/txt/geometry").c_str(), &st) != 0) {
		system(string("mkdir ./result/" + fname + "/txt/geometry").c_str());
	}
	if(stat(string("./result/" + fname + "/img").c_str(), &st) != 0) {
		system(string("mkdir ./result/" + fname + "/img").c_str());
	}

	for (i = 0; i < numOfSpecies; i++) {
		Species *s = los->get(i);
		variableInfo *sInfo = searchInfoById(varInfoList, s->getId().c_str());
		if (sInfo != 0) {
			if (sInfo->inVol) {
				if(stat(string("./result/" + fname + "/txt/volume").c_str(), &st) != 0) {
					system(string("mkdir ./result/" + fname + "/txt/volume").c_str());
				}
			} else {
				if(stat(string("./result/" + fname + "/txt/membrane").c_str(), &st) != 0) {
					system(string("mkdir ./result/" + fname + "/txt/membrane").c_str());
				}
			}
			if(stat(string("./result/" + fname + "/img/" + s->getId()).c_str(), &st) != 0) {
				system(string("mkdir ./result/" + fname + "/img/" + s->getId()).c_str());
			}
		}
	}

	//draw geometries
	variableInfo *xInfo = 0, *yInfo = 0, *zInfo = 0;
	if (dimension >= 1) {
		xInfo = searchInfoById(varInfoList, xaxis);
	}
	if (dimension >= 2) {
		yInfo = searchInfoById(varInfoList, yaxis);
	}
	if (dimension >= 3) {
		zInfo = searchInfoById(varInfoList, zaxis);
	}

	//set species' initial condition
	//boundary type(Xp, Xm, Yp, Ym, Zx, Zm)
	//parse dependence among species, compartments, parameters
	cout << "assigning initial conditions of species and parameters..." << endl;
	vector<variableInfo*> notOrderedInfo;
	for (i = 0; i < varInfoList.size(); i++) {
		variableInfo *info = varInfoList[i];
		ast = 0;
		if (model->getInitialAssignment(info->id) != 0) {//initial assignment
			if (info->value == 0) {//value is not set yet
				info->value = new double[numOfVolIndexes];
				fill_n(info->value, numOfVolIndexes, 0);
				if (info->sp != 0 && (!info->sp->isSetConstant() || !info->sp->getConstant())) {
					info->delta = new double[4 * numOfVolIndexes];
					fill_n(info->delta, 4 * numOfVolIndexes, 0);
				}
			}
			ast = const_cast<ASTNode*>((model->getInitialAssignment(info->id))->getMath());
		} else if (model->getRule(info->id) != 0 && model->getRule(info->id)->isAssignment()) {//assignment rule
			info->hasAssignmentRule = true;
			if (info->value == 0) {//value is not set yet
				info->value = new double[numOfVolIndexes];
				fill_n(info->value, numOfVolIndexes, 0);
				if (info->sp != 0 && (!info->sp->isSetConstant() || !info->sp->getConstant())) {
					//the species is variable
					info->delta = new double[4 * numOfVolIndexes];
					fill_n(info->delta, 4 * numOfVolIndexes, 0);
				}
			}
			ast = const_cast<ASTNode*>(((AssignmentRule*)model->getRule(info->id))->getMath());
		}
		if (ast != 0) {
			//cout << info->id << ": " << SBML_formulaToString(ast) << endl;
			rearrangeAST(ast);
			numOfASTNodes = 0;
			countAST(ast, numOfASTNodes);
			info->rpInfo = new reversePolishInfo();
			info->rpInfo->varList = new double*[numOfASTNodes];
			fill_n(info->rpInfo->varList, numOfASTNodes, reinterpret_cast<double*>(0));
			info->rpInfo->deltaList = 0;
			info->rpInfo->constList = new double*[numOfASTNodes];
			fill_n(info->rpInfo->constList, numOfASTNodes, reinterpret_cast<double*>(0));
			info->rpInfo->opfuncList = new int[numOfASTNodes];
			fill_n(info->rpInfo->opfuncList, numOfASTNodes, 0);
			info->rpInfo->listNum = numOfASTNodes;
			info->isResolved = false;
			parseDependence(ast, info->dependence, varInfoList);
			notOrderedInfo.push_back(info);
		}
	}
	//dependency of symbols
	unsigned int resolved_count = 0;
	vector<variableInfo*> orderedARule;
	if (notOrderedInfo.size() != 0) {
		for (i = 0;;i++) {
			variableInfo *info = notOrderedInfo[i];
			if (isResolvedAll(info->dependence) && info->isResolved == false) {
				if (model->getInitialAssignment(info->id) != 0) {//initial assignment
					ast = const_cast<ASTNode*>((model->getInitialAssignment(info->id))->getMath());
				} else if (model->getRule(info->id) != 0 && model->getRule(info->id)->isAssignment()) {//assignment rule
					ast = const_cast<ASTNode*>(((AssignmentRule*)model->getRule(info->id))->getMath());
				}
				parseAST(ast, info->rpInfo, varInfoList, info->rpInfo->listNum, freeConstList);
				cout << info->id << ": " << SBML_formulaToString(ast) << endl;
				bool isAllArea = (info->sp != 0)? false: true;
				if (info->sp != 0) info->geoi = searchAvolInfoByCompartment(geoInfoList, info->sp->getCompartment().c_str());
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
	cout << "finished" << endl;

	//calc normal unit vector of membrane (for mem diffusion and mem transport)
	normalUnitVector *nuVec = 0;
	voronoiInfo *vorI = 0;
	if (dimension >= 2) {
		//calc normalUnitVector at membrane
		nuVec = setNormalAngle(geoInfoList, Xsize, Ysize, Zsize, dimension, Xindex, Yindex, Zindex, numOfVolIndexes);
		//calc voronoi at membrane
		vorI = setVoronoiInfo(nuVec, xInfo, yInfo, zInfo, geoInfoList, Xsize, Ysize, Zsize, dimension, Xindex, Yindex, Zindex, numOfVolIndexes);
	}

	//set boundary type
	setBoundaryType(model, varInfoList, geoInfoList, Xindex, Yindex, Zindex, dimension);

	//numerical stability analysis of diffusion and advection
	double min_dt = dt;
	cout << endl << "checking numerical stability of diffusion and advection... " << endl;
	for (i = 0; i < numOfSpecies; i++) {
		variableInfo *sInfo = searchInfoById(varInfoList, los->get(i)->getId().c_str());
		//volume diffusion
		if (sInfo->diffCInfo != 0 && sInfo->geoi->isVol) {
			min_dt = min(min_dt, checkDiffusionStab(sInfo, deltaX, deltaY, deltaZ, Xindex, Yindex, dt));
		}
		//membane diffusion
		if (sInfo->diffCInfo != 0 && !sInfo->geoi->isVol) {
			min_dt = min(min_dt, checkMemDiffusionStab(sInfo, vorI, Xindex, Yindex, dt, dimension));
		}
		//advection
		if (sInfo->adCInfo != 0) {
			min_dt = min(min_dt, checkAdvectionStab(sInfo, deltaX, deltaY, deltaZ, dt, Xindex, Yindex, dimension));
		}
	}
	cout << "finished" << endl;
	if (dt > min_dt) {
		cout << "dt must be less than " << min_dt << endl;
		exit(1);
	}

	//reaction information
	setReactionInfo(model, varInfoList, rInfoList, fast_rInfoList, freeConstList, numOfVolIndexes);

	//rate rule information
	setRateRuleInfo(model, varInfoList, rInfoList, freeConstList, numOfVolIndexes);

	//output geometries
	cout << endl << "outputting geometries into text file... " << endl;
	int *geo_edge = new int[numOfVolIndexes];
	fill_n(geo_edge, numOfVolIndexes, 0);
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

	ofstream ofs;
	string geo_filename = "./result/" + fname + "/txt/geometry/all_membrane.csv";
	ofs.open(geo_filename.c_str());
	ofs << "# information about membrane domain" << endl;
	ofs << "# about values of domainTypeId column: 0 - not in the domain, 1 - in the domain, 2 - in the pseudo domain" << endl;
	//ofs << "# mesh num: [" << (Xindex + 1) / 2 - 1 << " x " << (Yindex + 1) / 2 - 1 << " x " << (Zindex + 1) / 2 - 1 << "]" << endl << endl;
	ofs << "# x";
	switch(dimension) {
	case 1:
		ofs << ", all";
		for (i = 0; i < memList.size(); i++) {
			ofs << ", " << memList[i];
		}
		ofs << endl;
		for (X = 0; X < Xindex; X++) {
			//all membrane
			ofs << xInfo->value[X] << ", " << geo_edge[X];
			//each membrane domain
			for (i = 0; i < memInfoList.size(); i++) {
				ofs << ", " << memInfoList[i]->isDomain[index];
			}
			ofs << endl;
		}
		break;
	case 2:
		ofs << ", y, all";
		for (i = 0; i < memList.size(); i++) {
			ofs << ", " << memList[i];
		}
		ofs << endl;
		for (Y = 0; Y < Yindex; Y++) {
			for (X = 0; X < Xindex; X++) {
				index = Y * Xindex + X;
				//all membrane
				ofs << xInfo->value[index] << ", " << yInfo->value[index] << ", " << geo_edge[index];
				//each membrane domain
				for (i = 0; i < memInfoList.size(); i++) {
					ofs << ", " << memInfoList[i]->isDomain[index];
				}
				ofs << endl;
			}
			ofs << endl;
		}
		break;
	case 3:
		ofs << ", y, z, all";
		for (i = 0; i < memList.size(); i++) {
			ofs << ", " << memList[i];
		}
		ofs << endl;
		for (Z = 0; Z < Zindex; Z++) {
			for (Y = 0; Y < Yindex; Y++) {
				for (X = 0; X < Xindex; X++) {
					index = Z * Yindex * Xindex + Y * Xindex + X;
					//all membrane
					ofs << xInfo->value[index] << ", " << yInfo->value[index] << ", " << zInfo->value[index] << ", " << geo_edge[index];
					//each membrane domain
					for (i = 0; i < memInfoList.size(); i++) {
						ofs << ", " << memInfoList[i]->isDomain[index];
					}
					ofs << endl;
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
	cout << "finished" << endl << endl;

	//simulation
	cout << "simulation starts" << endl;
	FILE *gp = popen("/opt/local/bin/gnuplot -persist", "w");
	clock_t diff_start, diff_end, boundary_start, boundary_end, out_start, out_end, re_start, re_end, ad_start, ad_end, assign_start, assign_end, update_start, update_end;
	clock_t re_time = 0, diff_time = 0, output_time = 0, ad_time = 0, update_time = 0, mem_time = 0, boundary_time = 0, assign_time = 0;
	clock_t sim_start = clock();
	cout << endl;
	for (t = 0; t <= static_cast<int>(end_time / dt); t++) {
		//cerr << static_cast<int>(100.0 * static_cast<double>(t) / (end_time / dt)) << endl;
		*sim_time = t * dt;
		//output
		out_start = clock();
		if (count % out_step == 0) {
      // funa
			outputTimeCource(gp, model, varInfoList, memList, xInfo, yInfo, zInfo, sim_time, end_time, dt, range_max, dimension, Xindex, Yindex, Zindex, Xsize, Ysize, Zsize, file_num, fname);
			file_num++;
		}
		out_end = clock();
		output_time += out_end - out_start;
		count++;

		//calculation
		//advection
		ad_start = clock();
		for (i = 0; i < numOfSpecies; i++) {
			variableInfo *sInfo = searchInfoById(varInfoList, los->get(i)->getId().c_str());
			//advection
			if (sInfo->adCInfo != 0) {
				cipCSLR(sInfo, deltaX, deltaY, deltaZ, dt, Xindex, Yindex, Zindex, dimension);
			}//end of advection
		}
		ad_end = clock();
		ad_time += ad_end - ad_start;
		//cout << "ad time: "<< ((ad_end - ad_start) / static_cast<double>(CLOCKS_PER_SEC)) << endl;

		//runge-kutta
		for (m = 0; m < 4; m++) {
			//diffusion
			for (i = 0; i < numOfSpecies; i++) {
				variableInfo *sInfo = searchInfoById(varInfoList, los->get(i)->getId().c_str());
				diff_start = clock();
				//volume diffusion
				if (sInfo->diffCInfo != 0 && sInfo->geoi->isVol) {
					calcDiffusion(sInfo, deltaX, deltaY, deltaZ, Xindex, Yindex, Zindex, m, dt);
				}
				//membane diffusion
				if (sInfo->diffCInfo != 0 && !sInfo->geoi->isVol) {
					calcMemDiffusion(sInfo, vorI, Xindex, Yindex, Zindex, m, dt, dimension);
				}
				diff_end = clock();
				diff_time += diff_end - diff_start;
				boundary_start = clock();
				//boundary condition
				if (sInfo->boundaryInfo != 0 && sInfo->geoi->isVol) {
					calcBoundary(sInfo, deltaX, deltaY, deltaZ, Xindex, Yindex, Zindex, m, dimension);
				}
				boundary_end = clock();
				boundary_time += (boundary_end - boundary_start);
			}
			//reaction
			re_start = clock();
			//for (i = 0; i < numOfReactions; i++) {
			//slow reaction
			for (i = 0; i < rInfoList.size(); i++) {
				//Reaction *r = model->getReaction(i);
				Reaction *r = rInfoList[i]->reaction;
				if (!rInfoList[i]->isMemTransport) {//normal reaction
					sr = r->getReactant(0);
					reversePolishRK(rInfoList[i], searchInfoById(varInfoList, sr->getSpecies().c_str())->geoi, Xindex, Yindex, Zindex, dt, m, r->getNumReactants(), true);
				} else {//membrane transport
					GeometryInfo *reactantGeo = searchInfoById(varInfoList, r->getReactant(0)->getSpecies().c_str())->geoi;
					GeometryInfo *productGeo = searchInfoById(varInfoList, r->getProduct(0)->getSpecies().c_str())->geoi;
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
					if (reactantGeo->isVol ^ productGeo->isVol) {
						if (!reactantGeo->isVol) {
							calcMemTransport(rInfoList[i], reactantGeo, nuVec, Xindex, Yindex, Zindex, dt, m, deltaX, deltaY, deltaZ, dimension, r->getNumReactants());
						}
						if (!productGeo->isVol) {
							calcMemTransport(rInfoList[i], productGeo, nuVec, Xindex, Yindex, Zindex, dt, m, deltaX, deltaY, deltaZ, dimension, r->getNumReactants());
						}
					}
				}
			}
			re_end = clock();
			re_time += (re_end - re_start);
			//rate rule
			for (i = 0; i < numOfRules; i++) {
				if (model->getRule(i)->isRate()) {
					RateRule *rrule = (RateRule*)model->getRule(i);
					variableInfo *sInfo = searchInfoById(varInfoList, rrule->getVariable().c_str());
 					reversePolishRK(rInfoList[i], sInfo->geoi, Xindex, Yindex, Zindex, dt, m, 1, false);
				}
			}
		}//end of runge-kutta
		//update values (advection, diffusion, slow reaction)
		update_start = clock();
		for (i = 0; i < numOfSpecies; i++) {
			s = los->get(i);
			variableInfo *sInfo = searchInfoById(varInfoList, s->getId().c_str());
			if (!s->isSetConstant() || !s->getConstant()) {
				for (j = 0; j < sInfo->geoi->domainIndex.size(); j++) {
					index = sInfo->geoi->domainIndex[j];
					Z = index / (Xindex * Yindex);
					Y = (index - Z * Xindex * Yindex) / Xindex;
					X = index - Z * Xindex * Yindex - Y * Xindex;
					divIndex = (Z / 2) * Ydiv * Xdiv + (Y / 2) * Xdiv + (X / 2);
					//update values for the next time
					sInfo->value[index] += dt * (sInfo->delta[index] + 2.0 * sInfo->delta[numOfVolIndexes + index] + 2.0 * sInfo->delta[2 * numOfVolIndexes + index] + sInfo->delta[3 * numOfVolIndexes + index]) / 6.0;
					for (k = 0; k < 4; k++) sInfo->delta[k * numOfVolIndexes + index] = 0.0;
				}
				//boundary condition
				if (sInfo->boundaryInfo != 0) {
					calcBoundary(sInfo, deltaX, deltaY, deltaZ, Xindex, Yindex, Zindex, 0, dimension);
				}
			}
		}
		update_end = clock();
		update_time += update_end - update_start;

		//fast reaction
		// 		for (i = 0; i < fast_rInfoList.size(); i++) {
		// 			Reaction *r = fast_rInfoList[i]->reaction;
		// 		}
		//assignment rule
		assign_start = clock();
		for (i = 0; i < orderedARule.size(); i++) {
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
		assign_end = clock();
		assign_time += assign_end - assign_start;
		//pseudo membrane
		clock_t mem_start = clock();
		for (i = 0; i < numOfSpecies; i++) {
			s = los->get(i);
			variableInfo *sInfo = searchInfoById(varInfoList, s->getId().c_str());
			if (!sInfo->geoi->isVol) {
				for (j = 0; j < sInfo->geoi->pseudoMemIndex.size(); j++) {
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
		clock_t mem_end = clock();
		mem_time += mem_end - mem_start;

		if (t == (static_cast<int>(end_time / dt) / 10) * percent) {
			cout << percent * 10 << "% finished" << endl;
			percent++;
		}
	}
	clock_t sim_end = clock();
	cerr << "simulation_time: "<< ((sim_end - sim_start) / static_cast<double>(CLOCKS_PER_SEC)) << endl;
	cerr << "reaction_time: "<< (re_time / static_cast<double>(CLOCKS_PER_SEC)) << endl;
	cerr << "diffusion_time: "<< (diff_time / static_cast<double>(CLOCKS_PER_SEC)) << endl;
	cerr << "advection_time: "<< (ad_time / static_cast<double>(CLOCKS_PER_SEC)) << endl;
	cerr << "update_time: "<< (update_time / static_cast<double>(CLOCKS_PER_SEC)) << endl;
	//cerr << "  mem_time: "<< (mem_time / static_cast<double>(CLOCKS_PER_SEC)) << endl;
	cerr << "assign_time: "<< (assign_time / static_cast<double>(CLOCKS_PER_SEC)) << endl;
	cerr << "output_time: "<< (output_time / static_cast<double>(CLOCKS_PER_SEC)) << endl;
	pclose(gp);

	//free
	freeVarInfo(varInfoList);
	freeAvolInfo(geoInfoList);
	freeRInfo(rInfoList);
	for (i = 0; i < freeConstList.size(); i++) {
		delete freeConstList[i];
		freeConstList[i] = 0;
	}
	delete sim_time;
	delete[] nuVec;
	delete[] vorI;
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
