
int old_main(int argc, const char *argv[])
{
  clock_t start = clock();

  if (argc < 2)
  {
    cerr << "Need one argument: sbml file with spatial simulation" << endl;
    exit(1);
  }
  string fileName (argv[1]);
  cout << "loading: "<< fileName << endl;
  SBMLDocument *doc = readSBML(fileName.c_str());

  if (doc->getModel() == NULL)
  {
    cerr << "encountered fatal errors while reading the file." << endl;
    doc->printErrors();
    exit(1);
  }

  //cout << doc->getNamespaces() << endl;
  //XMLNamespaces *x = doc->getNamespaces();
  if (doc->getModel()->getPlugin("spatial") != 0 && doc->getPkgRequired("spatial") && doc->getPkgRequired("req")) {//PDE
    cout << "found spatial model ... "<< endl;
    spatialSimulator(doc, argc, argv);
  } else {//ODE
    cout << "no spatial model ... " << endl;
  }
  SBMLDocument_free(doc);
  clock_t end = clock();
  //cout << endl << "time: " << ((end - start) / static_cast<double>(CLOCKS_PER_SEC)) << " sec" << endl;
  cerr << "time: " << ((end - start) / static_cast<double>(CLOCKS_PER_SEC)) << endl;
  return 0;
}

void spatialSimulator(SBMLDocument *doc, int argc, const char *argv[])
{
  struct stat st;
  unsigned int i, j, k, m;
  int X = 0, Y = 0, Z = 0, index = 0, t = 0, count = 0, file_num = 0, percent = 0;
  int Xdiv = 1, Ydiv = 1, Zdiv = 1;//cin
  double end_time = 0.0, dt = 0.0, sim_time = 0.0, out_time = 1.0;
  double deltaX = 0.0, deltaY = 0.0, deltaZ = 0.0;
  int numOfASTNodes = 0;
  //int numOfIndexes = 0;
  char *xaxis = 0, *yaxis = 0, *zaxis = 0;

#ifdef WIN32
  string gnuplotExecutable("C:/gnuplot/bin/gnuplot.exe");
#else
  string gnuplotExecutable("gnuplot");
#endif

  //sbml core
  Model *model = doc->getModel();
  ASTNode *ast;
  XMLNamespaces *xns = doc->getNamespaces();
  string spatialPrefix = xns->getPrefix("http://www.sbml.org/sbml/level3/version1/spatial/version1");
  string reqPrefix = xns->getPrefix("http://www.sbml.org/sbml/level3/version1/requiredElements/version1");
  ListOfSpecies *los = model->getListOfSpecies();
  ListOfCompartments *loc = model->getListOfCompartments();
  ListOfReactions *lor = model->getListOfReactions();
  ListOfParameters *lop = model->getListOfParameters();
  //ListOfRules *lorules = model->getListOfRules();

  //sbml spatial package
  SpatialModelPlugin *spPlugin = static_cast<SpatialModelPlugin*>(model->getPlugin("spatial"));
  Geometry *geometry = spPlugin->getGeometry();
  //ListOfCoordinateComponents *locc = geometry->getListOfCoordinateComponents();

  //size of list
  unsigned int numOfSpecies = static_cast<unsigned int>(model->getNumSpecies());
  unsigned int numOfReactions = static_cast<unsigned int>(model->getNumReactions());
  unsigned int numOfCompartments = static_cast<unsigned int>(model->getNumCompartments());
  unsigned int numOfParameters = static_cast<unsigned int>(model->getNumParameters());
  unsigned int numOfRules = static_cast<unsigned int>(model->getNumRules());

  vector<variableInfo*> varInfoList = vector<variableInfo*>();
  vector<bcOfSpeciesInfo*> bcOfSpeciesInfoList = vector<bcOfSpeciesInfo*>();
  vector<analyticVolInfo*> avolInfoList = vector<analyticVolInfo*>();
  vector<boundaryCInfo*> bcInfoList = vector<boundaryCInfo*>();
  vector<reactionInfo*> rInfoList = vector<reactionInfo*>();
  vector<adCInfo*> acInfoList = vector<adCInfo*>();
  vector<const char*> memList = vector<const char*>();
  unsigned int dimension = geometry->getNumCoordinateComponents();

  int Xplus = 0, Xminus = 0, Yplus = 0, Yminus = 0, Zplus = 0, Zminus = 0;

  if(stat("./result", &st) != 0) {    
    MKDIR("./result");
  }

  if(stat("./result/txt", &st) != 0) {    

    MKDIR("./result/txt");
  }
  if(stat("./result/png", &st) != 0) {
    MKDIR("./result/png");
  }

  //div
  if (dimension <= 1) {
    Ydiv = 1;
    Zdiv = 1;
  }
  if (dimension <= 2) {
    Zdiv = 1;
  }

  if (dimension >= 1) {
    cout << "x mesh num: ";
    cin >> Xdiv;
  }
  if (dimension >= 2) {
    cout << "y mesh num: ";
    cin >> Ydiv;
  }
  if (dimension >= 3) {
    cout << "z mesh num: ";
    cin >> Zdiv;
  }
  cout << "end time: ";
  cin >> end_time;
  cout << "dt: ";
  cin >> dt;

  int Xindex = 2 * Xdiv - 1, Yindex = 2 * Ydiv - 1, Zindex = 2 * Zdiv - 1;
  int numOfVolIndexes = Xindex * Yindex * Zindex;
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


  // 	//volume
  // 	for (Y = 0; Y < Yindex; Y += 2) {
  // 		for (X = 0; X < Xindex; X += 2) {
  // 			cout << (Y * Xindex + X) << " ";
  // 		}
  // 	}
  // 	cout << endl;
  // 	//membrane
  // 	for (Y = 0; Y < Yindex; Y++) {
  // 		for (X = 0; X < Xindex; X++) {
  // 			if (!(X % 2 == 0 && Y % 2 == 0)) {
  // 				cout << (Y * Xindex + X) << " ";
  // 			}
  // 		}
  // 	}
  // 	cout << endl;
  // 	for (Y = 0; Y < Yindex; Y++) {
  // 		for (X = 0; X < Xindex; X++) {
  // 			if (((Y * Xindex + X) % 2 != 0)) {
  // 				cout << (Y * Xindex + X) << " ";
  // 			}
  // 		}
  // 	}
  // 	cout << endl;
  // 	for (Y = 1; Y < Yindex; Y += 2) {
  // 		for (X = 1; X < Xindex; X += 2) {
  // 			cout << (Y * Xindex + X) << " ";
  // 		}
  // 	}
  // 	cout << endl;
  // 	exit(0);

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
      cout << c->getId() << endl;
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
    if (pPlugin == 0) {//normal parameter
      if (p->isSetValue()) {
        info->isResolved = true;
        for (Z = 0; Z < Zindex; Z++) {
          for (Y = 0; Y < Yindex; Y++) {
            for (X = 0; X < Xindex; X++) {
              info->value[Z * Yindex * Xindex + Y * Xindex + X] = p->getValue();
            }
          }
        }
      }
    } else {//spatial parameter plugin
      //???null?????
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
          INIT_DOUBLE(sInfo->dci[dc->getCoordinateIndex()]->value, numOfVolIndexes);
          for (Z = 0; Z < Zindex; Z++) {
            for (Y = 0; Y < Yindex; Y++) {
              for (X = 0; X < Xindex; X++) {
                sInfo->dci[dc->getCoordinateIndex()]->value[Z * Yindex * Xindex + Y * Xindex + X] = p->getValue();
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
          INIT_DOUBLE(sInfo->aci[ac->getCoordinateIndex()]->value, numOfVolIndexes);          

          for (Z = 0; Z < Zindex; Z++) {
            for (Y = 0; Y < Yindex; Y++) {
              for (X = 0; X < Xindex; X++) {
                sInfo->aci[ac->getCoordinateIndex()]->value[Z * Yindex * Xindex + Y * Xindex + X] = p->getValue();
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

  //time
  variableInfo *t_info = new variableInfo;
  InitializeVarInfo(t_info);
  varInfoList.push_back(t_info);
  t_info->id = (const char*)malloc(sizeof(char) * 1 + 1);
  strcpy(const_cast<char*>(t_info->id), "t");
  t_info->value = &sim_time;
  t_info->isResolved = true;

  //geometryDefinition
  for (i = 0; i < geometry->getNumGeometryDefinitions(); i++) {
    if (geometry->getGeometryDefinition(i)->isAnalyticGeometry()) {
      //AnalyticVolumes
      AnalyticGeometry *analyticGeo = static_cast<AnalyticGeometry*>(geometry->getGeometryDefinition(i));
      //gather information of compartment, domainType, analyticVolume
      for (j = 0; j < numOfCompartments; j++) {
        Compartment *c = loc->get(j);
        if (SpatialCompartmentPlugin *cPlugin = static_cast<SpatialCompartmentPlugin*>(c->getPlugin("spatial"))) {
          if (AnalyticVolume *analyticVol = analyticGeo->getAnalyticVolume(cPlugin->getCompartmentMapping()->getDomainType())) {
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
            avolInfo->isBoundary = new int[numOfVolIndexes];
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
            reversePolishInitial(avolInfo->rpInfo->varList, avolInfo->rpInfo->constList, avolInfo->rpInfo->opfuncList, avolInfo->isDomain, numOfASTNodes, Xindex, Yindex, Zindex, true);
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

  //merge external and internal analytic volumes and get boundary points of the geometry
  for (i = 0; i < geometry->getNumGeometryDefinitions(); i++) {
    AnalyticGeometry *analyticGeo = static_cast<AnalyticGeometry*>(geometry->getGeometryDefinition(i));
    if (geometry->getGeometryDefinition(i)->isAnalyticGeometry()) {
      for (j = 0; j < analyticGeo->getNumAnalyticVolumes(); j++) {
        AnalyticVolume *analyticVolEx = analyticGeo->getAnalyticVolume(j);
        analyticVolInfo *avolInfoEx = searchAvolInfoByDomainType(avolInfoList, analyticVolEx->getDomainType().c_str());
        for (k = 0; k < analyticGeo->getNumAnalyticVolumes(); k++) {
          AnalyticVolume *analyticVolIn = analyticGeo->getAnalyticVolume(k);
          if (analyticVolEx->getOrdinal() > analyticVolIn->getOrdinal()) {//ex's ordinal is larger than in's ordinal
            analyticVolInfo *avolInfoIn = searchAvolInfoByDomainType(avolInfoList, analyticVolIn->getDomainType().c_str());
            //merge
            for (Z = 0; Z < Zindex; Z += 2) {
              for (Y = 0; Y < Yindex; Y += 2) {
                for (X = 0; X < Xindex; X += 2) {
                  index = Z * Yindex * Xindex + Y * Xindex + X;
                  //if external == true && internal == true, then external = false;
                  if (static_cast<int>(avolInfoEx->isDomain[index]) == 1 && static_cast<int>(avolInfoIn->isDomain[index]) == 1) {
                    avolInfoEx->isDomain[index] = 0;
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

  //membrane compartment
  //adjacent??domain???????????????
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

  //????
  for (i = 0; i < avolInfoList.size(); i++) {
    analyticVolInfo *avolInfo = avolInfoList[i];
    if (avolInfo->isVol == false) {//avol is membrane
      analyticVolInfo *ad0 = 0, *ad1 = 0;
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
      avolInfo->isBoundary = new int[numOfVolIndexes];
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
        //??????
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
        //??????
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

  // 	vector<memArrayIndex*> memArrayIndexList = vector<memArrayIndex*>();
  // 	//info->value = new double[numOfVolIndexes];
  // 	vector<memArrayIndex*> mailZList = new vector<memArrayIndex*>[Zindex];
  // 	//membrane??????


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
        if (info->sp != 0 && (!info->sp->isSetConstant() || !info->sp->getConstant())) {
          info->delta = new double[4 * numOfVolIndexes];
        }
      }
      ast = const_cast<ASTNode*>((model->getInitialAssignment(info->id))->getMath());
    } else if (model->getRule(info->id) != 0 && model->getRule(info->id)->isAssignment()) {//assignment rule
      info->hasAssignmentRule = true;
      if (info->value == 0) {//value is not set yet
        info->value = new double[numOfVolIndexes];
        if (info->sp != 0 && (!info->sp->isSetConstant() || !info->sp->getConstant())) {
          //the species is variable
          info->delta = new double[4 * numOfVolIndexes];
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
  //????????????(????)
  while (!notOrderedInfo.empty()) {
    list<variableInfo*>::reverse_iterator it = notOrderedInfo.rbegin();
    while (it != notOrderedInfo.rend()) {
      variableInfo* current = (*it);
      if (isResolvedAll(current->dependence)) {
        reversePolishInitial(current->rpInfo->varList, current->rpInfo->constList, current->rpInfo->opfuncList, current->value, current->rpInfo->listNum, Xindex, Yindex, Zindex, current->inVol);
        current->isResolved = true;
        //orderedInfo.push_back(*it);
        notOrderedInfo.remove(current);
        //it++;
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
          } else {//membrane(??????)
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

  //simulation
#ifdef PRINT_DEBUG
  cout << endl;
#endif

  FILE *gp = popen((gnuplotExecutable + " -persist").c_str(),"w");
  clock_t re_time = 0, diff_time = 0, output_time = 0, ad_time = 0;
  clock_t sim_start = clock();
  for (t = 0; t <= static_cast<int>(end_time / dt); t++) {
    //cerr << static_cast<int>(100.0 * static_cast<double>(t) / (end_time / dt)) << endl;
    sim_time = t * dt;
    //output
    clock_t out_start = clock();
    if (count % static_cast<int>(out_time / dt) == 0) {
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
            fprintf(gp, "set cbrange[0.0:6.5]\n");
            fprintf(gp, "set palette defined (0 \"dark-blue\", 2 \"blue\", 4 \"green\", 8 \"yellow\", 10 \"red\")\n");
            fprintf(gp, "set xlabel \"x\"\n");
            fprintf(gp, "set ylabel \"y\"\n");
            fprintf(gp, "set title \"t=%lf\"\n", sim_time);
            fprintf(gp, "set terminal png truecolor\n");
#ifdef WIN32
            fprintf(gp, "set output \"/dev/null\"\n");
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
    }
    clock_t out_end = clock();
    output_time += out_end - out_start;
    count++;

    // 		stringstream ss;
    // 		ss << t;
    // 		string filename =  "./speed/" + ss.str() + ".txt";
    // 		ofstream ofs;
    // 		ofs.open(filename.c_str());
    //calculation
    //runge-kutta
    for (m = 0; m < 4; m++) {
      //reaction
      clock_t re_start = clock();
      for (i = 0; i < numOfReactions; i++) {
        Reaction *r = model->getReaction(i);
        //reactants
        for (j = 0; j < r->getNumReactants(); j++) {
          SpeciesReference *sr = r->getReactant(j);
          variableInfo *sInfo = searchInfoById(varInfoList, sr->getSpecies().c_str());
          Species *s = model->getSpecies(sr->getSpecies());
          if (!s->isSetConstant() || !s->getConstant()) {
            clock_t ast_start = clock();
            reversePolishRK(rInfoList[i], sInfo, Xindex, Yindex, Zindex, dt, m, reactants, sInfo->inVol, sInfo->dci);
            clock_t ast_end = clock();
            if (m == 0) {
              //ofs << "test" << endl;
              //cout << rInfoList[i]->rpInfo->listNum << " " << ((sim_end - sim_start) / static_cast<double>(CLOCKS_PER_SEC)) << endl;
              ofs << rInfoList[i]->rpInfo->listNum << " " << ((ast_end - ast_start) / static_cast<double>(CLOCKS_PER_SEC)) << endl;
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
      clock_t re_end = clock();
      re_time += (re_end - re_start);
      //rate rule
      for (i = 0; i < numOfRules; i++) {
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
      for (i = 0; i < numOfSpecies; i++) {
        variableInfo *sInfo = searchInfoById(varInfoList, los->get(i)->getId().c_str());
        double rk[4] = {0, 0.5, 0.5, 1.0};
        //diffusion
        clock_t diff_start = clock();
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
        clock_t diff_end = clock();
        diff_time += diff_end - diff_start;

        //advection
        clock_t ad_start = clock();
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

          //?????????
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

          //????????
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
          int index_Xm = 0, index_Ym = 0, /*index_Zm = 0,*/ index_XmYm = 0;
          double ux = 0.0, uy = 0.0,/* uz = 0.0,*/ Xsign = 0.0, Ysign = 0.0, /*Zsign = 0.0,*/ XX = 0.0, YY = 0.0/*, ZZ = 0.0*/;
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

        clock_t ad_end = clock();
        ad_time += ad_end - ad_start;
        //cout << "ad time: "<< ((ad_end - ad_start) / static_cast<double>(CLOCKS_PER_SEC)) << endl;
      }
    }//end of runge-kutta
    ofs.close();

    //update values
    for (i = 0; i < numOfSpecies; i++) {
      Species *s = los->get(i);
      variableInfo *sInfo = searchInfoById(varInfoList, s->getId().c_str());
      if (!s->isSetConstant() || !s->getConstant()) {
        for (Z = 0; Z < Zindex; Z++) {
          for (Y = 0; Y < Yindex; Y++) {
            for (X = 0; X < Xindex; X++) {
              if ((sInfo->inVol == true && X % 2 == 0 && Y % 2 == 0 && Z % 2 == 0) ||
                (sInfo->inVol == false && (Z * Yindex * Xindex + Y * Xindex + X) % 2 != 0)) {
                  index = Z * Yindex * Xindex + Y * Xindex + X;
                  //update values for the next time
                  if (static_cast<int>(sInfo->avi->isDomain[index] == 1) && sInfo->avi->isBoundary[index] == 0) {
                    sInfo->value[index] += dt * (sInfo->delta[index] + 2.0 * sInfo->delta[numOfVolIndexes + index] + 2.0 * sInfo->delta[2 * numOfVolIndexes + index] + sInfo->delta[3 * numOfVolIndexes + index]) / 6.0;
                    //reset delta
                    for (j = 0; j < 4; j++) {
                      sInfo->delta[j * numOfVolIndexes + index] = 0.0;
                    }
                  }
              }
            }
          }
        }
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
        }//boundary
        //at the edge of simulation area
        //parameter
        //Xp, Xm
        bcOfSpeciesInfo *bcos = searchBCOSInfoBySpeId(bcOfSpeciesInfoList, sInfo->id);
        if (dimension >= 1) {
          for (Z = 0; Z < Zindex; Z += 2) {
            for (Y = 0; Y < Yindex; Y += 2) {
              Xplus = Z * Yindex * Xindex + Y * Xindex + Xindex - 1;
              Xminus = Z * Yindex * Xindex + Y * Xindex;
              if (static_cast<int>(sInfo->avi->isDomain[Xplus]) == 1) {
                if (bcos->bcXpType == "Flux") {
                  sInfo->value[Xplus] = bcos->bcXp->value[Xplus] * deltaX + sInfo->value[Xplus - 2];
                } else if (bcos->bcXpType == "Value") {
                  sInfo->value[Xplus] = bcos->bcXp->value[Xplus];
                }
              }
              if (static_cast<int>(sInfo->avi->isDomain[Xminus]) == 1) {
                if (bcos->bcXmType == "Flux") {
                  sInfo->value[Xminus] = bcos->bcXm->value[Xminus] * deltaX + sInfo->value[Xminus + 2];
                } else if (bcos->bcXmType == "Value") {
                  sInfo->value[Xminus] = bcos->bcXm->value[Xminus];
                }
              }
            }
          }
        }
        if (dimension >= 2) {
          //Yp, Ym
          for (Z = 0; Z < Zindex; Z += 2) {
            for (X = 0; X < Xindex; X += 2) {
              Yplus = Z * Yindex * Xindex + (Yindex - 1) * Xindex + Xindex + X;
              Yminus = Z * Yindex * Xindex + X;
              if (static_cast<int>(sInfo->avi->isDomain[Yplus]) == 1) {
                if (bcos->bcYpType == "Flux") {
                  sInfo->value[Yplus] = bcos->bcYp->value[Yplus] * deltaY + sInfo->value[Yplus - 2];
                } else if (bcos->bcYpType == "Value") {
                  sInfo->value[Yplus] = bcos->bcYp->value[Yplus];
                }
              }
              if (static_cast<int>(sInfo->avi->isDomain[Yminus]) == 1 && sInfo->avi->isBoundary[Yminus] == 1) {
                if (bcos->bcYmType == "Flux") {
                  sInfo->value[Yminus] = bcos->bcYm->value[Yminus] * deltaY + sInfo->value[Yminus + 2];
                } else if (bcos->bcYmType == "Value") {
                  sInfo->value[Yminus] = bcos->bcYm->value[Yminus];
                }
              }
            }
          }
        }
        //Zp, Zm
        if (dimension == 3) {
          for (Y = 0; Y < Yindex; Y += 2) {
            for (X = 0; X < Xindex; X += 2) {
              Zplus = (Zindex - 1) * Yindex * Xindex + Y * Xindex + Xindex + X;
              Zminus = Y * Xindex + Xindex+ X;
              if (static_cast<int>(sInfo->avi->isDomain[Zplus]) == 1) {
                if (bcos->bcZpType == "Flux") {
                  sInfo->value[Zplus] = bcos->bcZp->value[Zplus] * deltaZ + sInfo->value[Zplus - 2];
                } else if (bcos->bcZpType == "Value") {
                  sInfo->value[Zplus] = bcos->bcZp->value[Zplus];
                }
              }
              if (static_cast<int>(sInfo->avi->isDomain[Zminus]) == 1) {
                if (bcos->bcZmType == "Flux") {
                  sInfo->value[Zminus] = bcos->bcZm->value[Zminus] * deltaZ + sInfo->value[Zminus + 2];
                } else if (bcos->bcZmType == "Value") {
                  sInfo->value[Zminus] = bcos->bcZm->value[Zminus];
                }
              }
            }
          }
        }
      }
    }
    //assignment rule
    for (i = 0; i < varInfoList.size(); i++) {
      variableInfo *info = varInfoList[i];
      if (info->hasAssignmentRule) {
        reversePolishInitial(info->rpInfo->varList, info->rpInfo->constList, info->rpInfo->opfuncList, info->value, info->rpInfo->listNum, Xindex, Yindex, Zindex, info->inVol);
        for (Z = 0; Z < Zindex; Z += 2) {
          for (Y = 0; Y < Yindex; Y += 2) {
            for (X = 0; X < Xindex; X += 2) {
              index = Z * Yindex * Xindex + Y * Xindex + X;
              if (static_cast<int>(info->avi->isDomain[index]) == 0) {
                info->value[index] = 0.0;
              }
            }
          }
        }
      }
    }
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
  cerr << "output_time: "<< (output_time / static_cast<double>(CLOCKS_PER_SEC)) << endl;
  pclose(gp);
  //cerr << 100 << endl;
  //free
  freeVarInfo(varInfoList, t_info->value);
  freeBcOfSpeciesInfo(bcOfSpeciesInfoList);
  freeAvolInfo(avolInfoList);
  freeRInfo(rInfoList);
}
