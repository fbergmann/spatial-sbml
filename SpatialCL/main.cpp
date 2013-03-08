//============================================================================
// Name        : SBMLSimulator.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================


#include <ctime>
#include <iostream>
#include <sbml/SBMLTypes.h>

#define PRINT_DEBUG

#include "spatialsimulator.h"

using namespace std;
LIBSBML_CPP_NAMESPACE_USE;



int main(int argc, const char *argv[])
{
  clock_t start = clock();

  if (argc < 2)
  {
    cerr << "Need one argument: sbml file with spatial simulation" << endl;
    exit(1);
  }
  string fileName (argv[1]);
  SBMLDocument *doc = readSBML(fileName.c_str());

  if (doc->getModel() == NULL)
  {
    cerr << "encountered fatal errors while reading the file." << endl;
    doc->printErrors();
    exit(1);
  }

  if (argc == 3)
  {
    return SpatialSimulator::runOldMain(argc, argv);
  }

  cout << "SpatialPlugin: " << (doc->getModel()->getPlugin("spatial") != NULL ? "not null": "null") <<endl;
  cout << "spatial required: "<< (doc->getPkgRequired("spatial") ? "yes": "no") << endl;
  cout << "req required: "<< (doc->getPkgRequired("req") ? "yes": "no") << endl;

  if (doc->getModel()->getPlugin("spatial") != NULL && 
    doc->getPkgRequired("spatial") && 
    doc->getPkgRequired("req")) {
      // probably a spatial simulation 
      SpatialSimulator sim(doc, 101, 101);
      sim.setGnuplotExecutable("C:/gnuplot/bin/gnuplot.exe");
      sim.run(1001, 0.01);
      //spatialSimulator(doc, argc, argv);
  } 
  else {
    cerr << "Please provide a spatial sbml file." << endl; 
  }
  delete doc;

  clock_t end = clock();

  cerr << "time: " << ((end - start) / static_cast<double>(CLOCKS_PER_SEC)) << endl;

  return 0;
}

