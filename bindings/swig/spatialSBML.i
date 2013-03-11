/**
 * \file    spatialSBML.i
 * \brief   Language-independent SWIG directives for wrapping spatialSBML
 * \author  Frank T. Bergmann
 * 
 */

%module spatialSBML

//#pragma SWIG nowarn=473,401,844

%{
#include "spatialSBML.h"
	
#include "local.cpp"
%}

/**
 *
 * Includes a language specific interface file.
 *
 */

%include local.i

/**
 * Ignore operator= and operator<< on all SBML objects.
 */
%ignore *::operator=;
%ignore *::operator<<;
%ignore operator==;
%ignore operator!=;


%include "std_string.i"

%include  SpatialSBML/spatialsimulator.h




