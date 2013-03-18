# Spatial SBML
This project contains a library, as well as a front end application for spatial simulation of SBML Models (using [SBML Level 3](http://sbml.org/Documents/Specifications) with the [Spatial package](http://sbml.org/Community/Wiki/SBML_Level_3_Proposals/Spatial_Geometries_and_Spatial_Processes "Spatial package") ). This project is made possible in cooperation with the [Funahashi Lab](http://fun.bio.keio.ac.jp/) who developed the simulator. (The initial presentation from COMBINE 2011 can still be found online [here](http://co.mbine.org/events/COMBINE_2011/agenda?q=system/files/2011-09-04-combine-funahashi-spatial-simulator.pdf).) 

## How to build
This project uses CMake, as such building ought to be as easy as running:

```
mkdir build  
cd build  
cmake -DBUILD_UI=ON <directory with sources>  
make  
make install
```

## Running the GUI
Once installed, simply double click SpatialUI, and the following interface will come up.   

![SpatialUI, the interface](https://raw.github.com/fbergmann/spatial-sbml/master/screenshots/2012-06-12_fish300x300.png)

From here, you could simply open one of the example models supplied. We use a custom annotation to store settings such as the last selected species, the colors used for them as well as the bounds. Of course you can open any model in the SBML Spatial format and define these settings yourself: 

- **Selecting the variables to plot**: from the *Palette Assignment* panel, simply select the Species to display, choose one of the palettes distributed with the program, and define the maximum concentration (will be used to scale the palette). Next click 'add'. While it is possible to add multiple species that way I recommend to choose < 3 (and those in the palettes 'black-red', 'black-green' and 'black-blue' as that maps best to the RGB color scheme. 
- **Displaying Concentrations at a certain point**: Simply click on a coordinate, and concentrations for all selected species will be displayed. 
- **Modifying the concentrations at a coordinate**: First select the species whose concentrations you would like to change, next specify the concentration and enable the 'Apply' button. Then simply click / drag the mouse along the canvas to change the concentrations at those points. 

Should you have further questions, or concerns please let me know. 
## SpatialSBML API 
I'm currently in the process of exposing the API through SWIG. A first prototype is available in the examples folder. To test it, compile with with the CMake option `WITH_JAVA=ON`, and then modify the example in `examples/bindings` to adjust the location of the sample file.  

The API consists of these functions: 

* SpatialSimulator ()
* virtual ~SpatialSimulator ()
* void initFromFile (const std::string &fileName, int xdim, int ydim, int zdim=1)
* void initFromString (const std::string &sbml, int xdim, int ydim, int zdim=1)
* void run (double endTime, double step)
* double oneStep (double initialTime, double step)
* void setParameterUniformly (const std::string &id, double value)
* void setParameter (const std::string &id, double value)
* int getXDim () const 
* int getYDim () const 
* int getZDim () const 
* double getVariableAt (const std::string &id, int x, int y, int z=0)

Currently I'm also exposing these functions that I'm not happy about as they reveal too much about the underlying implementation. In addition they are not wrapped by the language bindings: 

* SpatialSimulator (SBMLDocument *doc, int xdim, int ydim, int zdim=1)
* void initFromModel (SBMLDocument *doc, int xdim, int ydim, int zdim=1)
* const Model * getModel () const 
* int getIndexForPosition (double x, double y)
* int getVariableLength () const 
* int getGeometryLength () const 
* void setParameterUniformly (variableInfo *info, double value)
* double * getVariable (const std::string &id, int &length)
* double * getGeometry (const std::string &id, int &length)
* boundaryType * getBoundaryType (const std::string &id, int &length)
* int * getBoundary (const std::string &id, int &length)
* double * getX (int &length)
* double * getY (int &length)
* double * getZ (int &length)
* void setGnuplotExecutable (const std::string &location)
* const std::string & getGnuplotExecutable () const 

## Libraries
This project requires the following libraries: 

- [libSBML](http://sbml.org/Software/libSBML) along with a supported xml library (xml2, expat, xerces-c)
- [Qt4](http://qt-project.org/) for the frontend

## License
This project is licensed under the BSD license: 

```
Copyright (c) 2013, Frank T. Bergmann  
All rights reserved. 

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met: 

Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer. Redistributions in 
binary form must reproduce the above copyright notice, this list of 
conditions and the following disclaimer in the documentation and/or 
other materials provided with the distribution. THIS SOFTWARE IS 
PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

```
