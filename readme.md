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

![SpatialUI, the interface](http://)

From here, you could simply open one of the example models supplied. We use a custom annotation to store settings such as the last selected species, the colors used for them as well as the bounds. Of course you can open any model in the SBML Spatial format. 

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
