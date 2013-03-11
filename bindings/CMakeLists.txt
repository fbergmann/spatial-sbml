###############################################################################
#
# Description       : CMake build script for spatialSBML language bindings
# Original author(s): Frank Bergmann <fbergman@caltech.edu>
# Organization      : California Institute of Technology
#
###############################################################################

find_program(SWIG_EXECUTABLE
    NAMES swig
    PATHS
       c:/swigwin-2.0.8
       c:/swigwin-2.0.7
       c:/swigwin-2.0.4
       c:/swigwin-2.0.2
       c:/swigwin-2.0.1
       c:/swigwin-2.0.0
       /usr/local/bin
       /opt/local/bin
       /usr/bin
    DOC "The file name of the swig executable."
    )
set(SWIG_EXTRA_ARGS)


#
# and build selected language bindings
#

if(WITH_CSHARP)
add_subdirectory(csharp)
endif()

if(WITH_JAVA)
add_subdirectory(java)
endif()
