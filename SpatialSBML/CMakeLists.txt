####################################################################
#
# CMake Build Script for the Command line interface to the spatial simulator from akira
#
# Frank T. Bergmann (fbergman@caltech.edu)
#
project (SpatialSBML)

file(GLOB CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp )
file(GLOB   H_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h )

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${LIBSBML_INCLUDE_DIR})

add_library(SpatialSBML STATIC ${CPP_FILES} )
target_link_libraries(SpatialSBML ${LIBSBML_LIBRARY} ${EXTRA_LIBS})

