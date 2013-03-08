####################################################################
#
# CMake Build Script for the Command line interface to the spatial simulator from akira
#
# Frank T. Bergmann (fbergman@caltech.edu)
#
project (SpatialCL)

file(GLOB CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp )
file(GLOB   H_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h )

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR})
include_directories(../SpatialSBML)
include_directories(${LIBSBML_INCLUDE_DIR})


set(EXECUTABLE_FLAGS)

add_executable(SpatialCL main.cpp )
target_link_libraries(SpatialCL SpatialSBML)


install(TARGETS SpatialCL DESTINATION bin)