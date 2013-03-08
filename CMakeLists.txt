####################################################################
#
# CMake Build Script for spatial simulator from akira
#
# $Author$
# $Id$
# $HeadURL$
#

cmake_minimum_required(VERSION 2.6)
project(spatial)

set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_SOURCE_GENERATOR "ZIP")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_IGNORE_FILES
  "aa;aa.zip;screenshots;.DS_STORE;testmodels/.DS_Store;/build/;/.bzr/;~$;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)


find_package(OpenMP)
find_package(ZLIB)

find_library(LIBSBML_LIBRARY
  NAMES sbml-static.lib sbml-static sbml
  PATHS /usr/lib /usr/local/lib
        ${CMAKE_SOURCE_DIR}
        ${LIBSBML_DEPENDENCY_DIR}/lib
  DOC "The file name of the libsbml library."
)

find_path(LIBSBML_INCLUDE_DIR
  NAMES sbml/SBase.h
  PATHS /usr/include /usr/local/include
    ${CMAKE_SOURCE_DIR}/include
    ${LIBSBML_DEPENDENCY_DIR}/include
  DOC "The libsbml include directory."
)


set(EXTRA_LIBS "" CACHE STRING "List of additional libraries to link against." )

if(WIN32 AND NOT CYGWIN)
 add_definitions(-DLIBSBML_EXPORTS -DWIN32)
 if(MSVC)
        add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_DEPRECATE)
        set(BUILD_DEFINITIONS "${BUILD_DEFINITIONS} -D_CRT_SECURE_NO_WARNINGS")
        option(WITH_STATIC_RUNTIME "Compile using the static MSVC Runtime." OFF)
        if(WITH_STATIC_RUNTIME)
            foreach(flag_var
                CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
                CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO
                CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
                CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)

                if(${flag_var} MATCHES "/MD")
                    string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
                endif(${flag_var} MATCHES "/MD")
            endforeach(flag_var)
        add_definitions( -D_MT)
        endif(WITH_STATIC_RUNTIME)
  endif()
endif()

if (ZLIB_FOUND)
include_directories(${ZLIB_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${ZLIB_LIBRARIES})
endif(ZLIB_FOUND)

if(OPENMP_FOUND)
 set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
 set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(SpatialSBML)

option(BUILD_CL "Build Command line executable" ON)
if (BUILD_CL)
  add_subdirectory(SpatialCL)
endif()


option(BUILD_UI "Build Graphical User Interface (Qt based)" OFF)
if (BUILD_UI)
  add_subdirectory(SpatialUI)
endif()


file(GLOB EXAMPLE_FILES examples/*.xml )
install(FILES ${EXAMPLE_FILES} DESTINATION examples)


message(STATUS         
"""
 Project Spatial SBML ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}
 ==========================
 
 Install Prefix : ${CMAKE_INSTALL_PREFIX}
                
 SBML INCLUDE   : ${LIBSBML_INCLUDE_DIR}
 SBML Library   : ${LIBSBML_LIBRARY}
                
 C   FLAGS      : ${CMAKE_C_FLAGS}
 CXX FLAGS      : ${CMAKE_CXX_FLAGS}
 Extra Libs     : ${EXTRA_LIBS}
 
 
"""
)

