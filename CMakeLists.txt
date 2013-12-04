####################################################################
#
# CMake Build Script for SBML Spatial, a UI for the 
# spatial simulator from Akira
#
# Frank T. Bergmann
#

cmake_minimum_required(VERSION 2.6)
project(spatial)

set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "2")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_SOURCE_GENERATOR "ZIP")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/license.txt")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/readme.md")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_IGNORE_FILES
  "aa;aa.zip;screenshots;.DS_STORE;testmodels/.DS_Store;/build/;/.bzr/;~$;${CPACK_SOURCE_IGNORE_FILES}")
set(CPACK_NSIS_CREATE_ICONS "CreateShortCut '\$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\SBML Spatial.lnk' '\$INSTDIR\\\\bin\\\\SpatialUI.exe'")
include(CPack)


find_package(OpenMP)
find_package(ZLIB)

find_library(LIBSBML_LIBRARY
  NAMES sbml-static.lib sbml-static sbml
  PATHS /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64
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


###############################################################################
#
# Need some variables set up, such as the name for the libSBML
# library and the Path and file separator characters. The
# MISC_PREFIX variable will cause libsbml bindings, examples and
# documentation to be installed in PREFIX/${MISC_PREFIX}.
#

if(UNIX OR CYGWIN)
    set(PATH_SEP "/")
    set(FILE_SEP ":")
else()
    set(PATH_SEP "\\")
    set(FILE_SEP ";")
endif()

set(EXTRA_LIBS "" CACHE STRING "List of additional libraries to link against." )
set(EXTRA_FLAGS "" CACHE STRING "List of extra flags." )
set(EXTRA_INCLUDE "" CACHE STRING "List of include diretories." )

if (EXTRA_FLAGS)
foreach(${flag} ${EXTRA_FLAGS})
	add_definitions(-D${flag})
endforeach()
endif()

###############################################################################
#
# SBW Integration
#
option(ENABLE_SBW   "Enable SBW Integration."     ON)
if (ENABLE_SBW)

find_library(SBW_LIBRARY
  NAMES SBW-static.lib SBW-static SBW
  PATHS /usr/lib /usr/local/lib /usr/lib64 /usr/local/lib64
        ${CMAKE_SOURCE_DIR}
        ${LIBSBML_DEPENDENCY_DIR}/lib
  DOC "The file name of the SBW library."
)

find_path(SBW_INCLUDE_DIR
  NAMES SBW/SBW.h
  PATHS /usr/include /usr/local/include
    ${CMAKE_SOURCE_DIR}/include
    ${LIBSBML_DEPENDENCY_DIR}/include
  DOC "The SBW include directory."
)

endif()


option(WITH_CSHARP   "Generate C# language bindings."     OFF)
option(WITH_JAVA     "Generate Java language bindings."   OFF)


###############################################################################
#
# Find the C# compiler to use and set name for resulting library
#

if(WITH_CSHARP)
    find_program(CSHARP_COMPILER
        NAMES gmcs csc
        PATHS C:/Windows/Microsoft.NET/Framework/v2.0.50727/ /usr/bin /usr/local/bin
        DOC "The file name of the C# compiler."
        )
    if(UNIX)
    else()
        if(CMAKE_SIZEOF_VOID_P EQUAL 4)
            # mark libsbml library as x86
            set(CSHARP_EXTRA_ARGS -platform:x86 )
        elseif(CMAKE_SIZEOF_VOID_P EQUAL 8)
            # mark libsbml library as x64
            set(CSHARP_EXTRA_ARGS -platform:x64 )
        endif()
    endif()

endif(WITH_CSHARP)



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

option(ENABLE_OPENMP "Enable OpenMP if available" ON)
if(OPENMP_FOUND AND ENABLE_OPENMP)
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

if(WITH_CSHARP)
    message(STATUS "  Using C#                      = ${CSHARP_COMPILER}")
endif()

if(WITH_JAVA)
    message(STATUS "  Using Java                    = ${Java_JAVA_EXECUTABLE}")
endif()

add_subdirectory(bindings)


