####################################################################
#
# CMake Build Script for the UI for spatial simulator from akira
#
# Frank T. Bergmann (fbergman@caltech.edu)
#
project (SpatialUI)

find_package(Qt4 4.7.0 REQUIRED QtCore QtOpenGL QtGui Qt3Support)

file(GLOB CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp )
file(GLOB   H_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h )
file(GLOB  UI_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.ui )
file(GLOB RCC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.qrc )

qt4_wrap_cpp(HEADERS_MOC ${H_FILES})
qt4_wrap_ui(FORMS_HEADERS ${UI_FILES})
qt4_add_resources(RESOURCES_RCC ${RCC_FILES})

include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR})
include_directories(../SpatialSBML)
include_directories(${LIBSBML_INCLUDE_DIR})

set(EXECUTABLE_FLAGS)

if(WIN32 AND NOT CYGWIN)
  option(ENABLE_UI_CONSOLE "Show the console window (to see debug messages) when running Spatial UI" OFF)
  if (NOT ENABLE_UI_CONSOLE)
    set(EXECUTABLE_FLAGS WIN32)
    if(MSVC)
      set(CPP_FILES ${CPP_FILES} spatialUI.rc)
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ENTRY:mainCRTStartup ")
    endif(MSVC)
  endif(NOT ENABLE_UI_CONSOLE)
endif(WIN32 AND NOT CYGWIN)

if (APPLE)
 
set(EXECUTABLE_FLAGS MACOSX_BUNDLE)
set(MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/SpatialUI.plist )
set(MACOSX_BUNDLE_BUNDLE_VERSION  "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(MACOSX_BUNDLE_INFO_STRING "SpatialUI ${MACOSX_BUNDLE_BUNDLE_VERSION} is distributed under the BSD." )
set(MACOSX_BUNDLE_GUI_IDENTIFIER  "org.spatial.SpatialUI" )
set(MACOSX_BUNDLE_BUNDLE_NAME "SpatialUI")
set(MACOSX_BUNDLE_LONG_VERSION_STRING "${MACOSX_BUNDLE_BUNDLE_VERSION}")
set(MACOSX_BUNDLE_SHORT_VERSION_STRING  "${MACOSX_BUNDLE_BUNDLE_VERSION}")
set(MACOSX_BUNDLE_COPYRIGHT ${MACOSX_BUNDLE_INFO_STRING})

endif(APPLE)

add_executable(SpatialUI ${EXECUTABLE_FLAGS} ${CPP_FILES}
	${H_FILES}
    ${HEADERS_MOC}
    ${FORMS_HEADERS}
    ${RESOURCES_RCC}
	${SPATIAL_CPP_FILES}
	)
target_link_libraries(SpatialUI ${EXTRA_LIBS} SpatialSBML ${QT_LIBRARIES})


file(GLOB PALETTE_FILES palettes/*.txt )


SET(DIRS)
foreach(${lib} ${EXTRA_LIBS})
  if (EXISTS ${lib})
	get_filename_component(dir ${lib} PATH)	
	SET(DIRS  ${DIRS} ${dir})
  endif()
endforeach()
if (NOT DIRS)
   if (EXISTS "${EXTRA_LIBS}")
	get_filename_component(DIRS ${EXTRA_LIBS} PATH)
   endif()
endif()
if (QT_BINARY_DIR)
SET(DIRS ${QT_BINARY_DIR} ${DIRS})
endif()


if(NOT APPLE)

SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/SpatialUI${CMAKE_EXECUTABLE_SUFFIX}" CACHE INTERNAL "")
install(FILES ${PALETTE_FILES} DESTINATION bin)
install(TARGETS SpatialUI DESTINATION bin)
# ensure that next to the binaries also dependencies are copied
install(CODE "
    include(BundleUtilities)
    fixup_bundle(\"${APPS}\" \"\" \"${DIRS}\")
    " COMPONENT Runtime)

else()

install(FILES ${PALETTE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/SpatialUI.app/Contents/MacOS)

# on OS X everything is different, here we need to first fix the app bundle,
# and only then install it
SET(APPS "${CMAKE_CURRENT_BINARY_DIR}/SpatialUI.app" CACHE INTERNAL "")
# ensure that next to the binaries also dependencies are copied
install(CODE "
    include(BundleUtilities)
    fixup_bundle(\"${APPS}\" \"\" \"${DIRS}\")
    " COMPONENT Runtime)
install(TARGETS SpatialUI DESTINATION bin)
endif()

include(InstallRequiredSystemLibraries)
