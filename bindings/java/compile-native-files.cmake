###############################################################################
#
# Description       : CMake build script for native java files
# Original author(s): Frank Bergmann <fbergman@caltech.edu>
#
###############################################################################

message("Creating: spatialSBML.jar")

# find all sources
file(GLOB_RECURSE SOURCE_FILES RELATIVE ${BIN_DIRECTORY} ${BIN_DIRECTORY}/java-files/org/spatial/*.java)

# convert paths
set(NATIVE_FILES)
foreach(javaFile ${SOURCE_FILES})
	file(TO_NATIVE_PATH ${javaFile} temp)
	set(NATIVE_FILES ${NATIVE_FILES} ${temp})
endforeach()

# delete file if it exists
if (EXISTS ${BIN_DIRECTORY}/spatialSBML.jar)
	file(REMOVE ${BIN_DIRECTORY}/spatialSBML.jar)	
endif()

# compile files
execute_process(
	COMMAND "${Java_JAVAC_EXECUTABLE}"
		 -source 1.5
		 -target 1.5
		 -d java-files
		 ${NATIVE_FILES}	
	WORKING_DIRECTORY "${BIN_DIRECTORY}"
)

# enumerate class files
file(GLOB_RECURSE CLASS_FILES RELATIVE ${BIN_DIRECTORY}/java-files ${BIN_DIRECTORY}/java-files/org/spatial/*.class)
set(NATIVE_CLASS_FILES)
foreach(classFile ${CLASS_FILES})
	file(TO_NATIVE_PATH ${classFile} temp)
	set(NATIVE_CLASS_FILES ${NATIVE_CLASS_FILES} ${temp})
endforeach()

# create jar
execute_process(
	COMMAND "${Java_JAR_EXECUTABLE}"
		 -cvfm ..${PATH_SEP}spatialSBML.jar
		 ${SRC_DIRECTORY}/Manifest.txt
		 ${NATIVE_CLASS_FILES}	
	WORKING_DIRECTORY "${BIN_DIRECTORY}/java-files"
)
