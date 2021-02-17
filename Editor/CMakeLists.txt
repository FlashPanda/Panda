# get all header files and source files
FILE(GLOB EDITOR_SOURCES ./*.*)

IF (${WIN32})
	add_executable(Editor
		${EDITOR_SOURCES})
		
	target_link_libraries(Editor Panda)
	
ELSE (${WIN32})
	add_executable(Editor
		${EDITOR_SOURCES})
		
	target_link_libraries(Editor
		Panda
		)
ENDIF(${WIN32})
SOURCE_GROUP(TREE ${PANDA_ROOT} FILES ${EDITOR_SOURCES})
set_property(TARGET  Editor PROPERTY FOLDER "Editor")


