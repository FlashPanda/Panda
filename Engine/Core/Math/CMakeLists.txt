# get all header files and source files
FILE(GLOB MATH_HEADERS ./*.hpp)
FILE(GLOB MATH_SOURCES ./*.cpp)

# add to math lib
add_library(PandaMath
	${MATH_SOURCES}
	${MATH_HEADERS})