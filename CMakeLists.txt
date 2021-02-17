# This is for annotation
# The latest docs website is :https://cmake.org/cmake/help/latest/index.html
cmake_minimum_required (VERSION 3.2)	# Minimum limited version, 3.1

project(Panda) # project name

# initialize directories 
set (CMAKE_C_STNDARD 17)
set (CMAKE_CXX_STANDARD 17)
set (PANDA_ROOT ${PROJECT_SOURCE_DIR})
set (PANDA_ENGINE_SOURCE_ROOT ${PANDA_ROOT}/Engine)
set (PANDA_TEST_SOURCE_ROOT ${PANDA_ROOT}/Test)
set (PANDA_EDITOR_SOURCE_ROOT ${PANDA_ROOT}/Editor)
set (PANDA_BUILD_ROOT ${PANDA_ROOT}/Build)
set (PANDA_EXTERNAL_ROOT ${PANDA_ROOT}/3rdParty)
set (PANDA_EXTERNAL_INCLUDE ${PANDA_EXTERNAL_ROOT}/include)
set (PANDA_EXTERNAL_DEBUG_LIB ${PANDA_EXTERNAL_ROOT}/binaries/debug)
set (PANDA_EXTERNAL_RELEASE_LIB	${PANDA_EXTERNAL_ROOT}/binaries/release)

# use solution folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)  

IF(${UNIX})
    set(PANDA_TARGET_PLATFORM "Linux")
ELSEIF(${WIN32})
    set(PANDA_TARGET_PLATFORM "Windows")
ENDIF(${UNIX})

include_directories(${PANDA_ROOT})
include_directories(${PANDA_ENGINE_SOURCE_ROOT})
include_directories(${PANDA_EXTERNAL_INCLUDE})
include_directories(${PANDA_EXTERNAL_ROOT})

add_subdirectory(Engine)		# sub directory, must have CMakeLists.txt in it
add_subdirectory(Test)
add_subdirectory(Editor)

# ------------------------- Begin Generic CMake Variable Logging ------------------
# if you are building in-source, this is the same as CMAKE_SOURCE_DIR, otherwise
# this is the top level directory of your build tree
MESSAGE(STATUS "CMAKE_BINARY_DIR:			" ${CMAKE_BINARY_DIR})

# if you are building in-source, this is the same as CMAKE_CURRENT_SOURCE_DIR, otherwise this
# is the directory where the compiled or generated files from the current CMakeLists.txt will go to.
MESSAGE(STATUS "CMAKE_CURRENT_BINARY_DIR:	" ${CMAKE_CURRENT_BINARY_DIR})

# this is the directory, from which cmake was started, i.e. the top level source directory
MESSAGE(STATUS "CMAKE_SOURCE_DIR:			" ${CMAKE_SOURCE_DIR})

# this is the directory where the currently processed CMakeLists.txt is located in.
MESSAGE(STATUS "CMAKE_CURRENT_SOURCE_DIR:	" ${CMAKE_CURRENT_SOURCE_DIR})

# contains the full path to the top level directory of your build tree
MESSAGE(STATUS "PROJECT_BINARY_DIR:			" ${PROJECT_BINARY_DIR})

# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command
MESSAGE(STATUS "PROJECT_SOURCE_DIR:			" ${PROJECT_SOURCE_DIR})

# set this variable to specify a common place where CMake should put all executable files
# (instead of CMAKE_CURRENT_BINARY_DIR)
MESSAGE(STATUS "EXECUTABLE_OUTPUT_PATH: 	" ${EXECUTABLE_OUTPUT_PATH} )

# set this variable to specify a common place where CMake should put all executable files
# (instead of CMAKE_CURRENT_BINARY_DIR)
MESSAGE(STATUS "LIBRARY_OUTPUT_PATH:		" ${LIBRARY_OUTPUT_PATH})

# tell CMake to search first in directories listed in CMAKE_MODULE_PATH
# when you use FIND_PACKAGE() or INCLUDE()
MESSAGE(STATUS "CMAKE_MODULE_PATH:			" ${CMAKE_MODULE_PATH})

# this is the complete path of the cmake which runs currently (e.g. /usr/local/bin/cmake) 
MESSAGE(STATUS "CMAKE_COMMAND:				" ${CMAKE_COMMAND})

# this is the CMake installation directory
MESSAGE(STATUS "CMAKE_ROOT:					" ${CMAKE_ROOT})

# this is the filename including the complete path of the file where this variable is used. 
MESSAGE(STATUS "CMAKE_CURRENT_LIST_FILE:	" ${CMAKE_CURRENT_LIST_FILE})

# this is linenumber where the variable is used
MESSAGE(STATUS "CMAKE_CURRENT_LIST_LINE:	" ${CMAKE_CURRENT_LIST_LINE})

# this is used when searching for include files e.g. using the FIND_PATH() command.
MESSAGE(STATUS "CMAKE_INCLUDE_PATH:			" ${CMAKE_INCLUDE_PATH})

# this is used when searching for libraries e.g. using the FIND_LIBRARY() command.
MESSAGE(STATUS "CMAKE_LIBRARY_PATH:			" ${CMAKE_LIBRARY_PATH})

## HOST
# the complete system name, e.g. "Linux-2.4.22", "FreeBSD-5.4-RELEASE" or "Window 5.1"
MESSAGE(STATUS "CMAKE_HOST_SYSTEM:			" ${CMAKE_HOST_SYSTEM})

# the short system name, e.g. "Linux", "FreeBSD" or "Windows"
MESSAGE(STATUS "CMAKE_HOST_SYSTEM_NAME:		" ${CMAKE_HOST_SYSTEM_NAME})

# only the version part of CMAKE_SYSTEM
MESSAGE(STATUS "CMAKE_HOST_SYSTEM_VERSION:	" ${CMAKE_HOST_SYSTEM_VERSION})

# the process name (e.g. "Intel(R) Pentium(R) M processor 2.00 GHz")
MESSAGE(STATUS "CMAKE_HOST_SYSTEM_PROCESSOR: " ${CMAKE_HOST_SYSTEM_PROCESSOR})

## TARGET
# the complete system name, e.g. "Linux-2.4.22", "FreeBSD-5.4-RELEASE" or "Windows 5.1" 
MESSAGE(STATUS "CMAKE_SYSTEM:				" ${CMAKE_SYSTEM})

# the short system name, e.g. "Linux", "FreeBSD" or "Windows"
MESSAGE(STATUS "CMAKE_SYSTEM_NAME:			" ${CMAKE_SYSTEM_NAME})

# only the version part of CMAKE_SYSTEM
MESSAGE(STATUS "CMAKE_SYSTEM_VERSION:		" ${CMAKE_SYSTEM_VERSION})

# the processor name (e.g. "Intel(R) Pentium(R) M processor 2.00GHz") 
MESSAGE(STATUS "CMAKE_SYSTEM_PROCESSOR:		" ${CMAKE_SYSTEM_PROCESSOR})

# is TRUE on all UNIX-like OS's, including Apple OS X and CygWin
MESSAGE(STATUS "UNIX:						" ${UNIX})

# is TRUE on Windows, including CygWin 
MESSAGE(STATUS "WIN32:						" ${WIN32})

# is TRUE on Apple OS X
MESSAGE(STATUS "APPLE:						" ${APPLE})

# is TRUE when using the MinGW compiler in Windows
MESSAGE(STATUS "MINGW:						" ${MINGW})

# is TRUE on Windows when using the CygWin version of cmake
MESSAGE(STATUS "CYGWIN:						" ${CYGWIN})

# is TRUE on Windows when using a Borland compiler 
MESSAGE(STATUS "BROLAND:					" ${BORLAND})

# Microsoft compiler 
MESSAGE(STATUS "MSVC:						" ${MSVC})
MESSAGE(STATUS "MSVC_IDE:					" ${MSVC_IDE})
MESSAGE(STATUS "MSVC60:						" ${MSVC60})
MESSAGE(STATUS "MSVC70:						" ${MSVC70})
MESSAGE(STATUS "MSVC71:						" ${MSVC71})
MESSAGE(STATUS "MSVC80:						" ${MSVC80})
MESSAGE(STATUS "CMAKE_COMPILER_2005:		" ${CMAKE_COMPILER_2005})

# set this to true if you don't want to rebuild the object files if the rules have changed, 
# but not the actual source files or headers (e.g. if you changed the some compiler switches) 
MESSAGE(STATUS "CMAKE_SKIP_RULE_DEPENDENCY:	" ${CMAKE_SKIP_RULE_DEPENDENCY})

# since CMake 2.1 the install rule depends on all, i.e. everything will be built before installing. 
# If you don't like this, set this one to true.
MESSAGE( STATUS "CMAKE_SKIP_INSTALL_ALL_DEPENDENCY: " ${CMAKE_SKIP_INSTALL_ALL_DEPENDENCY} )

# If set, runtime paths are not added when using shared libraries. Default it is set to OFF
MESSAGE(STATUS "CMAKE_SKIP_RPATH:			" ${CMAKE_SKIP_RPATH})

# set this to true if you are using makefiles and want to see the full compile and link 
# commands instead of only the shortened ones
MESSAGE(STATUS "CMAKE_VERBOSE_MAKEFILE:		" ${CMAKE_VERBOSE_MAKEFILE})

# this will cause CMake to not put in the rules that re-run CMake. This might be useful if 
# you want to use the generated build files on another machine. 
MESSAGE( STATUS "CMAKE_SUPPRESS_REGENERATION: " ${CMAKE_SUPPRESS_REGENERATION} )

# A simple way to get switches to the compiler is to use ADD_DEFINITIONS(). 
# But there are also two variables exactly for this purpose: 

# the compiler flags for compiling C sources 
MESSAGE(STATUS "CMAKE_C_FLAGS:				" ${CMAKE_C_FLAGS})
MESSAGE(STATUS "CMAKE_C_FLAGS_DEBUG:		" ${CMAKE_C_FLAGS_DEBUG})
MESSAGE(STATUS "CMAKE_C_FLAGS_RELEASE:		" ${CMAKE_C_FLAGS_RELEASE})
MESSAGE(STATUS "CMAKE_C_FLAGS_RELWITHEBINFO: " ${CMAKE_C_FLAGS_RELWITHDEBINFO})
MESSAGE(STATUS "CMAKE_C_FLAGS_MINSIZEREL:	" ${CMAKE_C_FLAGS_MINSIZEREL})

# the compiler flags for compiling C++ sources 
MESSAGE(STATUS "CMAKE_CXX_FLAGS:			" ${CMAKE_CXX_FLAGS})
MESSAGE(STATUS "CMAKE_CXX_FLAGS_DEBUG:		" ${CMAKE_CXX_FLAGS_DEBUG})
MESSAGE(STATUS "CMAKE_CXX_FLAGS_RELEASE:	" ${CMAKE_CXX_FLAGS_RELEASE})
MESSAGE(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO: " ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
MESSAGE(STATUS "CMAKE_CXX_FLAGS_MINSIZEREL: " ${CMAKE_CXX_FLAGS_MINSIZEREL})

# Choose the type of build.  Example: SET(CMAKE_BUILD_TYPE Debug) 
MESSAGE(STATUS "CMAKE_BUILD_TYPE:			" ${CMAKE_BUILD_TYPE})

# if this is set to ON, then all libraries are built as shared libraries by default.
MESSAGE(STATUS "BUILD_SHARED_LIBS:			" ${BUILD_SHARED_LIBS})

# the compiler used for C files 
MESSAGE(STATUS "CMAKE_C_COMPILER:			" ${CMAKE_C_COMPILER})

# the compiler used for C++ files 
MESSAGE(STATUS "CMAKE_CXX_COMPILER:			" ${CMAKE_CXX_COMPILER})

# if the compiler is a variant of gcc, this should be set to 1 
MESSAGE(STATUS "CMAKE_COMPILER_IS_GNUCC:	" ${CMAKE_COMPILER_IS_GNUCC})

# if the compiler is a variant of g++, this should be set to 1 
MESSAGE(STATUS "CMAKE_COMPILER_IS_GNUCXX:	" ${CMAKE_COMPILER_IS_GNUCXX})

# the tools for creating libraries 
MESSAGE(STATUS "CMAKE_AR:					" ${CMAKE_AR})
MESSAGE(STATUS "CMAKE_RANLIB:				" ${CMAKE_RANLIB})

MESSAGE(STATUS "CMAKE_CROSSCOMPILING:		" ${CMAKE_CROSSCOMPILING})

MESSAGE(STATUS "TARGET_ARCH:				" ${TARGET_ARCH})

MESSAGE(STATUS "PANDA_EXTERNAL_ROOT: " ${PANDA_EXTERNAL_ROOT})
MESSAGE(STATUS "PANDA_EXTERNAL_DEBUG_LIB: " ${PANDA_EXTERNAL_DEBUG_LIB})
MESSAGE(STATUS "PANDA_EXTERNAL_RELEASE_LIB: " ${PANDA_EXTERNAL_RELEASE_LIB})


# ------------------------- End of Generic CMake Variable Logging ------------------

