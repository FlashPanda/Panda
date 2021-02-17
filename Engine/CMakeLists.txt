
# Core
FILE(GLOB CORE_SOURCE 
	./Core/*.*
	./Core/Interface/*.*
	./Core/Math/*.*
	./Core/Parser/*.*
	)

# Visual Studio source folders
SOURCE_GROUP(TREE ${PANDA_ROOT} FILES ${CORE_SOURCE})

# RHI
IF(${WIN32})
	FILE(GLOB RHI_SOURCE	
		./RHI/*.*
		./RHI/D3D/*.*)
ELSE (${WIN32})
	FILE(GLOB RIH_SOURCE
		./RHI/*.*
		./RHI/OpenGL/*.*)
ENDIF (${WIN32})
SOURCE_GROUP(TREE ${PANDA_ROOT} FILES ${RHI_SOURCE})

# SceneManager
FILE(GLOB SCENE_MANAGER_SOURCE	./SceneManager/*.*)
SOURCE_GROUP(TREE ${PANDA_ROOT} FILES ${SCENE_MANAGER_SOURCE})

# DebugManager
FILE(GLOB DEBUG_MANAGER_SOURCE	./DebugManager/*.*)
SOURCE_GROUP(TREE ${PANDA_ROOT} FILES ${DEBUG_MANAGER_SOURCE})

# InputManager
FILE(GLOB INPUT_MANAGER_SOURCE	./InputManager/*.*)
SOURCE_GROUP(TREE ${PANDA_ROOT} FILES ${INPUT_MANAGER_SOURCE})

#Common
FILE(GLOB COMMON_SOURCE ./Common/*.*)
IF(${WIN32})
	set (COMMON_WINDOWS 
		./Common/Entry/main.cpp
		./Common/Entry/BaseApplication.hpp
		./Common/Entry/BaseApplication.cpp
		./Common/Entry/WindowsApplication.hpp
		./Common/Entry/WindowsApplication.cpp
		./Common/Entry/WindowsD3D12Application.hpp
		./Common/Entry/WindowsD3D12Application.cpp
		./Common/Entry/WindowsD3D12ApplicationConfig.cpp)
ELSE (${WIN32})
	set (COMMON_LINUX
		./Common/Entry/main.cpp
		./Common/Entry/BaseApplication.hpp
		./Common/Entry/BaseApplication.cpp
		./Common/Entry/LinuxOpenGLApplication.hpp
		./Common/Entry/LinuxOpenGLApplication.cpp
		./Common/Entry/LinuxOpenGLApplicationConfig.cpp

		./RHI/OpenGL/glad/src/glad_glx.c)
	include_directories(./RHI/OpenGL/glad/include)
ENDIF(${WIN32})
SOURCE_GROUP(TREE ${PANDA_ROOT} FILES ${COMMON_SOURCE})
IF(${WIN32})
	SOURCE_GROUP(TREE ${PANDA_ROOT} FILES ${COMMON_WINDOWS})
ELSE (${WIN32})
	SOURCE_GROUP(TREE ${PANDA_ROOT} FILES ${COMMON_LINUX})
ENDIF (${WIN32})

# package lib
IF(${WIN32})
	add_library(Panda STATIC
		${CORE_SOURCE}
		${RHI_SOURCE}
		${SCENE_MANAGER_SOURCE}
		${DEBUG_MANAGER_SOURCE}
		${INPUT_MANAGER_SOURCE}
		${COMMON_SOURCE}
		${COMMON_WINDOWS})
		
	# find debug libs
	find_library(CROSSGUID_DEBUG crossguid PATHS ${PANDA_EXTERNAL_DEBUG_LIB} NO_DEFAULT_PATH)
	find_library(ASSIMP_DEBUG assimp PATHS ${PANDA_EXTERNAL_DEBUG_LIB} NO_DEFAULT_PATH)
	find_library(LIBPNG_DEBUG libpng PATHS ${PANDA_EXTERNAL_DEBUG_LIB} NO_DEFAULT_PATH)
	find_library(OPENDDL_DEBUG OpenDDL PATHS ${PANDA_EXTERNAL_DEBUG_LIB} NO_DEFAULT_PATH)
	find_library(OPENGEX_DEBUG OpenGEX PATHS ${PANDA_EXTERNAL_DEBUG_LIB} NO_DEFAULT_PATH)
	find_library(TURBOJPEG_DEBUG turbojpeg PATHS ${PANDA_EXTERNAL_DEBUG_LIB} NO_DEFAULT_PATH)
	find_library(ZLIB_DEBUG zlib PATHS ${PANDA_EXTERNAL_DEBUG_LIB} NO_DEFAULT_PATH)

	target_link_libraries(Panda debug ${CROSSGUID_DEBUG})
	target_link_libraries(Panda debug ${ASSIMP_DEBUG})
	target_link_libraries(Panda debug ${LIBPNG_DEBUG})
	target_link_libraries(Panda debug ${OPENDDL_DEBUG})
	target_link_libraries(Panda debug ${OPENGEX_DEBUG})
	target_link_libraries(Panda debug ${TURBOJPEG_DEBUG})
	target_link_libraries(Panda debug ${ZLIB_DEBUG})
		
	# find release libs
	find_library(CROSSGUID_RELEASE crossguid PATHS ${PANDA_EXTERNAL_RELEASE_LIB} NO_CMAKE_FIND_ROOT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
	find_library(ASSIMP_RELEASE assimp PATHS ${PANDA_EXTERNAL_RELEASE_LIB} NO_CMAKE_FIND_ROOT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
	find_library(LIBPNG_RELEASE libpng PATHS ${PANDA_EXTERNAL_RELEASE_LIB} NO_CMAKE_FIND_ROOT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
	find_library(OPENDDL_RELEASE OpenDDL PATHS ${PANDA_EXTERNAL_RELEASE_LIB} NO_CMAKE_FIND_ROOT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
	find_library(OPENGEX_RELEASE OpenGEX PATHS ${PANDA_EXTERNAL_RELEASE_LIB} NO_CMAKE_FIND_ROOT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
	find_library(TURBOJPEG_RELEASE turbojpeg PATHS ${PANDA_EXTERNAL_RELEASE_LIB} NO_CMAKE_FIND_ROOT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
	find_library(ZLIB_RELEASE zlib PATHS ${PANDA_EXTERNAL_RELEASE_LIB} NO_CMAKE_FIND_ROOT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
	
	target_link_libraries(Panda optimized ${CROSSGUID_RELEASE})
	target_link_libraries(Panda optimized ${ASSIMP_RELEASE})
	target_link_libraries(Panda optimized ${LIBPNG_RELEASE})
	target_link_libraries(Panda optimized ${OPENDDL_RELEASE})
	target_link_libraries(Panda optimized ${OPENGEX_RELEASE})
	target_link_libraries(Panda optimized ${TURBOJPEG_RELEASE})
	target_link_libraries(Panda optimized ${ZLIB_RELEASE})
	
	
ELSE (${WIN32})

ENDIF(${WIN32})

set_property(TARGET  Panda PROPERTY FOLDER "Engine")