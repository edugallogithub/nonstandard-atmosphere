cmake_minimum_required(VERSION 3.6)			# Cmake version check
project(atm)								# create project "phd"
set(CMAKE_CXX_STANDARD 11)            		# enable c++11 standard

# ONLY ONE LINE ACTIVE BELOW
#set(CMAKE_BUILD_TYPE debug CACHE STRING "Choose debug or release." FORCE)
set(CMAKE_BUILD_TYPE release CACHE STRING "Choose debug or release." FORCE)
IF (WIN32)
	set(CMAKE_CXX_FLAGS "/MP /EHsc" CACHE STRING "Value given by jmake. Posible override though env variable." FORCE)
	set(CMAKE_CXX_FLAGS_DEBUG "/MDd /ZI /Od" CACHE STRING "Value given by jmake. Posible override though env variable." FORCE)
	set(CMAKE_CXX_FLAGS_RELEASE "/MD" CACHE STRING "Value given by jmake. Posible override though env variable." FORCE)
ELSE()
	set(CMAKE_CXX_FLAGS "-fPIC -fvisibility=hidden -fvisibility-inlines-hidden " CACHE STRING "Add if required.." FORCE)
	set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3" CACHE STRING "Debug options." FORCE)
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g0" CACHE STRING "Release options." FORCE)
ENDIF()
set(DEBUG_TARGET_SUFFIX "d" CACHE STRING "Add to targets as suffix.")

set(PROJECTS_SRC_ROOT ${CMAKE_SOURCE_DIR}) # full path to top level of current CMake source tree
set(PROJECTS_BIN_ROOT ${CMAKE_BINARY_DIR}) # full path to top level of current CMake build tree
set(BUILD_SANDBOX ${PROJECTS_BIN_ROOT}/sandbox) # sandbox path one level below build tree
file(MAKE_DIRECTORY ${BUILD_SANDBOX}) # create sandbox directory (and its parents if required)
file(WRITE ${BUILD_SANDBOX}/hook "this file is for cmake to find this directory") # create "hook" file within sandbox
find_path(EXECUTABLE_OUTPUT_PATH hook ${BUILD_SANDBOX}) # find path where "hook" file is located
find_path(LIBRARY_OUTPUT_PATH hook ${BUILD_SANDBOX}) # find path where "hook" file is located
link_directories(${LIBRARY_OUTPUT_PATH}) # specify paths where linker will look for libraries
set(CMAKE_MODULE_PATH $ENV{CMAKE_MODULE_PATH}) # set cache variable equal to environment variable

#=================================================================================================================
# Build Detection 32/64-bit
#=================================================================================================================
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(BUILD_ARCH 64)
    set(BUILD_ARCH_SUFFIX "")
    message("-- Detected 64-bit Build")
elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(BUILD_ARCH 32)
    set(BUILD_ARCH_SUFFIX "32")
    set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB32_PATHS TRUE)
    message("-- Notice! Detected 32-bit Build")
else()
    message(FATAL_ERROR "-- Unknown Build Architecture: contact with your administrator")
endif()

add_subdirectory(atm) # adds subdirectory to look for source CMakeLists.txt files










