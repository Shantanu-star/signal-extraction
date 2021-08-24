project(SignalExtraction CXX)
cmake_minimum_required(VERSION 3.10)

if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
endif()

# in DEBUG mode make verbose Makefile
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif ()

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -DDEBUG -D__DEBUG -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ftree-vectorize -ffast-math -DNODEBUG")
message(STATUS "Using CXX flags for ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")

find_package(ROOT REQUIRED)

message(STATUS "Using ROOT: ${ROOT_VERSION} <${ROOT_CONFIG}>")
include_directories(${ROOT_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR} src)
include(${ROOT_USE_FILE})

add_subdirectory(src)
add_subdirectory(macro)