set(PROJECT_MAJOR_VERSION 0)
set(PROJECT_MINOR_VERSION 0)
set(PROJECT_PATCH_VERSION 0)
#===========================================================================

if(POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW) # CMake 3.0.0
endif()

PROJECT (ALG VERSION ${PROJECT_MAJOR_VERSION}.${PROJECT_MINOR_VERSION}.${PROJECT_PATCH_VERSION} )
cmake_minimum_required(VERSION 3.0)

ENABLE_LANGUAGE(CXX)

MESSAGE(STATUS "CXX ${CMAKE_CXX_COMPILER_ID}")

MESSAGE(STATUS "version = ${ALG_VERSION}")

#===========================================================================
# Compiler options -DNDEBUG
# -fsanitize=address leak  null (to debug memory leak) undefined
#===========================================================================

set( CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++11 -O3 -march=native" CACHE STRING "C++ compiler flags used during builds." FORCE )

set(HEADERS alg.h alg_iter.h)
set(SOURCES alg_cg_dir.cpp)

ADD_EXECUTABLE(exe main.cpp ${HEADERS} ${SOURCES})

#===========================================================================
# unit tests, by default it is not built. you have to switch it ON while cmaking to build the unit tests and to be able to call 'make test' afterward:
#  cmake . -DENABLE_UTESTS=ON
# NB: ENABLE_UTESTS value is written in CMakeCache.txt
#===========================================================================

option(ENABLE_UTESTS "Enable tests" OFF)

message( STATUS "enable unit tests: ENABLE_UTESTS = "  ${ENABLE_UTESTS} )
if(${ENABLE_UTESTS})
    include(CTest)
    add_subdirectory(tests)
endif()
