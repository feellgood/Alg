set(PROJECT_MAJOR_VERSION 1)
set(PROJECT_MINOR_VERSION 4)
set(PROJECT_PATCH_VERSION 0)
#===========================================================================

if(POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW) # CMake 3.0.0
endif()

PROJECT (ALG VERSION ${PROJECT_MAJOR_VERSION}.${PROJECT_MINOR_VERSION}.${PROJECT_PATCH_VERSION} )
cmake_minimum_required(VERSION 3.0)

#  LANGUAGES CXX CUDA
ENABLE_LANGUAGE(CXX CUDA)

MESSAGE(STATUS "CXX ${CMAKE_CXX_COMPILER_ID}")
MESSAGE(STATUS "version = ${ALG_VERSION}")

add_subdirectory(${PROJECT_SOURCE_DIR}/src)

set( CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++17 -O3 -march=native" CACHE STRING "C++ compiler flags used during builds." FORCE )

set(HEADERS fem.h)
set(SOURCES main.cpp lecture.cpp femutil.cpp chapeaux.cpp affichage.cpp solve.cpp savesol.cpp savevtk.cpp integrales.cpp)


ADD_EXECUTABLE(exe ${SOURCES} ${HEADERS})
add_dependencies(exe alg)
target_include_directories(exe PUBLIC src)
target_link_libraries(exe -ltbb alg)

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

#===========================================================================
# GPU cuda code, by default it is not built. you have to switch it ON while cmaking to build some executables using GPU :
#  cmake . -DENABLE_GPU=ON
# NB: ENABLE_GPU value is written in CMakeCache.txt
#===========================================================================

option(ENABLE_GPU "Enable GPU" OFF)

message( STATUS "enable GPU: ENABLE_GPU = "  ${ENABLE_GPU} )
if(${ENABLE_GPU})
    add_subdirectory(src/GPU)
endif()

