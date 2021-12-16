set(PROJECT_MAJOR_VERSION 1)
set(PROJECT_MINOR_VERSION 4)
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

set( CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++17 -O3 -march=native" CACHE STRING "C++ compiler flags used during builds." FORCE )

set(HEADERS alg.h alg_core.h alg_iter.h alg_coeff.h alg_sparseVect.h alg_sparseMat.h alg_denseMat.h alg_lu.h alg_cg.h alg_cg_dir.h alg_cg_ilu.h alg_cg_ilu_dir.h alg_bicg.h alg_bicg_dir.h alg_bicg_ilu.h alg_bicg_ilu_dir.h)
set(SOURCES main.cpp lecture.cpp femutil.cpp chapeaux.cpp affichage.cpp solve.cpp savesol.cpp savevtk.cpp integrales.cpp fem.h)

ADD_EXECUTABLE(exe ${SOURCES})
target_include_directories(exe PUBLIC src)
target_link_libraries(exe -ltbb)

install(FILES src/${HEADERS} 
	DESTINATION include/alg
	COMPONENT headers)

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
