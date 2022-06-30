PROJECT(ALG LANGUAGES CXX )
cmake_minimum_required(VERSION 3.16)

set(ALG_HEADERS alg.h alg_core.h alg_iter.h alg_coeff.h alg_sparseVect.h alg_sparseMat.h alg_denseMat.h alg_lu.h alg_cg.h alg_cg_dir.h alg_cg_ilu.h alg_cg_ilu_dir.h alg_bicg.h alg_bicg_dir.h alg_bicg_ilu.h alg_bicg_ilu_dir.h)

file(GLOB ALG_SOURCES *.cpp )

add_library(alg STATIC ${ALG_HEADERS} ${ALG_SOURCES})
	target_compile_features ( alg PUBLIC cxx_std_17 )
	set_target_properties( alg PROPERTIES
						LINKER_LANGUAGE CXX
						#POSITION_INDEPENDENT_CODE ON # for -fPIC option on
						)

install(FILES src/${ALG_HEADERS} 
	DESTINATION include/alg
	COMPONENT headers)

