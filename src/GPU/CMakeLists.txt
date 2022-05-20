PROJECT(ALG LANGUAGES CXX CUDA)
cmake_minimum_required(VERSION 3.16)

find_package( CUDA 10.1 REQUIRED) # not sure if we really need version 10.1
if (CUDA_FOUND)

	set(CMAKE_CXX_STANDARD 11)
	set(CMAKE_CUDA_STANDARD 11)

	#message("CMAKE_CUDA_ARCHITECTURES : ${CMAKE_CUDA_ARCHITECTURES} ") #does not work : always return 30 (default value)

	#here we detect what is the installed GPU, this is undocumented, and might change with future version of cmake.
	include(FindCUDA/select_compute_arch)
	CUDA_DETECT_INSTALLED_GPUS(STR_INSTALLED_GPU)
	
	message("gpu found : ${STR_INSTALLED_GPU}")
	string(STRIP "${STR_INSTALLED_GPU}" INSTALLED_GPU)
	
	
#	string(FIND "${INSTALLED_GPU}" ";" best_arch REVERSE )
#	string(FIND "${INSTALLED_GPU}" "+PTX" nestor )
#	math(EXPR my_len  "${nestor} - ${best_arch}")
#	string(SUBSTRING "${INSTALLED_GPU}" ${best_arch} ${my_len} str_arch  )
#	string(REPLACE ";" "" archi ${str_arch})
		
	string(REPLACE "." "" CUDA_ARCH ${INSTALLED_GPU})
	message("CUDA_ARCH = ${CUDA_ARCH}")

	set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH}) # do i really need that ?

	add_executable( add add.cu )
	set_property( TARGET add PROPERTY CUDA_STANDARD 11)
	set_target_properties( add PROPERTIES CUDA_ARCHITECTURES ${CUDA_ARCH} )

	add_executable( pscal pscal.cu )
	set_property( TARGET pscal PROPERTY CUDA_STANDARD 11)
	set_target_properties( pscal PROPERTIES CUDA_ARCHITECTURES ${CUDA_ARCH} )
	target_compile_options(pscal PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -std=c++11; --expt-extended-lambda   > )

	add_executable( sparseVect sparseVect.cu )
	set_property( TARGET sparseVect PROPERTY CUDA_STANDARD 11)
	set_target_properties( sparseVect PROPERTIES CUDA_ARCHITECTURES ${CUDA_ARCH} )
	target_compile_options(sparseVect PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -std=c++11; --expt-extended-lambda   > )


	add_library(alg_GPU STATIC infos.cu cg.cu precond_cg.cu alg_GPU.h)
	target_compile_features ( alg_GPU PUBLIC cxx_std_11 )
	set_target_properties( alg_GPU PROPERTIES
						LINKER_LANGUAGE CUDA
						CUDA_ARCHITECTURES ${CUDA_ARCH}
						#POSITION_INDEPENDENT_CODE ON # for -fPIC option on
						CUDA_SEPARABLE_COMPILATION ON
						#CUDA_RESOLVE_DEVICE_SYMBOLS ON
						#IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/libalg_GPU.a"
						)
	target_link_libraries( alg_GPU ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY} )


	add_executable( test_cg test_cg.cu alg_GPU.h )
	set_target_properties( test_cg PROPERTIES
						LINKER_LANGUAGE CUDA
						CUDA_ARCHITECTURES ${CUDA_ARCH} 
						#POSITION_INDEPENDENT_CODE ON # for -fPIC option on
						CUDA_SEPARABLE_COMPILATION ON
						)
	target_link_libraries( test_cg PRIVATE alg_GPU )
else()
	message("Error : CUDA not found.")
endif(CUDA_FOUND)

