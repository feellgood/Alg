PROJECT(ALG LANGUAGES CXX CUDA)
cmake_minimum_required(VERSION 3.16)

find_package( CUDA 10.1 REQUIRED)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)

include(FindCUDA/select_compute_arch)
CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU)
string(STRIP ${INSTALLED_GPU} INSTALLED_GPU2)
string(REPLACE " " ";" INSTALLED_GPU3 "${INSTALLED_GPU2}")
string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU3}" )

message("CUDA_ARCH = ${CUDA_ARCH_LIST}")

add_executable( hello hello.cu )
set_property( TARGET hello PROPERTY CUDA_STANDARD 11)
set_target_properties( hello PROPERTIES CUDA_ARCHITECTURES ${CUDA_ARCH_LIST} )

add_executable( add add.cu )
set_property( TARGET add PROPERTY CUDA_STANDARD 11)
set_target_properties( add PROPERTIES CUDA_ARCHITECTURES ${CUDA_ARCH_LIST} )

add_executable( pscal pscal.cu )
set_property( TARGET pscal PROPERTY CUDA_STANDARD 11)
set_target_properties( pscal PROPERTIES CUDA_ARCHITECTURES ${CUDA_ARCH_LIST} )
target_compile_options(pscal PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -std=c++11; --expt-extended-lambda   > )

add_executable( sparseVect sparseVect.cu )
set_property( TARGET sparseVect PROPERTY CUDA_STANDARD 11)
set_target_properties( sparseVect PROPERTIES CUDA_ARCHITECTURES ${CUDA_ARCH_LIST} )
target_compile_options(sparseVect PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -std=c++11; --expt-extended-lambda   > )


add_library(alg_GPU STATIC infos.cu cg.cu alg_GPU.h)
target_compile_features ( alg_GPU PUBLIC cxx_std_11 )
set_target_properties( alg_GPU PROPERTIES
					LINKER_LANGUAGE CUDA
					CUDA_ARCHITECTURES ${CUDA_ARCH_LIST}
					#POSITION_INDEPENDENT_CODE ON # for -fPIC option on
					CUDA_SEPARABLE_COMPILATION ON
					#CUDA_RESOLVE_DEVICE_SYMBOLS ON
					#IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/libalg_GPU.a"
					)
target_link_libraries( alg_GPU ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY} )


#link_directories( ${CMAKE_SOURCE_DIR}/src/GPU )
add_executable( test_cg test_cg.cu alg_GPU.h )
set_target_properties( test_cg PROPERTIES
					LINKER_LANGUAGE CUDA
					CUDA_ARCHITECTURES ${CUDA_ARCH_LIST} 
					#POSITION_INDEPENDENT_CODE ON # for -fPIC option on
					CUDA_SEPARABLE_COMPILATION ON
					
					)
target_link_libraries( test_cg PRIVATE alg_GPU )

