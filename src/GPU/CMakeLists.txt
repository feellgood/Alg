PROJECT(ALG)
cmake_minimum_required(VERSION 3.16)

find_package(CUDA REQUIRED)

enable_language(CUDA)

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
target_compile_options(pscal PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
-std=c++11;
--expt-extended-lambda   > )

add_executable( sparseVect sparseVect.cu )
set_property( TARGET sparseVect PROPERTY CUDA_STANDARD 11)
set_target_properties( sparseVect PROPERTIES CUDA_ARCHITECTURES ${CUDA_ARCH_LIST} )
target_compile_options(sparseVect PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
-std=c++11;
--expt-extended-lambda   > )

add_executable( test_cg test_cg.cu )
set_property( TARGET test_cg PROPERTY CUDA_STANDARD 11)
set_target_properties( test_cg PROPERTIES CUDA_ARCHITECTURES ${CUDA_ARCH_LIST} )
target_link_libraries( test_cg ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY} )
target_compile_options( test_cg PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
-std=c++11;
--expt-extended-lambda   > )

