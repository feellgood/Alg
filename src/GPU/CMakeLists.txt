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
