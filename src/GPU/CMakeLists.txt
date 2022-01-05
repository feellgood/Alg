PROJECT(ALG)
cmake_minimum_required(VERSION 3.8.2)

find_package(CUDA REQUIRED)

enable_language(CUDA)

add_executable( hello hello.cu )
set_property( TARGET hello PROPERTY CUDA_STANDARD 11)

add_executable( add add.cu )
set_property( TARGET add PROPERTY CUDA_STANDARD 11)
