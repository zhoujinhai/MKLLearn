cmake_minimum_required(VERSION 3.10)

# set the project name
project(Tutorial VERSION 1.0)

# set verbose makefile
set(CMAKE_VERBOSE_MAKEFILE ON)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# MKL
set(MKL_DIR "D:/work/heygears/code/algorithm--master/thirdParty/mkl")
include_directories(${MKL_DIR}/include)
link_directories(${MKL_DIR}/lib/WIN64)

# add the executable
add_executable(${PROJECT_NAME} linearAlgebra.cpp)
target_link_libraries(${PROJECT_NAME} mkl_core_dll mkl_intel_lp64_dll mkl_tbb_thread_dll)


