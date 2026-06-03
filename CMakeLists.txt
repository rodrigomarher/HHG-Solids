cmake_minimum_required(VERSION 3.15)
project(wannier_project LANGUAGES C CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)
find_package(FFTW REQUIRED)
find_package(CUDA QUIET)
find_package(CUDAToolkit QUIET)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -g -Wall -Wextra -ffast-math  -march=native -fopenmp")
set(CMAKE_C_FLAGS "-O3 -g -Wall -Wextra -ffast-math  -march=native -fopenmp")
#set(CMAKE_CXX_FLAGS "-fsanitize=address -fsanitize=undefined -g -O1 -fno-omit-frame-pointer")
#set(CMAKE_C_FLAGS "-fsanitize=address -fsanitize=undefined -g -O1 -fno-omit-frame-pointer")
#set(LAPACK_DIR "/opt/homebrew/Cellar/lapack/3.12.1/")
#set(OPENBLAS_DIR "/opt/homebrew/Cellar/openblas/0.3.32/")
#set(FFTW_DIR "/opt/homebrew/Cellar/fftw/3.3.10_3/")

find_path(
    LAPACKE_INCLUDE_DIR
    lapacke.h
    REQUIRED
)
find_path(
    LAPACK_INCLUDE_DIR
    lapack.h
    REQUIRED
)

include_directories(include
		    ${LAPACK_INCLUDE_DIR}
		    ${LAPACKE_INCLUDE_DIR}
	   	    ${OPENBLAS_DIR}/include
		    ${FFTW_DIR}/include)
#
#link_directories(${LAPACK_DIR}/lib)
link_directories(${OPENBLAS_DIR}/lib)
link_directories(${FFTW_DIR}/lib)

set(wannier_src
    src/settings.cpp
    src/vec3_util.cpp
    src/wannier_tb.cpp
    src/grid.cpp
    src/matrix_field.cpp
    src/hamiltonian.cpp
    src/swe.cpp
    src/operator.cpp
    src/rdm.cpp
    src/efield.cpp
    src/berry_connection.cpp
    src/velocity.cpp
    src/solver.cpp
    src/solver_kspace.cpp
    src/fftw_helper.cpp
    src/observable.cpp
    src/interpolation.cpp
    src/cwrapper.cpp
)

if(CUDA_FOUND)
    enable_language(CUDA)
    add_definitions(-DHAVE_CUDA)
    set(cuda_src
	src/cuda_kernels.cu
	src/solver_cuda.cpp)
    add_library(wannier_cuda ${cuda_src})
    include_directories(${CUDA_INCLUDE_DIRS})
    add_library(wannier SHARED ${wannier_src} ${cuda_src})
else()
    add_library(wannier SHARED ${wannier_src})
endif()


set_property(TARGET wannier PROPERTY POSITION_INDEPENDENT_CODE 1)

if(CUDA_FOUND)
    target_link_libraries(wannier PUBLIC 
			          BLAS::BLAS
				  CUDA::cufft_static
				  fftw3)
    
    set_target_properties(wannier PROPERTIES
		          CUDA_SEPARABLE_COMPILATION ON)

else()
    target_link_libraries(wannier PUBLIC 
			          BLAS::BLAS				  
				  fftw3)
endif()

add_executable(main src/main.cpp)

target_link_libraries(main
		      wannier) 

make_directory(tmp)
