cmake_minimum_required(VERSION 3.15)
project(wannier_project)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "PREFIX/findFFTW")

set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
#find_package(FFTW REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -g -Wall -Wextra -ffast-math  -march=native -fopenmp")
set(CMAKE_C_FLAGS "-O3 -g -Wall -Wextra -ffast-math  -march=native -fopenmp")
set(LAPACK_DIR "/opt/homebrew/Cellar/lapack/3.12.1/")
set(OPENBLAS_DIR "/opt/homebrew/Cellar/openblas/0.3.32/")
set(FFTW_DIR "/opt/homebrew/Cellar/fftw/3.3.10_3/")

include_directories(include
	   	    ${LAPACK_DIR}/include
	   	    ${OPENBLAS_DIR}/include
		    ${FFTW_DIR}/include)

link_directories(${LAPACK_DIR}/lib)
link_directories(${OPENBLAS_DIR}/lib)
link_directories(${FFTW_DIR}/lib)
add_library(wannier_src OBJECT
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

set_property(TARGET wannier_src PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(wannier SHARED $<TARGET_OBJECTS:wannier_src>)
target_link_libraries(wannier PUBLIC lapacke
			          LAPACK::LAPACK
			          BLAS::BLAS
				  fftw3)

add_executable(main src/main.cpp)
target_link_libraries(main
		      wannier) 

make_directory(tmp)
