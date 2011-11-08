PROJECT(TARANG_MPI)
cmake_minimum_required(VERSION 2.6)

SUBDIRS(lib src)

INCLUDE_DIRECTORIES(
${CMAKE_SOURCE_DIR}/lib/basis_basicfn
${CMAKE_SOURCE_DIR}/lib/scft
${CMAKE_SOURCE_DIR}/lib/four
${CMAKE_SOURCE_DIR}/lib/plainfields
${CMAKE_SOURCE_DIR}/lib/fields
${CMAKE_SOURCE_DIR}/lib/IncFlow
${CMAKE_SOURCE_DIR}/lib/IncFluid
)


LINK_DIRECTORIES(
${CMAKE_SOURCE_DIR}/lib/basis_basicfn
${CMAKE_SOURCE_DIR}/lib/scft
${CMAKE_SOURCE_DIR}/lib/four
${CMAKE_SOURCE_DIR}/lib/plainfields
${CMAKE_SOURCE_DIR}/lib/fields
${CMAKE_SOURCE_DIR}/lib/IncFlow
${CMAKE_SOURCE_DIR}/lib/IncFluid
${CMAKE_SOURCE_DIR}/src
${CMAKE_SOURCE_DIR}/src/IFluid
${CMAKE_SOURCE_DIR}/src/IMHD
${CMAKE_SOURCE_DIR}/src/IScalar
${CMAKE_SOURCE_DIR}/src/RBslip
${CMAKE_SOURCE_DIR}/src/RBslipMHD
)


ADD_EXECUTABLE(tarangMPI  src/main.h src/main.cc
	src/read_prog_para.cc
	src/read_para.cc
	src/read_diag_para.cc
	src/IFluid/Ifluid_main.h
	src/IFluid/Ifluid_main.cc
	src/IFluid/Ifluid_diag_main.cc
	src/IMHD/IMHD_main.h
	src/IMHD/IMHD_main.cc
	src/IMHD/IMHD_diag_main.cc
	src/IScalar/Iscalar_main.h
	src/IScalar/Iscalar_main.cc
	src/IScalar/Iscalar_diag_main.cc
	src/NB/NonBoussinesq_main.cc
	src/NB/NonBoussinesq_main.h
	src/RBslip/RB_slip_main.h
	src/RBslip/RB_slip_main.cc
	src/RBslip/RB_slip_diag_main.cc
	src/RBslipMHD/RB_slipMHD_main.h
	src/RBslipMHD/RB_slipMHD_main.cc
	)


find_package(MPI REQUIRED)

TARGET_LINK_LIBRARIES(tarangMPI IncFluid IncFlow fields plainfields four scft basis_basic m blitz fftw3_mpi fftw3 mpi)

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

SET(CMAKE_BUILD_TYPE distribution)
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -g")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3")
