cmake_minimum_required(VERSION 3.16)
project(CFDSolutions LANGUAGES C CXX)

# Enable Fortran only when a compiler is available; otherwise keep the C/C++ build path usable.
set(ENABLE_FORTRAN OFF)
if(CMAKE_Fortran_COMPILER)
  set(FORTRAN_COMPILER "${CMAKE_Fortran_COMPILER}")
  set(ENABLE_FORTRAN ON)
elseif(CMAKE_HOST_WIN32)
  set(_fortran_hint_paths)
  if(DEFINED ENV{FORTRAN_ROOT} AND NOT "$ENV{FORTRAN_ROOT}" STREQUAL "")
    list(APPEND _fortran_hint_paths "$ENV{FORTRAN_ROOT}")
  endif()
  list(APPEND _fortran_hint_paths
    "C:/msys64/ucrt64/bin"
    "C:/msys64/mingw64/bin"
    "C:/msys64/clang64/bin"
    "C:/msys64/mingw32/bin"
    "C:/msys64/usr/bin"
    "C:/Program Files/WinLibs/bin"
    "C:/Program Files/WinLibs/mingw64/bin"
  )
  find_program(FORTRAN_COMPILER NAMES gfortran flang-new flang ifort f95 f77 HINTS ${_fortran_hint_paths})
  if(FORTRAN_COMPILER)
    set(ENABLE_FORTRAN ON)
  else()
    message(STATUS "No Fortran compiler detected; Fortran sources will be skipped")
  endif()
else()
  set(ENABLE_FORTRAN ON)
endif()

if(FORTRAN_COMPILER)
  message(STATUS "Using Fortran compiler: ${FORTRAN_COMPILER}")
endif()

if(ENABLE_FORTRAN)
  enable_language(Fortran)
endif()

option(ENABLE_MPI "Build with MPI support" OFF)
option(ENABLE_NETCDF "Enable NetCDF support" OFF)
option(ENABLE_CGNS "Enable CGNS support" OFF)
option(BUILD_SHARED_LIBS "Build library dependencies as shared libraries" OFF)

if(ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS C CXX)
endif()

if(ENABLE_NETCDF)
  find_package(NetCDF REQUIRED)
endif()

if(ENABLE_CGNS)
  find_path(CGNS_INCLUDE_DIR cgnslib.h
    HINTS
      ${CMAKE_PREFIX_PATH}
      /usr/include
      /usr/local/include
      "C:/Program Files/cgns/include"
      "C:/Program Files (x86)/cgns/include"
  )
  find_library(CGNS_LIBRARY NAMES cgns cg
    HINTS
      ${CMAKE_PREFIX_PATH}
      /usr/lib
      /usr/local/lib
      /usr/lib64
      /usr/local/lib64
      "C:/Program Files/cgns/lib"
      "C:/Program Files (x86)/cgns/lib"
  )

  if(NOT CGNS_INCLUDE_DIR OR NOT CGNS_LIBRARY)
    message(FATAL_ERROR "ENABLE_CGNS is ON, but CGNS headers or library were not found.")
  endif()
endif()

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# All compiler/optimization/warning flags live in one place.
include(cmake/CompilerFlags.cmake)

# Export compile commands for editor/IDE tooling
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(UTILS)
add_subdirectory(MATH)
add_subdirectory(CFDDB)

if(ENABLE_FORTRAN)
  add_subdirectory(NISTThermo/NISTThermo)
  add_subdirectory(EOS)
  add_subdirectory(MPFLOW3D)
  message(STATUS "Enabled optional Fortran-backed thermodynamics and solver targets")
else()
  message(STATUS "Skipping optional Fortran-backed thermodynamics and solver targets because Fortran support is disabled")
endif()

add_subdirectory(VERDICT)
add_subdirectory(Euler2D)
add_subdirectory(Euler3D)
add_subdirectory(Mesh2D)
add_subdirectory(MESHQA)
add_subdirectory(ParEuler3D)
add_subdirectory(PostCFD)
add_subdirectory(UGRID2SU2)
