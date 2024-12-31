# CMAKE_MSVC_RUNTIME_LIBRARY requires 3.15
cmake_minimum_required(VERSION 3.20...3.26)

project(
  PhreeqcRM
  VERSION 3.8.0
  LANGUAGES CXX C
)

if(DEFINED SKBUILD)
  include(CheckPIESupported)
  check_pie_supported()
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# check if this is the root project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(STANDALONE_BUILD 1)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/MP>)
else()
  set(STANDALONE_BUILD 0)
endif()

# overide docdir on windows
if(WIN32 AND NOT CMAKE_INSTALL_DOCDIR)
  set(CMAKE_INSTALL_DOCDIR "" CACHE PATH "documentation root (doc)")
  set(CMAKE_INSTALL_DOCDIR "doc")
endif()

include(GNUInstallDirs)

option(PHREEQCRM_FORTRAN_TESTING "Build Fortran test" OFF)
if(PHREEQCRM_FORTRAN_TESTING)
  enable_language(Fortran)
endif()

if(DEFINED SKBUILD)
  # SKBUILD is defined when scikit-build-core is used to build python wheels (see pyproject.toml)

  # Include FetchContent module
  include(FetchContent)

  # Fetch yaml-cpp
  FetchContent_Declare(
    yaml-cpp
    URL https://github.com/jbeder/yaml-cpp/archive/refs/tags/0.8.0.tar.gz
  )

  # Make yaml-cpp available
  FetchContent_MakeAvailable(yaml-cpp)


  # Fetch zlib
  FetchContent_Declare(
    zlib
    URL https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz
  )

  # Don't install zlib in wheel
  set(SKIP_INSTALL_ALL ON)
  # Make zlib available
  FetchContent_MakeAvailable(zlib)
endif()

option(PHREEQCRM_WITH_YAML_CPP "Build with yaml-cpp support" OFF)
if(DEFINED SKBUILD)
  # SKBUILD is defined when scikit-build-core is used to build python wheels (see pyproject.toml)
  set(PHREEQCRM_WITH_YAML_CPP ON)
endif()
if(PHREEQCRM_WITH_YAML_CPP AND NOT DEFINED SKBUILD)
  find_package(yaml-cpp REQUIRED)
  if(yaml-cpp_FOUND)
    message(STATUS "Found yaml-cpp: TRUE (found version \"${yaml-cpp_VERSION}\")")
  endif()
endif()

option(PHREEQCRM_DISABLE_OPENMP "Disable openmp" OFF)
if(NOT PHREEQCRM_DISABLE_OPENMP)
  if(STANDALONE_BUILD EQUAL 1)
    # check for openmp
    find_package(OpenMP)
  endif()
endif()

# check for MPI
option(PHREEQCRM_BUILD_MPI "Build MPI configuration" OFF)
if(PHREEQCRM_BUILD_MPI)
  find_package(MPI REQUIRED)
endif()

option(PHREEQCRM_BUILD_PYTHON "Build python wrapper using SWIG" OFF)
if(DEFINED SKBUILD)
  # SKBUILD is defined when scikit-build-core is used to build python wheels (see pyproject.toml)
  set(PHREEQCRM_BUILD_PYTHON ON)
endif()
if(PHREEQCRM_BUILD_PYTHON)
  if(NOT PHREEQCRM_WITH_YAML_CPP)
    message(FATAL_ERROR "Building python wrapper requires PHREEQCRM_WITH_YAML_CPP.")
  endif()
  find_package(SWIG REQUIRED)
  ###find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development Development.Module NumPy)
  find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
  # Set Python_ROOT_DIR to specify which python to use
  # ie -DPython_ROOT_DIR=C:/Python38-x64
endif()

if(DEFINED SKBUILD)
  # see https://scikit-build-core.readthedocs.io/en/latest/cmakelists.html
  set(Python_SOABI ${SKBUILD_SOABI})
endif()

option(PHREEQCRM_NO_UTF8_ENCODING "Build without utf-8 encoding" OFF)

set(PhreeqcRM_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BMI_interface_C.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BMI_interface_C.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BMI_interface_F.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BMI_interface_F.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/bmi.hxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BMIPhreeqcRM.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BMIPhreeqcRM.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BMIVariant.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BMIVariant.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/CSelectedOutput.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/CSelectedOutput.hxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/CVar.hxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/Debug.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/ErrorReporter.hxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface_F.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface_F.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqcCallbacks.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqcLib.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/advection.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/basicsubs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cl1.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/Parser.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/Parser.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_base.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_base.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_exports.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_io.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_io.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/phrqtype.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/Utils.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/Utils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cvdense.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cvdense.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cvode.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cvode.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxKinetics.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxKinetics.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxMix.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxMix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dense.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dense.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Dictionary.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Dictionary.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dumper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dumper.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Exchange.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Exchange.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ExchComp.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ExchComp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasComp.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasComp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/gases.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasPhase.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasPhase.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/global_structures.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/input.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/integrate.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/inverse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolution.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolution.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolutionComp.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolutionComp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/isotopes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/kinetics.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/KineticsComp.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/KineticsComp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/mainsubs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NA.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NumKeyword.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NumKeyword.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector_serial.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector_serial.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/parse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PBasic.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PBasic.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/phqalloc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/phqalloc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Phreeqc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Phreeqc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io_output.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/pitzer_structures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/pitzer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblage.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblage.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblageComp.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblageComp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/prep.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Pressure.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Pressure.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/print.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Reaction.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Reaction.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/read.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ReadClass.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/readtr.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/runner.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/runner.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SelectedOutput.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SelectedOutput.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Serializer.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Serializer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sit.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/smalldense.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/smalldense.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Solution.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Solution.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SolutionIsotope.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SolutionIsotope.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/spread.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SS.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SS.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SSassemblage.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SSassemblage.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SScomp.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SScomp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/step.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBinList.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBinList.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/structures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialsmath.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialsmath.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialstypes.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Surface.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Surface.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceCharge.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceCharge.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceComp.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceComp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/System.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/System.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/tally.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Temperature.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Temperature.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/tidy.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/transport.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Use.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Use.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/UserPunch.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/UserPunch.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/thread.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/Var.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/Var.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/Version.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqcPhast.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqcPhast.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqcPhastLib.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqcPhastLib.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/IrmResult.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/PhreeqcRM.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/PhreeqcRM.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/RM_interface_C.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/RM_interface_C.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/RM_interface_F.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/RM_interface_F.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/RMVARS.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/VarManager.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/VarManager.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/YAML_interface_C.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/YAML_interface_C.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/YAML_interface_F.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/YAML_interface_F.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/YAMLPhreeqcRM.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/YAMLPhreeqcRM.h
)

# dll/shared objects
option(BUILD_SHARED_LIBS "Build Shared Libraries" OFF)

add_library(PhreeqcRM ${PhreeqcRM_SOURCES})
add_library(PhreeqcRM::PhreeqcRM ALIAS PhreeqcRM)

target_include_directories(PhreeqcRM
  PUBLIC
    $<BUILD_INTERFACE:${PhreeqcRM_BINARY_DIR}>
    $<BUILD_INTERFACE:${PhreeqcRM_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${PhreeqcRM_SOURCE_DIR}/src/IPhreeqcPhast>
    $<BUILD_INTERFACE:${PhreeqcRM_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc>
    $<BUILD_INTERFACE:${PhreeqcRM_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp>
    $<BUILD_INTERFACE:${PhreeqcRM_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common>
    $<BUILD_INTERFACE:${PhreeqcRM_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# generate "${PhreeqcRM_BINARY_DIR}/irm_dll_export.h"
# based on how the PhreeqcRM library is built (BUILD_SHARED_LIBS)
include(GenerateExportHeader)
generate_export_header(PhreeqcRM BASE_NAME IRM_DLL)

# c++11
target_compile_features(PhreeqcRM PUBLIC cxx_std_11)

# iphreeqc defs
target_compile_definitions(PhreeqcRM PRIVATE SWIG_SHARED_OBJ)
target_compile_definitions(PhreeqcRM PRIVATE USE_PHRQ_ALLOC)
target_compile_definitions(PhreeqcRM PRIVATE IPhreeqc_EXPORTS)

if (PHREEQCRM_NO_UTF8_ENCODING)
  target_compile_definitions(PhreeqcRM PRIVATE NO_UTF8_ENCODING)
endif()

# Disable deprecation warnings for standard C and STL functions in VS2005
# and later
if(MSVC_VERSION EQUAL 1400 OR MSVC_VERSION GREATER 1400)
  target_compile_definitions(PhreeqcRM PRIVATE _CRT_SECURE_NO_DEPRECATE)
  target_compile_definitions(PhreeqcRM PRIVATE _SCL_SECURE_NO_WARNINGS)
endif()

if (WIN32 AND BUILD_SHARED_LIBS)
  # Required to build IMPLIB
  # (Seems to be automatically set when using Visual Studio as the generator)
  # This is only necessary to export iphreeqc methods
  # iphreeqc should be updated to use generate_export_header
  target_compile_definitions(PhreeqcRM PRIVATE _WINDLL)
endif()

# init exportPackageDependencies
set(exportPackageDependencies "")

# yaml-cpp
if(PHREEQCRM_WITH_YAML_CPP AND yaml-cpp_FOUND)
  target_compile_definitions(PhreeqcRM PUBLIC USE_YAML)
  if(TARGET yaml-cpp::yaml-cpp)
    target_link_libraries(PhreeqcRM PUBLIC yaml-cpp::yaml-cpp)
  else()
    target_link_libraries(PhreeqcRM PUBLIC yaml-cpp)
  endif()
  string(APPEND exportPackageDependencies "find_dependency(yaml-cpp)\n")
  # set(phreeqcrm_requires "yaml-cpp")
  set(useyaml_define "-DUSE_YAML")
endif()
if(DEFINED SKBUILD)
  # SKBUILD is defined when scikit-build-core is used to build python wheels (see pyproject.toml)
  target_compile_definitions(PhreeqcRM PUBLIC USE_YAML)
  target_link_libraries(PhreeqcRM PUBLIC yaml-cpp)
endif()

# check for zlib
if(NOT ZLIB_ROOT)
  if(EXISTS $ENV{ZLIB_ROOT})
    set(ZLIB_ROOT $ENV{ZLIB_ROOT})
  endif()
endif()
option(PHREEQCRM_USE_ZLIB "Use zlib compression for dump files" OFF)
if(DEFINED SKBUILD)
  # SKBUILD is defined when scikit-build-core is used to build python wheels (see pyproject.toml)
  set(PHREEQCRM_USE_ZLIB ON)
endif()
if(PHREEQCRM_USE_ZLIB AND NOT DEFINED SKBUILD)
  find_package(ZLIB REQUIRED)
endif()

if(PHREEQCRM_USE_ZLIB AND ZLIB_FOUND)
  target_compile_definitions(PhreeqcRM PRIVATE USE_GZ)
  target_link_libraries(PhreeqcRM PRIVATE ZLIB::ZLIB)
  string(APPEND exportPackageDependencies "find_dependency(ZLIB)\n")
endif()
if(DEFINED SKBUILD)
  # SKBUILD is defined when scikit-build-core is used to build python wheels (see pyproject.toml)
  target_compile_definitions(PhreeqcRM PRIVATE USE_GZ)
  # for the wheel link to static zlib
  target_link_libraries(PhreeqcRM PRIVATE zlibstatic)
endif()

# Allow user to override POSTFIX properties (but mandate them so that
# all builds can be installed to the same directory)
if(NOT CMAKE_DEBUG_POSTFIX)
  set_target_properties(PhreeqcRM PROPERTIES DEBUG_POSTFIX "d")
endif()

# # At least with MSMPI PhreeqcRM.dll requires MPI_CXX_LIBRARIES to link
# if(MSVC AND BUILD_SHARED_LIBS AND PHREEQCRM_BUILD_MPI)
#   target_link_libraries(PhreeqcRM ${MPI_CXX_LIBRARIES})
# endif()

# if(MSVC AND BUILD_SHARED_LIBS AND PHREEQCRM_USE_ZLIB AND ZLIB_FOUND)
#   target_link_libraries(PhreeqcRM ${ZLIB_LIBRARIES})
# endif()

set(PhreeqcRM_Headers
  ${PhreeqcRM_BINARY_DIR}/irm_dll_export.h
  ${PROJECT_SOURCE_DIR}/src/BMI_interface_C.h
  ${PROJECT_SOURCE_DIR}/src/BMI_interface_F.h
  ${PROJECT_SOURCE_DIR}/src/BMIVariant.h
  ${PROJECT_SOURCE_DIR}/src/bmi.hxx
  ${PROJECT_SOURCE_DIR}/src/BMIPhreeqcRM.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqcCallbacks.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/Parser.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_base.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_exports.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_io.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/phrqtype.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/System.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/Var.h
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqcPhast.h
  ${PROJECT_SOURCE_DIR}/src/IrmResult.h
  ${PROJECT_SOURCE_DIR}/src/PhreeqcRM.h
  ${PROJECT_SOURCE_DIR}/src/RM_interface_C.h
  ${PROJECT_SOURCE_DIR}/src/RMVARS.h
  ${PROJECT_SOURCE_DIR}/src/VarManager.h
  ${PROJECT_SOURCE_DIR}/src/YAML_interface_C.h
  ${PROJECT_SOURCE_DIR}/src/YAML_interface_F.h
  ${PROJECT_SOURCE_DIR}/src/YAMLPhreeqcRM.h
)

set(PhreeqcRM_SRC
  ${PROJECT_SOURCE_DIR}/src/BMI_interface.F90
  ${PROJECT_SOURCE_DIR}/src/BMI_not_implemented.inc
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.F90
  ${PROJECT_SOURCE_DIR}/src/RM_interface.F90
  ${PROJECT_SOURCE_DIR}/src/YAML_interface.F90
)

if(PHREEQCRM_BUILD_MPI)
  # USE_MPI needs to be public in order to compile the
  # fortran interface files
  target_compile_definitions(PhreeqcRM PUBLIC USE_MPI)
  target_link_libraries(PhreeqcRM PUBLIC MPI::MPI_CXX)
  set(usempi_define "-DUSE_MPI")
  string(APPEND exportPackageDependencies "find_dependency(MPI)\n")
else()
  if(NOT PHREEQCRM_DISABLE_OPENMP)
    if(OPENMP_FOUND)
      target_compile_definitions(PhreeqcRM PRIVATE USE_OPENMP)
      target_link_libraries(PhreeqcRM PUBLIC OpenMP::OpenMP_CXX)
      string(APPEND exportPackageDependencies "find_dependency(OpenMP)\n")
      if(PHREEQCRM_FORTRAN_TESTING)
        target_link_libraries(PhreeqcRM PUBLIC OpenMP::OpenMP_Fortran)
      endif()
    endif()
  endif()
endif()

# pkg-config
# Skip imported link libraries
#
# set(PKGCONFIG_LIBS "")
# set(PKGCONFIG_INCS "")
# set(PHREEQCRM_OpenMP_C_FLAGS "")
# if(NOT PHREEQCRM_DISABLE_OPENMP AND OpenMP_FOUND)
#   set(PHREEQCRM_OpenMP_C_FLAGS ${OpenMP_C_FLAGS})
# endif()
# get_target_property(_link_libs PhreeqcRM LINK_LIBRARIES)
# if (_link_libs)
#   foreach(_link_lib ${_link_libs})
#     cmake_print_variables(_link_lib)
#     get_target_property(_imported ${_link_lib} IMPORTED)
#     if(_imported)
#       get_target_property(_incs ${_link_lib} INTERFACE_INCLUDE_DIRECTORIES)
#       if(_incs)
#         foreach(_inc ${_incs})
#           string(APPEND PKGCONFIG_INCS " -I${_inc}")
#         endforeach()
#       endif()
#       get_target_property(_loc ${_link_lib} LOCATION)
#       if (_loc)
#         string(APPEND PKGCONFIG_LIBS " ${_loc}")
#       else()
#         get_target_property(_libs ${_link_lib} INTERFACE_LINK_LIBRARIES)
#         if(_libs)
#           foreach(_lib ${_libs})
#             string(APPEND PKGCONFIG_LIBS " ${_lib}")
#           endforeach()
#         endif()
#       endif()
#     endif()
#   endforeach()
#   if(PKGCONFIG_INCS)
#     string(STRIP ${PKGCONFIG_INCS} PKGCONFIG_INCS)
#   endif()
#   if(PKGCONFIG_LIBS)
#     string(STRIP ${PKGCONFIG_LIBS} PKGCONFIG_LIBS)
#   endif()
# endif()
configure_file(phreeqcrm.cmake.in phreeqcrm.pc @ONLY)

# install
if(STANDALONE_BUILD EQUAL 1 AND NOT DEFINED SKBUILD)
  # SKBUILD is defined when scikit-build-core is used to build python wheels (see pyproject.toml)

  include(CMakePackageConfigHelpers)

  configure_package_config_file(PhreeqcRMConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/PhreeqcRMConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PhreeqcRM
  )

  write_basic_package_version_file("PhreeqcRMConfigVersion.cmake"
    VERSION ${PhreeqcRM_VERSION}
    COMPATIBILITY SameMajorVersion
  )

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/PhreeqcRMConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/PhreeqcRMConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PhreeqcRM
  )

  install(TARGETS PhreeqcRM
    EXPORT   PhreeqcRMTargets
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  install(EXPORT PhreeqcRMTargets
    FILE PhreeqcRMTargets.cmake
    NAMESPACE PhreeqcRM::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PhreeqcRM
  )

  install(FILES ${PhreeqcRM_Headers}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  install(FILES ${PhreeqcRM_SRC}
    # @todo consider installing into INCLUDE directory
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PhreeqcRM/src
  )

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/phreeqcrm.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
  )
else()
  # option(PHREEQCRM_INSTALL_TARGETS "" OFF)
  if(PHREEQCRM_INSTALL_TARGETS)
    install(TARGETS PhreeqcRM DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()
endif()

if(STANDALONE_BUILD EQUAL 1 AND MSVC)
  # PHREEQCRM_STATIC_RUNTIME (VCOMP???.dll) cannot be made static
  option(PHREEQCRM_STATIC_RUNTIME "Build with a static runtime (excluding VCOMP???.DLL)" OFF)
  if(PHREEQCRM_STATIC_RUNTIME)
    # compile with static runtime

    set(CompilerFlags
      CMAKE_CXX_FLAGS
      CMAKE_CXX_FLAGS_DEBUG
      CMAKE_CXX_FLAGS_RELEASE
      CMAKE_CXX_FLAGS_MINSIZEREL
      CMAKE_CXX_FLAGS_RELWITHDEBINFO
      CMAKE_C_FLAGS
      CMAKE_C_FLAGS_DEBUG
      CMAKE_C_FLAGS_RELEASE
      CMAKE_C_FLAGS_MINSIZEREL
      CMAKE_C_FLAGS_RELWITHDEBINFO
    )

    foreach(CompilerFlag ${CompilerFlags})
      if(${CompilerFlag} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
      endif()
    endforeach()

    set(FortranCompilerFlags CMAKE_Fortran_FLAGS)

    foreach(FortranCompilerFlags ${FortranCompilerFlags})
      if(${FortranCompilerFlags} MATCHES "/libs:dll")
        string(REGEX REPLACE "/libs:dll" "/libs:static" ${FortranCompilerFlags} "${${FortranCompilerFlags}}")
      endif()
    endforeach()
  endif()
endif()

if(STANDALONE_BUILD EQUAL 1)
  # tests
  enable_testing()

  # subdirs
  if(NOT DEFINED SKBUILD)
    # SKBUILD is defined when scikit-build-core is used to build python wheels (see pyproject.toml)
    add_subdirectory(database)
    add_subdirectory(doc)
    add_subdirectory(Doxygen)
    add_subdirectory(Tests)
  endif()

  if(PHREEQCRM_BUILD_PYTHON)
    add_subdirectory(swig)
  endif()
endif()

# Turn off Incremental Linking for Debug and RelWithDebInfo
if(MSVC)
  target_link_options(PhreeqcRM PRIVATE $<$<CONFIG:Debug,RelWithDebInfo>:/INCREMENTAL:NO>)
endif()

if(DEFINED SKBUILD)
  # SKBUILD is defined when scikit-build-core is used to build python wheels (see pyproject.toml)
  message(STATUS "SKBUILD IS defined")
  #set_target_properties(PhreeqcRM TestRM TestBMIdtor TestRMdtor PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
  set_target_properties(PhreeqcRM PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
  message(STATUS "SKIP_INSTALL_ALL=${SKIP_INSTALL_ALL}")
  if (SKIP_INSTALL_ALL)
    message(STATUS "SKIP_INSTALL_ALL TRUE")
  else()
    message(STATUS "SKIP_INSTALL_ALL TRUE")
  endif()
else()
  message(STATUS "SKBUILD is NOT defined")
endif()

# get_cmake_property(_variableNames VARIABLES)
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
#   message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()
