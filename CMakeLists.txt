cmake_minimum_required (VERSION 2.8)

##include(CheckLanguage)
##check_language(Fortran)
##if (CMAKE_Fortran_COMPILER)
##  enable_language(Fortran)
##endif()

if (POLICY CMP0054)
  cmake_policy(PUSH)
  cmake_policy(SET CMP0054 OLD)
endif()
project (PhreeqcRM CXX C )
if (POLICY CMP0054)
  cmake_policy(POP)
endif()

# overide docdir on windows
if(WIN32 AND NOT CMAKE_INSTALL_DOCDIR)
  set(CMAKE_INSTALL_DOCDIR "" CACHE PATH "documentation root (doc)")
  set(CMAKE_INSTALL_DOCDIR "doc")
endif()

### Fortran tests option
##if (NOT PHREEQCRM_FORTRAN_TESTING)
##  if (NOT CMAKE_Fortran_COMPILER OR CMAKE_Fortran_COMPILER STREQUAL "CMAKE_Fortran_COMPILER-NOTFOUND")
##    option (PHREEQCRM_FORTRAN_TESTING "Build FORTRAN support" OFF)
##  else()
##    option (PHREEQCRM_FORTRAN_TESTING "Build FORTRAN support" ON)
##  endif()
##else()
##  if (PHREEQCRM_FORTRAN_TESTING)
##    enable_language(Fortran)
##    ##include(FortranCInterface)
##    ##FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
##    ##FortranCInterface_VERIFY(CXX QUIET)
##  endif()  
##endif()


option(PHREEQCRM_FORTRAN_TESTING "Build Fortran test" OFF)
if (PHREEQCRM_FORTRAN_TESTING)
  enable_language(Fortran)
endif()

include(GNUInstallDirs)

# Disable deprecation warnings for standard C and STL functions in VS2005
# and later
IF(MSVC_VERSION EQUAL 1400 OR MSVC_VERSION GREATER 1400)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
ENDIF(MSVC_VERSION EQUAL 1400 OR MSVC_VERSION GREATER 1400) 


### Copy Release to openmp
##SET(CMAKE_CXX_FLAGS_OPENMP
##    ${CMAKE_CXX_FLAGS_RELEASE} CACHE STRING
##    "Flags used by the C++ compiler during openmp builds."
##    FORCE)
##SET(CMAKE_C_FLAGS_OPENMP
##    ${CMAKE_C_FLAGS_RELEASE} CACHE STRING
##    "Flags used by the C compiler during openmp builds."
##    FORCE)
##if (PHREEQCRM_FORTRAN_TESTING)
##  SET(CMAKE_Fortran_FLAGS_OPENMP
##      ${CMAKE_Fortran_FLAGS_RELEASE} CACHE STRING
##      "Flags used by the Fortran compiler during openmp builds."
##      FORCE)
##endif ()
##SET(CMAKE_EXE_LINKER_FLAGS_OPENMP
##    ${CMAKE_EXE_LINKER_FLAGS_RELEASE} CACHE STRING
##    "Flags used for linking binaries during openmp builds."
##    FORCE)
##SET(CMAKE_SHARED_LINKER_FLAGS_OPENMP
##    ${CMAKE_SHARED_LINKER_FLAGS_RELEASE} CACHE STRING
##    "Flags used by the shared libraries linker during openmp builds."
##    FORCE)
##MARK_AS_ADVANCED(
##    CMAKE_CXX_FLAGS_OPENMP
##    CMAKE_C_FLAGS_OPENMP
##    CMAKE_EXE_LINKER_FLAGS_OPENMP
##    CMAKE_SHARED_LINKER_FLAGS_OPENMP)
##if (PHREEQCRM_FORTRAN_TESTING)
##  MARK_AS_ADVANCED(CMAKE_Fortran_FLAGS_OPENMP)
##endif()

    
### Copy Debug to openmp_debug
##SET(CMAKE_CXX_FLAGS_OPENMP_DEBUG
##    ${CMAKE_CXX_FLAGS_DEBUG} CACHE STRING
##    "Flags used by the C++ compiler during debug openmp builds."
##    FORCE)
##SET(CMAKE_C_FLAGS_OPENMP_DEBUG
##    ${CMAKE_C_FLAGS_DEBUG} CACHE STRING
##    "Flags used by the C compiler during debug openmp builds."
##    FORCE)
##if (PHREEQCRM_FORTRAN_TESTING)
##  SET(CMAKE_Fortran_FLAGS_OPENMP_DEBUG
##      ${CMAKE_Fortran_FLAGS_DEBUG} CACHE STRING
##      "Flags used by the Fortran compiler during debug openmp builds."
##      FORCE)
##endif ()
##SET(CMAKE_EXE_LINKER_FLAGS_OPENMP_DEBUG
##    ${CMAKE_EXE_LINKER_FLAGS_DEBUG} CACHE STRING
##    "Flags used for linking binaries during debug openmp builds."
##    FORCE)
##SET(CMAKE_SHARED_LINKER_FLAGS_OPENMP_DEBUG
##    ${CMAKE_SHARED_LINKER_FLAGS_DEBUG} CACHE STRING
##    "Flags used by the shared libraries linker during debug openmp builds."
##    FORCE)
##MARK_AS_ADVANCED(
##    CMAKE_CXX_FLAGS_OPENMP_DEBUG
##    CMAKE_C_FLAGS_OPENMP_DEBUG
##    CMAKE_EXE_LINKER_FLAGS_OPENMP_DEBUG
##    CMAKE_SHARED_LINKER_FLAGS_OPENMP_DEBUG)
##if (PHREEQCRM_FORTRAN_TESTING)
##  MARK_AS_ADVANCED(CMAKE_Fortran_FLAGS_OPENMP_DEBUG)
##endif()


### Copy Release to mpi
##SET(CMAKE_CXX_FLAGS_MPI
##    ${CMAKE_CXX_FLAGS_RELEASE} CACHE STRING
##    "Flags used by the C++ compiler during mpi builds."
##    FORCE)
##SET(CMAKE_C_FLAGS_MPI
##    ${CMAKE_C_FLAGS_RELEASE} CACHE STRING
##    "Flags used by the C compiler during mpi builds."
##    FORCE)
##if (PHREEQCRM_FORTRAN_TESTING)
##  SET(CMAKE_Fortran_FLAGS_MPI
##      ${CMAKE_Fortran_FLAGS_RELEASE} CACHE STRING
##      "Flags used by the Fortran compiler during mpi builds."
##      FORCE)
##endif ()
##SET(CMAKE_EXE_LINKER_FLAGS_MPI
##    ${CMAKE_EXE_LINKER_FLAGS_RELEASE} CACHE STRING
##    "Flags used for linking binaries during mpi builds."
##    FORCE)
##SET(CMAKE_SHARED_LINKER_FLAGS_MPI
##    ${CMAKE_SHARED_LINKER_FLAGS_RELEASE} CACHE STRING
##    "Flags used by the shared libraries linker during mpi builds."
##    FORCE)
##MARK_AS_ADVANCED(
##    CMAKE_CXX_FLAGS_MPI
##    CMAKE_C_FLAGS_MPI
##    CMAKE_EXE_LINKER_FLAGS_MPI
##    CMAKE_SHARED_LINKER_FLAGS_MPI)
##if (PHREEQCRM_FORTRAN_TESTING)
##  MARK_AS_ADVANCED(CMAKE_Fortran_FLAGS_MPI)
##endif()


### Copy Debug to mpi_debug
##SET(CMAKE_CXX_FLAGS_MPI_DEBUG
##    ${CMAKE_CXX_FLAGS_DEBUG} CACHE STRING
##    "Flags used by the C++ compiler during debug mpi builds."
##    FORCE)
##SET(CMAKE_C_FLAGS_MPI_DEBUG
##    ${CMAKE_C_FLAGS_DEBUG} CACHE STRING
##    "Flags used by the C compiler during debug mpi builds."
##    FORCE)
##if (PHREEQCRM_FORTRAN_TESTING)
##  SET(CMAKE_Fortran_FLAGS_MPI_DEBUG
##      ${CMAKE_Fortran_FLAGS_DEBUG} CACHE STRING
##      "Flags used by the Fortran compiler during debug mpi builds."
##      FORCE)
##endif ()
##SET(CMAKE_EXE_LINKER_FLAGS_MPI_DEBUG
##    ${CMAKE_EXE_LINKER_FLAGS_DEBUG} CACHE STRING
##    "Flags used for linking binaries during debug mpi builds."
##    FORCE)
##SET(CMAKE_SHARED_LINKER_FLAGS_MPI_DEBUG
##    ${CMAKE_SHARED_LINKER_FLAGS_DEBUG} CACHE STRING
##    "Flags used by the shared libraries linker during debug mpi builds."
##    FORCE)
##MARK_AS_ADVANCED(
##    CMAKE_CXX_FLAGS_MPI_DEBUG
##    CMAKE_C_FLAGS_MPI_DEBUG
##    CMAKE_EXE_LINKER_FLAGS_MPI_DEBUG
##    CMAKE_SHARED_LINKER_FLAGS_MPI_DEBUG)
##if (PHREEQCRM_FORTRAN_TESTING)
##  MARK_AS_ADVANCED(CMAKE_Fortran_FLAGS_MPI_DEBUG)
##endif()


# check for openmp
find_package(OpenMP)
if (OPENMP_FOUND)
##  # openmp release
##  SET(CMAKE_C_FLAGS_OPENMP
##      "-DUSE_OPENMP ${CMAKE_C_FLAGS_OPENMP} ${OpenMP_C_FLAGS}" CACHE STRING
##      "Flags used by the compiler during openmp builds."
##      FORCE)
##  SET(CMAKE_CXX_FLAGS_OPENMP
##      "-DUSE_OPENMP ${CMAKE_CXX_FLAGS_OPENMP} ${OpenMP_CXX_FLAGS}" CACHE STRING
##      "Flags used by the compiler during openmp builds."
##      FORCE)
##  if (PHREEQCRM_FORTRAN_TESTING)
##    SET(CMAKE_Fortran_FLAGS_OPENMP
##        "-DUSE_OPENMP ${CMAKE_Fortran_FLAGS_OPENMP} ${OpenMP_Fortran_FLAGS}" CACHE STRING
##        "Flags used by the compiler during openmp builds."
##        FORCE)
##  endif ()
##
##  # openmp debug
##  SET(CMAKE_C_FLAGS_OPENMP_DEBUG
##      "-DUSE_OPENMP ${CMAKE_C_FLAGS_OPENMP_DEBUG} ${OpenMP_C_FLAGS}" CACHE STRING
##      "Flags used by the compiler during debug openmp builds."
##      FORCE)
##  SET(CMAKE_CXX_FLAGS_OPENMP_DEBUG
##      "-DUSE_OPENMP ${CMAKE_CXX_FLAGS_OPENMP_DEBUG} ${OpenMP_CXX_FLAGS}" CACHE STRING
##      "Flags used by the compiler during debug openmp builds."
##      FORCE)
##  if (PHREEQCRM_FORTRAN_TESTING)
##    SET(CMAKE_Fortran_FLAGS_OPENMP_DEBUG
##        "-DUSE_OPENMP ${CMAKE_Fortran_FLAGS_OPENMP_DEBUG} ${OpenMP_Fortran_FLAGS}" CACHE STRING
##        "Flags used by the compiler during debug openmp builds."
##        FORCE)
##  endif ()
endif (OPENMP_FOUND)


# check for MPI
option(PHREEQCRM_BUILD_MPI "Build MPI configuration" "OFF")
if (PHREEQCRM_BUILD_MPI)
  cmake_policy(SET CMP0054 OLD)
  find_package(MPI)
  cmake_policy(SET CMP0054 NEW)
  
##  if (MPI_CXX_FOUND)
##    # mpi release
##    SET(CMAKE_C_FLAGS_MPI
##        "-DUSE_MPI ${CMAKE_C_FLAGS_MPI} ${MPI_C_COMPILE_FLAGS}" CACHE STRING
##        "Flags used by the compiler during mpi builds."
##        FORCE)
##    SET(CMAKE_CXX_FLAGS_MPI
##        "-DUSE_MPI ${CMAKE_CXX_FLAGS_MPI} ${MPI_CXX_COMPILE_FLAGS}" CACHE STRING
##        "Flags used by the compiler during mpi builds."
##        FORCE)
##    if (PHREEQCRM_FORTRAN_TESTING)
##      SET(CMAKE_Fortran_FLAGS_MPI
##          "-DUSE_MPI ${CMAKE_Fortran_FLAGS_MPI} ${MPI_Fortran_COMPILE_FLAGS}" CACHE STRING
##          "Flags used by the Fortran compiler during mpi builds."
##          FORCE)
##    endif ()
##    if (MPI_C_INCLUDE_PATH)
##      SET(CMAKE_C_FLAGS_MPI
##          "-I\"${MPI_C_INCLUDE_PATH}\" ${CMAKE_C_FLAGS_MPI}" CACHE STRING
##          "Flags used by the compiler during mpi builds."
##          FORCE)
##    endif (MPI_C_INCLUDE_PATH)
##    if (MPI_CXX_INCLUDE_PATH)
##      SET(CMAKE_CXX_FLAGS_MPI
##          "-I\"${MPI_CXX_INCLUDE_PATH}\" ${CMAKE_CXX_FLAGS_MPI}" CACHE STRING
##          "Flags used by the compiler during mpi builds."
##          FORCE)
##    endif (MPI_CXX_INCLUDE_PATH)
##    if (PHREEQCRM_FORTRAN_TESTING AND MPI_Fortran_INCLUDE_PATH)
##      if("${CMAKE_SIZEOF_VOID_P}" EQUAL 8)
##        set(MS_MPI_ARCH_DIR amd64)
##      else()
##        set(MS_MPI_ARCH_DIR i386)
##      endif()
##      SET(CMAKE_Fortran_FLAGS_MPI
##          "-I\"${MPI_Fortran_INCLUDE_PATH}\" -I\"${MPI_Fortran_INCLUDE_PATH}/${MS_MPI_ARCH_DIR}\" ${CMAKE_Fortran_FLAGS_MPI}" CACHE STRING
##          "Flags used by the Fortran compiler during mpi builds."
##          FORCE)
##    endif ()    
##    
##    # mpi debug
##    SET(CMAKE_C_FLAGS_MPI_DEBUG
##        "-DUSE_MPI ${CMAKE_C_FLAGS_MPI_DEBUG} ${MPI_C_COMPILE_FLAGS}" CACHE STRING
##        "Flags used by the compiler during debug mpi builds."
##        FORCE)
##    SET(CMAKE_CXX_FLAGS_MPI_DEBUG
##        "-DUSE_MPI ${CMAKE_CXX_FLAGS_MPI_DEBUG} ${MPI_CXX_COMPILE_FLAGS}" CACHE STRING
##        "Flags used by the compiler during debug mpi builds."
##        FORCE)
##    if (PHREEQCRM_FORTRAN_TESTING)
##      SET(CMAKE_Fortran_FLAGS_MPI_DEBUG
##          "-DUSE_MPI ${CMAKE_Fortran_FLAGS_MPI_DEBUG} ${MPI_Fortran_COMPILE_FLAGS}" CACHE STRING
##          "Flags used by the Fortran compiler during debug mpi builds."
##          FORCE)
##    endif ()
##    if (MPI_C_INCLUDE_PATH)
##      SET(CMAKE_C_FLAGS_MPI_DEBUG
##          "-I\"${MPI_C_INCLUDE_PATH}\" ${CMAKE_C_FLAGS_MPI_DEBUG}" CACHE STRING
##          "Flags used by the compiler during debug mpi builds."
##          FORCE)
##    endif (MPI_C_INCLUDE_PATH)
##    if (MPI_CXX_INCLUDE_PATH)
##      SET(CMAKE_CXX_FLAGS_MPI_DEBUG
##          "-I\"${MPI_CXX_INCLUDE_PATH}\" ${CMAKE_CXX_FLAGS_MPI_DEBUG}" CACHE STRING
##          "Flags used by the compiler during debug mpi builds."
##          FORCE)
##    endif (MPI_CXX_INCLUDE_PATH)    
##    if (PHREEQCRM_FORTRAN_TESTING AND MPI_Fortran_INCLUDE_PATH)
##      if("${CMAKE_SIZEOF_VOID_P}" EQUAL 8)
##        set(MS_MPI_ARCH_DIR amd64)
##      else()
##        set(MS_MPI_ARCH_DIR i386)
##      endif()
##      SET(CMAKE_Fortran_FLAGS_MPI_DEBUG
##          "-I\"${MPI_Fortran_INCLUDE_PATH}\" -I\"${MPI_Fortran_INCLUDE_PATH}/${MS_MPI_ARCH_DIR}\" ${CMAKE_Fortran_FLAGS_MPI_DEBUG}" CACHE STRING
##          "Flags used by the Fortran compiler during debug mpi builds."
##          FORCE)
##    endif ()    
##  endif (MPI_CXX_FOUND)
endif (PHREEQCRM_BUILD_MPI)


##if (PHREEQCRM_BUILD_MPI)
##  if (CMAKE_CONFIGURATION_TYPES)
##    SET(CMAKE_CONFIGURATION_TYPES "openmp_debug;openmp;mpi;mpi_debug" CACHE STRING 
##    "Semicolon separated list of supported configuration types."
##    FORCE)
##  endif (CMAKE_CONFIGURATION_TYPES)
##else (PHREEQCRM_BUILD_MPI)
##  if (CMAKE_CONFIGURATION_TYPES)
##    SET(CMAKE_CONFIGURATION_TYPES "openmp_debug;openmp" CACHE STRING 
##    "Semicolon separated list of supported configuration types."
##    FORCE)
##  endif (CMAKE_CONFIGURATION_TYPES)
##endif (PHREEQCRM_BUILD_MPI)


SET(PhreeqcRM_SOURCES
src/IrmResult.h
src/PhreeqcRM.cpp
src/PhreeqcRM.h
src/RM_interface_C.cpp
src/RM_interface_C.h
src/RM_interface_F.cpp
src/RM_interface_F.h
src/IPhreeqcPhast/IPhreeqcPhast.cxx
src/IPhreeqcPhast/IPhreeqcPhast.h
src/IPhreeqcPhast/IPhreeqcPhastLib.cpp
src/IPhreeqcPhast/IPhreeqcPhastLib.h
src/IPhreeqcPhast/IPhreeqc/thread.h
src/IPhreeqcPhast/IPhreeqc/CSelectedOutput.cpp
src/IPhreeqcPhast/IPhreeqc/CSelectedOutput.hxx
src/IPhreeqcPhast/IPhreeqc/CVar.hxx
src/IPhreeqcPhast/IPhreeqc/Debug.h
src/IPhreeqcPhast/IPhreeqc/ErrorReporter.hxx
src/IPhreeqcPhast/IPhreeqc/IPhreeqc.cpp
src/IPhreeqcPhast/IPhreeqc/IPhreeqc.h
src/IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp
src/IPhreeqcPhast/IPhreeqc/IPhreeqcCallbacks.h
src/IPhreeqcPhast/IPhreeqc/IPhreeqcLib.cpp
src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface_F.cpp
src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface_F.h
src/IPhreeqcPhast/IPhreeqc/Var.c
src/IPhreeqcPhast/IPhreeqc/Var.h
src/IPhreeqcPhast/IPhreeqc/Version.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartHandler.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartHandler.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartObject.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartObject.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/CurveObject.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/CurveObject.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ExchComp.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ExchComp.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Exchange.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Exchange.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasComp.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasComp.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasPhase.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasPhase.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolution.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolution.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolutionComp.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolutionComp.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/KineticsComp.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/KineticsComp.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NA.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NumKeyword.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NumKeyword.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PBasic.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PBasic.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_base.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_base.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io_output.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblage.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblage.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblageComp.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblageComp.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Parser.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Parser.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Phreeqc.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Phreeqc.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Pressure.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Pressure.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Reaction.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Reaction.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ReadClass.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SS.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SS.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SSassemblage.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SSassemblage.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SScomp.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SScomp.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SelectedOutput.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SelectedOutput.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Solution.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Solution.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SolutionIsotope.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SolutionIsotope.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBinList.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBinList.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Surface.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Surface.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceCharge.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceCharge.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceComp.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceComp.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/System.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/System.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Temperature.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Temperature.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Use.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Use.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/UserPunch.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/UserPunch.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Utils.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Utils.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/advection.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/basicsubs.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cl1.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cvdense.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cvdense.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cvode.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cvode.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxKinetics.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxKinetics.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxMix.cxx
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxMix.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dense.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dense.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dumper.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dumper.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dw.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/gases.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/global_structures.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/input.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/integrate.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/inverse.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/isotopes.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/kinetics.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/mainsubs.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/model.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector_serial.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector_serial.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/parse.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/phqalloc.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/phqalloc.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/phrqtype.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/pitzer.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/pitzer_structures.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/prep.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/print.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/read.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/readtr.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/runner.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/runner.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sit.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/smalldense.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/smalldense.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/spread.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/step.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/structures.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialsmath.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialsmath.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialstypes.h
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/tally.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/tidy.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/transport.cpp
src/IPhreeqcPhast/IPhreeqc/phreeqcpp/utilities.cpp
)

##if (PHREEQCRM_FORTRAN_TESTING AND MSVC)
##  SET(PhreeqcRM_SOURCES
##  ${PhreeqcRM_SOURCES}
##  src/IPhreeqcPhast/IPhreeqc/fwrap.cpp
##  src/IPhreeqcPhast/IPhreeqc/fwrap1.cpp
##  src/IPhreeqcPhast/IPhreeqc/fwrap2.cpp
##  src/IPhreeqcPhast/IPhreeqc/fwrap3.cpp
##  src/IPhreeqcPhast/IPhreeqc/fwrap4.cpp
##  src/IPhreeqcPhast/IPhreeqc/fwrap5.cpp
##  src/IPhreeqcPhast/IPhreeqc/fwrap6.cpp      
##  src/IPhreeqcPhast/IPhreeqc/fwrap7.cpp
##  src/IPhreeqcPhast/IPhreeqc/fwrap8.cpp      
##  )
##endif ()

# setup include directories
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast")
include_directories("${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc")
include_directories("${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp")
include_directories("${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords")

# dll/shared objects
OPTION(BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
SET(LIB_TYPE STATIC)
if (BUILD_SHARED_LIBS)
  SET(LIB_TYPE SHARED)
endif()
add_library(PhreeqcRM ${LIB_TYPE} ${PhreeqcRM_SOURCES})
set_target_properties(PhreeqcRM PROPERTIES DEBUG_POSTFIX "d")

# At least with MSMPI PhreeqcRM.dll requires MPI_CXX_LIBRARIES to link
if (MSVC AND BUILD_SHARED_LIBS AND PHREEQCRM_BUILD_MPI)
  target_link_libraries(PhreeqcRM ${MPI_CXX_LIBRARIES})
endif()

SET(PhreeqcRM_Headers
${PROJECT_SOURCE_DIR}/src/IrmResult.h
${PROJECT_SOURCE_DIR}/src/PhreeqcRM.h
${PROJECT_SOURCE_DIR}/src/RM_interface_C.h
${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqcPhast.h
${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.h
${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp
${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqcCallbacks.h
${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/Var.h
)

if (PHREEQCRM_FORTRAN_TESTING)
  SET(PhreeqcRM_Headers
  ${PhreeqcRM_Headers}
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.f90.inc
  ${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.f.inc
  )
endif()

SET(PhreeqcRM_SRC
${PROJECT_SOURCE_DIR}/src/RM_interface.F90
${PROJECT_SOURCE_DIR}/src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.F90
)

SET(PhreeqcRM_ROOT
${PROJECT_SOURCE_DIR}/doc/NOTICE.TXT
${PROJECT_SOURCE_DIR}/doc/README.TXT
)

if (PHREEQCRM_BUILD_MPI)
  ADD_DEFINITIONS(-DUSE_MPI)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  if (PHREEQCRM_FORTRAN_TESTING AND MPI_Fortran_INCLUDE_PATH)
    if("${CMAKE_SIZEOF_VOID_P}" EQUAL 8)
      set(MS_MPI_ARCH_DIR amd64)
    else()
      set(MS_MPI_ARCH_DIR i386)
    endif()
    include_directories(${MPI_Fortran_INCLUDE_PATH}/${MS_MPI_ARCH_DIR})
  endif()
  set(LINK_LIBS ${LINK_LIBS} ${MPI_CXX_LIBRARIES})
  if (PHREEQCRM_FORTRAN_TESTING AND MPI_Fortran_LIBRARIES MATCHES "msmpi.lib$")
    STRING(REGEX REPLACE "msmpi.lib$" "msmpifec.lib" EXTRA_FORTRAN_LIB ${MPI_Fortran_LIBRARIES})
    SET(LINK_LIBS ${LINK_LIBS} ${EXTRA_FORTRAN_LIB})
  endif ()  
else()
  if (OPENMP_FOUND)
    ADD_DEFINITIONS(-DUSE_OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")
    if (PHREEQCRM_FORTRAN_TESTING)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    endif()
  endif()
endif()  

# install
install(TARGETS PhreeqcRM DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${PhreeqcRM_Headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${PhreeqcRM_SRC} DESTINATION src)
install(FILES ${PhreeqcRM_ROOT} DESTINATION .)

# tests
enable_testing()

# subdirs
add_subdirectory(Doxygen)
add_subdirectory(Tests)
