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
IrmResult.h
PhreeqcRM.cpp
PhreeqcRM.h
RM_interface_C.cpp
RM_interface_C.h
RM_interface_F.cpp
RM_interface_F.h
IPhreeqcPhast/IPhreeqcPhast.cxx
IPhreeqcPhast/IPhreeqcPhast.h
IPhreeqcPhast/IPhreeqcPhastLib.cpp
IPhreeqcPhast/IPhreeqcPhastLib.h
IPhreeqcPhast/IPhreeqc/thread.h
IPhreeqcPhast/IPhreeqc/CSelectedOutput.cpp
IPhreeqcPhast/IPhreeqc/CSelectedOutput.hxx
IPhreeqcPhast/IPhreeqc/CVar.hxx
IPhreeqcPhast/IPhreeqc/Debug.h
IPhreeqcPhast/IPhreeqc/ErrorReporter.hxx
IPhreeqcPhast/IPhreeqc/IPhreeqc.cpp
IPhreeqcPhast/IPhreeqc/IPhreeqc.h
IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp
IPhreeqcPhast/IPhreeqc/IPhreeqcCallbacks.h
IPhreeqcPhast/IPhreeqc/IPhreeqcLib.cpp
IPhreeqcPhast/IPhreeqc/IPhreeqc_interface_F.cpp
IPhreeqcPhast/IPhreeqc/IPhreeqc_interface_F.h
IPhreeqcPhast/IPhreeqc/Var.c
IPhreeqcPhast/IPhreeqc/Var.h
IPhreeqcPhast/IPhreeqc/Version.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartHandler.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartHandler.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartObject.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartObject.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/CurveObject.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/CurveObject.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/ExchComp.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/ExchComp.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/Exchange.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/Exchange.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/GasComp.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/GasComp.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/GasPhase.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/GasPhase.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolution.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolution.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolutionComp.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolutionComp.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/KineticsComp.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/KineticsComp.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/NA.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/NumKeyword.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/NumKeyword.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/PBasic.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/PBasic.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_base.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_base.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io_output.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblage.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblage.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblageComp.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblageComp.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/Parser.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/Parser.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/Phreeqc.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/Phreeqc.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/Pressure.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/Pressure.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/Reaction.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/Reaction.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/ReadClass.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/SS.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/SS.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/SSassemblage.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/SSassemblage.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/SScomp.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/SScomp.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/SelectedOutput.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/SelectedOutput.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/Solution.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/Solution.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/SolutionIsotope.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/SolutionIsotope.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBinList.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBinList.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/Surface.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/Surface.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceCharge.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceCharge.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceComp.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceComp.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/System.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/System.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/Temperature.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/Temperature.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/Use.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/Use.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/UserPunch.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/UserPunch.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/Utils.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/Utils.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/advection.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/basicsubs.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/cl1.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/cvdense.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/cvdense.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/cvode.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/cvode.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxKinetics.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxKinetics.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxMix.cxx
IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxMix.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/dense.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/dense.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/dumper.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/dumper.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/dw.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/gases.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/global_structures.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/input.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/integrate.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/inverse.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/isotopes.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/kinetics.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/mainsubs.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/model.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector_serial.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector_serial.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/parse.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/phqalloc.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/phqalloc.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/phrqtype.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/pitzer.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/pitzer_structures.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/prep.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/print.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/read.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/readtr.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/runner.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/runner.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/sit.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/smalldense.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/smalldense.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/spread.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/step.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/structures.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialsmath.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialsmath.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialstypes.h
IPhreeqcPhast/IPhreeqc/phreeqcpp/tally.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/tidy.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/transport.cpp
IPhreeqcPhast/IPhreeqc/phreeqcpp/utilities.cpp
)

##if (PHREEQCRM_FORTRAN_TESTING AND MSVC)
##  SET(PhreeqcRM_SOURCES
##  ${PhreeqcRM_SOURCES}
##  IPhreeqcPhast/IPhreeqc/fwrap.cpp
##  IPhreeqcPhast/IPhreeqc/fwrap1.cpp
##  IPhreeqcPhast/IPhreeqc/fwrap2.cpp
##  IPhreeqcPhast/IPhreeqc/fwrap3.cpp
##  IPhreeqcPhast/IPhreeqc/fwrap4.cpp
##  IPhreeqcPhast/IPhreeqc/fwrap5.cpp
##  IPhreeqcPhast/IPhreeqc/fwrap6.cpp      
##  IPhreeqcPhast/IPhreeqc/fwrap7.cpp
##  IPhreeqcPhast/IPhreeqc/fwrap8.cpp      
##  )
##endif ()

# setup include directories
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/IPhreeqcPhast")
include_directories("${PROJECT_SOURCE_DIR}/IPhreeqcPhast/IPhreeqc")
include_directories("${PROJECT_SOURCE_DIR}/IPhreeqcPhast/IPhreeqc/phreeqcpp")
include_directories("${PROJECT_SOURCE_DIR}/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords")

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
${PROJECT_SOURCE_DIR}/IrmResult.h
${PROJECT_SOURCE_DIR}/PhreeqcRM.h
${PROJECT_SOURCE_DIR}/RM_interface_C.h
${PROJECT_SOURCE_DIR}/IPhreeqcPhast/IPhreeqcPhast.h
${PROJECT_SOURCE_DIR}/IPhreeqcPhast/IPhreeqc/IPhreeqc.h
${PROJECT_SOURCE_DIR}/IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp
${PROJECT_SOURCE_DIR}/IPhreeqcPhast/IPhreeqc/IPhreeqcCallbacks.h
${PROJECT_SOURCE_DIR}/IPhreeqcPhast/IPhreeqc/Var.h
)

if (PHREEQCRM_FORTRAN_TESTING)
  SET(PhreeqcRM_Headers
  ${PhreeqcRM_Headers}
  ${PROJECT_SOURCE_DIR}/IPhreeqcPhast/IPhreeqc/IPhreeqc.f90.inc
  ${PROJECT_SOURCE_DIR}/IPhreeqcPhast/IPhreeqc/IPhreeqc.f.inc
  )
endif()

SET(PhreeqcRM_SRC
${PROJECT_SOURCE_DIR}/RM_interface.F90
${PROJECT_SOURCE_DIR}/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.F90
)

##{{
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
##}}

# install
install(TARGETS PhreeqcRM DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${PhreeqcRM_Headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${PhreeqcRM_SRC} DESTINATION src)

# tests
enable_testing()

add_subdirectory(Tests)