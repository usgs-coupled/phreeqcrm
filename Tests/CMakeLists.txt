### tests
if (PHREEQCRM_FORTRAN_TESTING)
  project(FortranAdvect Fortran)
  
  # TestRM Fortran library
  add_library(FortranAdvect STATIC)

  # TestRM Fortran source
  target_sources(FortranAdvect PRIVATE
    ../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.F90
    ../src/RM_interface.F90
    advection_f90.F90
    gas_f90.F90
    species_f90.F90
  )

  # place [*.F[90]] files into "Source Files"
  source_group("Source Files" FILES "../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.F90")
  source_group("Source Files" FILES "../src/RM_interface.F90")
  source_group("Source Files" FILES "advection_f90.F90")
  source_group("Source Files" FILES "gas_f90.F90")
  source_group("Source Files" FILES "species_f90.F90")
  
  # The next line sets defines/includes/libraries
  # for compiling and linking the FortranAdvect library
  target_link_libraries(FortranAdvect PhreeqcRM::PhreeqcRM)
endif()

project(TestRM CXX)

# Tests/advect.pqi
configure_file(advect.pqi.in advect.pqi COPYONLY)

# Tests/phreeqc.dat
configure_file(phreeqc.dat.in phreeqc.dat COPYONLY)
 
# Tests/units.pqi
configure_file(units.pqi.in units.pqi COPYONLY)
  
# Tests/units.pqi
configure_file(gas.pqi.in gas.pqi COPYONLY)

# test executable
add_executable(TestRM)

# TestRM C/CXX source
target_sources(TestRM PRIVATE
  advection_bmi_cpp.cpp
  advection_c.c
  advection_cpp.cpp
  gas_c.c
  gas_cpp.cpp
  main.cpp
  SimpleAdvection.cpp
  species_c.c
  species_cpp.cpp
  WriteYAMLFile.cpp
  YAMLPhreeqcRM.cpp
  YAMLPhreeqcRM.h
)

if(PHREEQCRM_FORTRAN_TESTING)
  target_compile_definitions(TestRM PRIVATE TEST_FORTRAN)
endif()

# library dependencies
if(PHREEQCRM_FORTRAN_TESTING)
  # Override bug reported on:
  # http://www.cmake.org/pipermail/cmake/2009-July/030954.html
  if(WIN32 AND ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")  
    if(NOT CMAKE_Fortran_FLAGS_DEBUG MATCHES "/Od")
      set(CMAKE_Fortran_FLAGS_DEBUG "/Od ${CMAKE_Fortran_FLAGS_DEBUG}")
    endif()
  endif()

  target_link_libraries(TestRM PhreeqcRM::PhreeqcRM FortranAdvect)
else()
  target_link_libraries(TestRM PhreeqcRM::PhreeqcRM)
endif()

# test compile and run
add_test(TestCompileAndRun TestRM)

if(MSVC AND BUILD_SHARED_LIBS)
  # copy dll
  add_custom_command(TARGET TestRM POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PhreeqcRM> $<TARGET_FILE_DIR:TestRM>
  )
endif()

# MPI test
# get_property(exe_location TARGET testrm PROPERTY LOCATION)
# message (STATUS "exe_location = ${exe_location}")
if(PHREEQCRM_BUILD_MPI)
  add_test(NAME TestMPI
           COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3
           ${MPIEXEC_PREFLAGS}
           $<TARGET_FILE:TestRM>
           ${MPIEXEC_POSTFLAGS}
           )
endif()
