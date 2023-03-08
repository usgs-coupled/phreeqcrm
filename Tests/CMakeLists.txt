### tests
if (PHREEQCRM_FORTRAN_TESTING)
  project(FortranAdvect Fortran)

  # TestRM Fortran library
  add_library(FortranAdvect STATIC)

  # TestRM Fortran source
  target_sources(FortranAdvect PRIVATE
    ../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.F90
    ../src/RM_interface.F90
    Advect_f90.F90
    Gas_f90.F90
	  SimpleAdvect_f90.F90
    Species_f90.F90
  )

  # place [*.F[90]] files into "Source Files"
  source_group("Source Files" FILES "../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.F90")
  source_group("Source Files" FILES "../src/RM_interface.F90")
  source_group("Source Files" FILES "Advect_f90.F90")
  source_group("Source Files" FILES "Gas_f90.F90")
  source_group("Source Files" FILES "SimpleAdvect_f90.F90")
  source_group("Source Files" FILES "Species_f90.F90")

  # The next line sets defines/includes/libraries
  # for compiling and linking the FortranAdvect library
  target_link_libraries(FortranAdvect PhreeqcRM::PhreeqcRM)
endif()

project(TestRM CXX)

# Tests/advect.pqi
configure_file(advect.pqi advect.pqi COPYONLY)

# Tests/phreeqc.dat
configure_file(phreeqc.dat phreeqc.dat COPYONLY)

# Tests/units.pqi
configure_file(units.pqi units.pqi COPYONLY)

# Tests/units.pqi
configure_file(gas.pqi gas.pqi COPYONLY)

# test executable
add_executable(TestRM)

# TestRM C/CXX source
target_sources(TestRM PRIVATE
  AdvectBMI_cpp.cpp
  Advect_c.c
  Advect_cpp.cpp
  Gas_c.c
  Gas_cpp.cpp
  main.cpp
  SimpleAdvect_c.c
  SimpleAdvect_cpp.cpp
  Species_c.c
  Species_cpp.cpp
  WriteYAMLFile_cpp.cpp
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

set_property(TARGET TestRM APPEND
  PROPERTY ADDITIONAL_CLEAN_FILES
    $<TARGET_FILE_DIR:TestRM>/AdvectBMI_cpp.chem.txt
    $<TARGET_FILE_DIR:TestRM>/AdvectBMI_cpp.dmp
    $<TARGET_FILE_DIR:TestRM>/AdvectBMI_cpp.log.txt
    $<TARGET_FILE_DIR:TestRM>/AdvectBMI_cpp_utility.txt
    $<TARGET_FILE_DIR:TestRM>/AdvectBMI_cpp.yaml
    $<TARGET_FILE_DIR:TestRM>/Advect_c.chem.txt
    $<TARGET_FILE_DIR:TestRM>/Advect_c.dmp
    $<TARGET_FILE_DIR:TestRM>/Advect_c.log.txt
    $<TARGET_FILE_DIR:TestRM>/Advect_cpp.chem.txt
    $<TARGET_FILE_DIR:TestRM>/Advect_cpp.dmp
    $<TARGET_FILE_DIR:TestRM>/Advect_cpp.log.txt
    $<TARGET_FILE_DIR:TestRM>/Advect_cpp_units_utility.txt
    $<TARGET_FILE_DIR:TestRM>/Advect_cpp_units_worker.chem.txt
    $<TARGET_FILE_DIR:TestRM>/Advect_cpp_units_worker.log.txt
    $<TARGET_FILE_DIR:TestRM>/Advect_cpp_utility.txt
    $<TARGET_FILE_DIR:TestRM>/Advect_c_utility.txt
    $<TARGET_FILE_DIR:TestRM>/Gas_c.chem.txt
    $<TARGET_FILE_DIR:TestRM>/Gas_c.log.txt
    $<TARGET_FILE_DIR:TestRM>/Gas_cpp.chem.txt
    $<TARGET_FILE_DIR:TestRM>/Gas_cpp.log.txt
    $<TARGET_FILE_DIR:TestRM>/SimpleAdvect_c.chem.txt
    $<TARGET_FILE_DIR:TestRM>/SimpleAdvect_c.log.txt
    $<TARGET_FILE_DIR:TestRM>/SimpleAdvect_cpp.chem.txt
    $<TARGET_FILE_DIR:TestRM>/SimpleAdvect_cpp.log.txt
    $<TARGET_FILE_DIR:TestRM>/Species_c.chem.txt
    $<TARGET_FILE_DIR:TestRM>/Species_c.dmp
    $<TARGET_FILE_DIR:TestRM>/Species_c.log.txt
    $<TARGET_FILE_DIR:TestRM>/Species_cpp.chem.txt
    $<TARGET_FILE_DIR:TestRM>/Species_cpp.dmp
    $<TARGET_FILE_DIR:TestRM>/Species_cpp.log.txt
    $<TARGET_FILE_DIR:TestRM>/Species_cpp_utility.txt
    $<TARGET_FILE_DIR:TestRM>/Species_c_utility.txt
)

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
