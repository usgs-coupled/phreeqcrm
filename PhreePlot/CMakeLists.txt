project(PhreePlot Fortran)

configure_file (
  ../database/wateq4f.dat
  wateq4f.dat
  COPYONLY
  )

configure_file (
  preloop2.pqi
  preloop2.pqi
  COPYONLY
  )

# source
SET(PhreePlot_SOURCES
  PhreePlotOpenMP.f90
  ../src/RM_interface.F90
  ../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.F90
)

add_executable(PhreePlot ${PhreePlot_SOURCES})

# library dependencies
if (PHREEQCRM_BUILD_MPI)
  SET(EXTRA_LIBS ${EXTRA_LIBS} PhreeqcRM ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})
else ()
  SET(EXTRA_LIBS ${EXTRA_LIBS} PhreeqcRM)
endif ()

# link (zlib must be last on linux)
if (PHREEQCRM_USE_ZLIB AND ZLIB_FOUND)
  set(EXTRA_LIBS ${EXTRA_LIBS} ${ZLIB_LIBRARIES})
endif()
target_link_libraries(PhreePlot ${EXTRA_LIBS})
