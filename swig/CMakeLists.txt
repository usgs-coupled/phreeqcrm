include(${SWIG_USE_FILE})

if (PHREEQCRM_BUILD_PYTHON OR DEFINED SKBUILD)
  add_subdirectory(python)
endif()
