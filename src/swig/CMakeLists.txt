include(${SWIG_USE_FILE})

if (PHREEQCRM_BUILD_PYTHON)
  add_subdirectory(python)
endif()
