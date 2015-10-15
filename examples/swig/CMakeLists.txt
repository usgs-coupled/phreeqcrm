if (SWIG_ENABLE_JAVA)
  add_subdirectory(java)
endif()

if (SWIG_ENABLE_PYTHON)
  add_subdirectory(python)
endif()

if (SWIG_ENABLE_RUBY)
  add_subdirectory(ruby)
endif()
