set(PYTHON_TARGET_NAME phreeqcrmpy)

# copy database to build directory for testing
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../../database/phreeqc.dat phreeqc.dat COPYONLY)

# copy advect.pqi to build directory for testing
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../../Tests/advect.pqi advect.pqi COPYONLY)

# copy database to build directory for testing
configure_file(SimpleAdvect.py SimpleAdvect.py COPYONLY)

get_target_property(PhreeqcRM_SOURCES PhreeqcRM SOURCES)

set(phreeqcrm_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../../PhreeqcRM.i)

set_source_files_properties(${phreeqcrm_INTERFACE} PROPERTIES SWIG_FLAGS "-Wall")
set_source_files_properties(${phreeqcrm_INTERFACE} PROPERTIES CPLUSPLUS ON)

swig_add_library(
  ${PYTHON_TARGET_NAME}
  LANGUAGE python
  TYPE MODULE
  SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../../PhreeqcRM.i ${PhreeqcRM_SOURCES}
)

target_include_directories(${PYTHON_TARGET_NAME}
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../..
)

target_link_libraries(${PYTHON_TARGET_NAME}
  PRIVATE
    Python::Python
    Python::NumPy
)

add_test(NAME SimpleAdvect
  COMMAND ${Python_EXECUTABLE} SimpleAdvect.py
)

set_tests_properties(SimpleAdvect PROPERTIES ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:${PYTHON_TARGET_NAME}>")
