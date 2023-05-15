set(PYTHON_TARGET_NAME phreeqcrmpy)

set(phreeqcrm_INTERFACE PhreeqcRM.i)

set(phreeqcrm_INTERFACE_FULLPATH ${CMAKE_CURRENT_SOURCE_DIR}/../../${phreeqcrm_INTERFACE})

# copy PhreeqcRM.i in order to not have pyfragments.swg be automatically included
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../${phreeqcrm_INTERFACE} ${phreeqcrm_INTERFACE} COPYONLY)

# @todo
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../pyfragments.swg.in pyfragments.swg COPYONLY)

# copy database to build directory for testing
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../../database/phreeqc.dat phreeqc.dat COPYONLY)

# copy advect.pqi to build directory for testing
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../../Tests/advect.pqi advect.pqi COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../../Tests/all_reactants.pqi all_reactants.pqi COPYONLY)

# copy SimpleAdvect.py to build directory for testing
configure_file(SimpleAdvect.py SimpleAdvect.py COPYONLY)

# copy yamlphreeqcrm.py to build directory for testing
configure_file(yamlphreeqcrm.py yamlphreeqcrm.py COPYONLY)

# copy AdvectBMI_py.py to build directory for testing
configure_file(AdvectBMI_py.py AdvectBMI_py.py COPYONLY)
configure_file(WriteYAMLFile_py.py WriteYAMLFile_py.py COPYONLY)
configure_file(TestAllMethods_py.py TestAllMethods_py.py COPYONLY)

get_target_property(PhreeqcRM_SOURCES PhreeqcRM SOURCES)

set(CMAKE_SWIG_FLAGS "-Wall" )
set_source_files_properties(${phreeqcrm_INTERFACE_FULLPATH} PROPERTIES CPLUSPLUS ON)

if(NOT PHREEQCRM_DISABLE_OPENMP)
  if(OPENMP_FOUND)
    list(APPEND CMAKE_SWIG_FLAGS "-DUSE_OPENMP")
  endif()
endif()

if(PHREEQCRM_WITH_YAML_CPP AND yaml-cpp_FOUND)
  list(APPEND CMAKE_SWIG_FLAGS "-DUSE_YAML")
endif()

# CMAKE_SWIG_FLAGS must be fully set before swig_add_library is called

swig_add_library(
  ${PYTHON_TARGET_NAME}
  LANGUAGE python
  TYPE MODULE
  SOURCES ${phreeqcrm_INTERFACE_FULLPATH} ${PhreeqcRM_SOURCES}
)

target_compile_features(${PYTHON_TARGET_NAME} PUBLIC cxx_std_11)

if(NOT PHREEQCRM_DISABLE_OPENMP)
  if(OPENMP_FOUND)
    target_compile_definitions(${PYTHON_TARGET_NAME} PRIVATE USE_OPENMP)
    target_link_libraries(${PYTHON_TARGET_NAME} PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

if(PHREEQCRM_WITH_YAML_CPP AND yaml-cpp_FOUND)
  target_compile_definitions(${PYTHON_TARGET_NAME} PUBLIC USE_YAML)
  target_link_libraries(${PYTHON_TARGET_NAME} PUBLIC yaml-cpp)
endif()

if (PHREEQCRM_NO_UTF8_ENCODING)
  target_compile_definitions(${PYTHON_TARGET_NAME} PRIVATE NO_UTF8_ENCODING)
endif()

target_link_libraries(${PYTHON_TARGET_NAME}
  PRIVATE
    Python::Python
    Python::NumPy
)

get_target_property(PhreeqcRM_INC PhreeqcRM INCLUDE_DIRECTORIES)
target_include_directories(${PYTHON_TARGET_NAME}
  PUBLIC
    ${PhreeqcRM_INC}
)

get_property(IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

if(NOT CMAKE_DEBUG_POSTFIX)
  # SWIG BUG -- PyInit__phreeqcrmpy is misnamed when either of these properties is used:
  # set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES DEBUG_POSTFIX "_d")
  # set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES OUTPUT_NAME_DEBUG "${PYTHON_TARGET_NAME}_d")
  if(WIN32)
    # copy _phreeqcrmpy.pyd to _phreeqcrmpy_d.pyd
    add_custom_command(TARGET ${PYTHON_TARGET_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:${PYTHON_TARGET_NAME}>
        $<TARGET_FILE_DIR:${PYTHON_TARGET_NAME}>/_$<TARGET_FILE_BASE_NAME:${PYTHON_TARGET_NAME}>$<$<CONFIG:Debug>:_d>$<TARGET_FILE_SUFFIX:${PYTHON_TARGET_NAME}>
    )
  endif()
endif()

if(IS_MULTI_CONFIG)
  add_test(NAME SimpleAdvect.py
    COMMAND ${Python_EXECUTABLE} SimpleAdvect.py
    CONFIGURATIONS Release RelWithDebInfo MinSizeRel
  )
  set_tests_properties(SimpleAdvect.py PROPERTIES ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:${PYTHON_TARGET_NAME}>")
else()
  if("${CMAKE_BUILD_TYPE}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    add_test(NAME SimpleAdvect.py
      COMMAND ${Python_EXECUTABLE} SimpleAdvect.py
    )
  endif()
endif()

if(WIN32)
  string(REPLACE python.exe python_d.exe Python_d_EXECUTABLE ${Python_EXECUTABLE})

  if(IS_MULTI_CONFIG)
    add_test(NAME SimpleAdvect_d.py
      COMMAND ${Python_d_EXECUTABLE} SimpleAdvect.py
      CONFIGURATIONS Debug
    )
    set_tests_properties(SimpleAdvect_d.py PROPERTIES ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:${PYTHON_TARGET_NAME}>")
  else()
    if("${CMAKE_BUILD_TYPE}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
      add_test(NAME SimpleAdvect_d.py
        COMMAND ${Python_d_EXECUTABLE} SimpleAdvect.py
      )
    endif()
  endif()
  if(NOT EXISTS ${Python_d_EXECUTABLE})
    message(STATUS "Debug python interpreter not found")
    if(IS_MULTI_CONFIG)
      set_tests_properties(SimpleAdvect_d.py PROPERTIES DISABLED TRUE)
    else()
      if("${CMAKE_BUILD_TYPE}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
        set_tests_properties(SimpleAdvect_d.py PROPERTIES DISABLED TRUE)
      endif()
    endif()
  endif()
endif()

if(MSVC)
  # setup python Debugging
  # Configuration Properties->Debugging

  # Get the base name
  get_filename_component(base_name ${Python_EXECUTABLE} NAME_WLE)

  # Get the extension
  get_filename_component(ext ${Python_EXECUTABLE} LAST_EXT)

  # Get the directory of the input file
  get_filename_component(dir_name ${Python_EXECUTABLE} DIRECTORY)

  set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES VS_DEBUGGER_COMMAND "${dir_name}/${base_name}$<$<CONFIG:Debug>:_d>${ext}")
  # set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES VS_DEBUGGER_ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:${PYTHON_TARGET_NAME}>")
  set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES VS_DEBUGGER_ENVIRONMENT "PYTHONPATH=$(OutDirFullPath)")
  set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES VS_DEBUGGER_COMMAND_ARGUMENTS "SimpleAdvect.py")
endif()

# include(CMakePrintHelpers)
# cmake_print_variables(Python_EXECUTABLE)
# cmake_print_variables(Python_d_EXECUTABLE)
# cmake_print_variables(IS_MULTI_CONFIG)
