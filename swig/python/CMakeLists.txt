set(PYTHON_TARGET_NAME swig-python)

set(phreeqcrm_INTERFACE PhreeqcRM.i)

set(phreeqcrm_INTERFACE_FULLPATH ${CMAKE_CURRENT_SOURCE_DIR}/../${phreeqcrm_INTERFACE})

# This fixes the numpy integer swig incompatibility bug
# see https://numpy.org/doc/stable/reference/swig.interface-file.html#numpy-array-scalars-and-swig
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../pyfragments.swg.in pyfragments.swg COPYONLY)

# copy database to build directory for testing
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../database/phreeqc.dat phreeqc.dat COPYONLY)
if(DEFINED SKBUILD)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../database/phreeqc.dat ${CMAKE_CURRENT_SOURCE_DIR}/phreeqc.dat COPYONLY)
endif()

# copy advect.pqi to build directory for testing
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../Tests/advect.pqi advect.pqi COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../Tests/all_reactants.pqi all_reactants.pqi COPYONLY)
if(DEFINED SKBUILD)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../../Tests/advect.pqi ${CMAKE_CURRENT_SOURCE_DIR}/advect.pqi COPYONLY)
endif()

# copy SimpleAdvect.py to build directory for testing
configure_file(SimpleAdvect.py SimpleAdvect.py COPYONLY)

# copy yamlphreeqcrm.py to build directory for testing
configure_file(yamlphreeqcrm.py yamlphreeqcrm.py COPYONLY)

# copy AdvectBMI_py.py to build directory for testing
configure_file(AdvectBMI_py.py AdvectBMI_py.py COPYONLY)
configure_file(WriteYAMLFile_py.py WriteYAMLFile_py.py COPYONLY)
configure_file(TestAllMethods_py.py TestAllMethods_py.py COPYONLY)

# copy pytest files to build directory for testing
configure_file(conftest.py conftest.py COPYONLY)
configure_file(constants.py constants.py COPYONLY)
configure_file(test_basic.py test_basic.py COPYONLY)
configure_file(test_get_value.py test_get_value.py COPYONLY)
configure_file(test_init.py test_init.py COPYONLY)
configure_file(test_irf.py test_irf.py COPYONLY)
configure_file(test_set_value.py test_set_value.py COPYONLY)

get_target_property(PhreeqcRM_SOURCES PhreeqcRM SOURCES)

set_property(SOURCE ${phreeqcrm_INTERFACE_FULLPATH} APPEND PROPERTY COMPILE_OPTIONS "-w401")
set_source_files_properties(${phreeqcrm_INTERFACE_FULLPATH} PROPERTIES CPLUSPLUS ON)

# if any of these files change force swig to be re-run
set_property(SOURCE ${phreeqcrm_INTERFACE_FULLPATH} APPEND PROPERTY DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/BMIPhreeqcRM_docstrings.swg
  ${CMAKE_CURRENT_SOURCE_DIR}/PhreeqcRM_docstrings.swg
  ${PhreeqcRM_SOURCE_DIR}/src/BMIVariant.h
  ${PhreeqcRM_SOURCE_DIR}/src/IrmResult.h
  ${PhreeqcRM_SOURCE_DIR}/src/PhreeqcRM.h
  ${PhreeqcRM_SOURCE_DIR}/src/BMIPhreeqcRM.h
)

if(NOT PHREEQCRM_DISABLE_OPENMP)
  if(OPENMP_FOUND)
    set_property(SOURCE ${phreeqcrm_INTERFACE_FULLPATH} APPEND PROPERTY COMPILE_DEFINITIONS USE_OPENMP)
  endif()
endif()

if((PHREEQCRM_WITH_YAML_CPP AND yaml-cpp_FOUND) OR (DEFINED SKBUILD))
  set_property(SOURCE ${phreeqcrm_INTERFACE_FULLPATH} APPEND PROPERTY COMPILE_DEFINITIONS USE_YAML)
  set_property(SOURCE ${phreeqcrm_INTERFACE_FULLPATH} APPEND PROPERTY DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/YAMLPhreeqcRM.swg
  )
endif()

# CMAKE_SWIG_FLAGS must be fully set before swig_add_library is called

swig_add_library(
  ${PYTHON_TARGET_NAME}
  TYPE MODULE
  LANGUAGE python
  SOURCES ${phreeqcrm_INTERFACE_FULLPATH} ${PhreeqcRM_SOURCES}
)

if(WIN32)
  set_property(TARGET ${PYTHON_TARGET_NAME} PROPERTY SUFFIX ".${Python_SOABI}.pyd")
else()
  set_property(TARGET ${PYTHON_TARGET_NAME}
               PROPERTY SUFFIX ".${Python_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}")
endif()

set_property(TARGET ${PYTHON_TARGET_NAME} PROPERTY OUTPUT_NAME "phreeqcrm")

target_compile_features(${PYTHON_TARGET_NAME} PUBLIC cxx_std_11)

target_compile_definitions(${PYTHON_TARGET_NAME} PRIVATE IRM_DLL_STATIC_DEFINE)

# iphreeqc defs
target_compile_definitions(${PYTHON_TARGET_NAME} PRIVATE SWIG_SHARED_OBJ)
target_compile_definitions(${PYTHON_TARGET_NAME} PRIVATE USE_PHRQ_ALLOC)
target_compile_definitions(${PYTHON_TARGET_NAME} PRIVATE IPhreeqc_EXPORTS)

# numpy defs
target_compile_definitions(${PYTHON_TARGET_NAME} PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

# Turn off Incremental Linking for Debug and RelWithDebInfo
if(MSVC)
  target_link_options(${PYTHON_TARGET_NAME} PRIVATE $<$<CONFIG:Debug,RelWithDebInfo>:/INCREMENTAL:NO>)
endif()

if(NOT PHREEQCRM_DISABLE_OPENMP)
  if(OPENMP_FOUND)
    target_compile_definitions(${PYTHON_TARGET_NAME} PRIVATE USE_OPENMP)
    target_link_libraries(${PYTHON_TARGET_NAME} PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

if(PHREEQCRM_WITH_YAML_CPP AND yaml-cpp_FOUND)
  target_compile_definitions(${PYTHON_TARGET_NAME} PUBLIC USE_YAML)
  if(TARGET yaml-cpp::yaml-cpp)
    # The CD scikit_build_core build uses the latest (as of 2023-07-19) commit
    # b8882652fcbeba4c00dec019a39da91e702e474e of yaml-cpp which uses the
    # yaml-cpp:: namespace
    target_link_libraries(${PYTHON_TARGET_NAME} PUBLIC yaml-cpp::yaml-cpp)
  else()
    # yaml-cpp <= 0.7.0 doesn't use a namespace
    # as of 2023-07-19:
    #   macos        brew    0.7.0
    #   windows      vcpkg   0.7.0#1
    #   ubuntu-20.04 apt-get 0.6
    #   ubuntu-22.04 apt-get 0.7
    target_link_libraries(${PYTHON_TARGET_NAME} PUBLIC yaml-cpp)
  endif()
endif()

if(DEFINED SKBUILD)
  target_compile_definitions(${PYTHON_TARGET_NAME} PUBLIC USE_YAML)
  if(TARGET yaml-cpp::yaml-cpp)
    message(STATUS "SKBUILD found yaml-cpp::yaml-cpp")
    target_link_libraries(${PYTHON_TARGET_NAME} PUBLIC yaml-cpp::yaml-cpp)
  else()
    message(STATUS "SKBUILD didn't find yaml-cpp::yaml-cpp")
    target_link_libraries(${PYTHON_TARGET_NAME} PUBLIC yaml-cpp)
  endif()
endif()

if(DEFINED SKBUILD)
  target_compile_definitions(${PYTHON_TARGET_NAME} PUBLIC USE_GZ)
  target_link_libraries(${PYTHON_TARGET_NAME} PUBLIC zlibstatic)
endif()

if(PHREEQCRM_NO_UTF8_ENCODING)
  target_compile_definitions(${PYTHON_TARGET_NAME} PRIVATE NO_UTF8_ENCODING)
endif()

if(WIN32 AND Python_LIBRARY_DEBUG)
  # Since we're using SWIG_PYTHON_INTERPRETER_NO_DEBUG in PhreeqcRM.i we only
  # link against the Release version of the Python library
  #
  # Attempting to run python tests with phreeqcrmpy Debug build on windows yields:
  #   Fatal Python error: _PyInterpreterState_GET: the function must be called with the GIL held, but the GIL is released (the current Python thread state is NULL)
  message(STATUS "Force ${PYTHON_TARGET_NAME} to link against the Release version of the Python library on Windows")

  # Python include directory
  target_include_directories(${PYTHON_TARGET_NAME}
    PRIVATE
      ${Python_INCLUDE_DIRS}
  )

  # Numpy include directory
  get_target_property(NumPy_INC_DIRS Python::NumPy INTERFACE_INCLUDE_DIRECTORIES)
  target_include_directories(${PYTHON_TARGET_NAME}
    PRIVATE
      ${NumPy_INC_DIRS}
  )

  # Python link library
  target_link_libraries(${PYTHON_TARGET_NAME}
    PRIVATE
      ${Python_LIBRARY_RELEASE}
  )
else()
  target_link_libraries(${PYTHON_TARGET_NAME}
    PRIVATE
      Python::Module
      Python::NumPy
  )
endif()


get_target_property(PhreeqcRM_INC PhreeqcRM INCLUDE_DIRECTORIES)
target_include_directories(${PYTHON_TARGET_NAME}
  PUBLIC
    ${PhreeqcRM_INC}
)

get_property(IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

# if(NOT CMAKE_DEBUG_POSTFIX)
#   # SWIG BUG -- PyInit__phreeqcrmpy is misnamed when either of these properties is used:
#   # set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES DEBUG_POSTFIX "_d")
#   # set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES OUTPUT_NAME_DEBUG "${PYTHON_TARGET_NAME}_d")
#   if(WIN32)
#     # copy _phreeqcrmpy.pyd to _phreeqcrmpy_d.pyd
#     add_custom_command(TARGET ${PYTHON_TARGET_NAME} POST_BUILD
#       COMMAND ${CMAKE_COMMAND} -E copy_if_different
#         $<TARGET_FILE:${PYTHON_TARGET_NAME}>
#         $<TARGET_FILE_DIR:${PYTHON_TARGET_NAME}>/_$<TARGET_FILE_BASE_NAME:${PYTHON_TARGET_NAME}>$<$<CONFIG:Debug>:_d>$<TARGET_FILE_SUFFIX:${PYTHON_TARGET_NAME}>
#     )
#   endif()
# endif()

# SimpleAdvect.py test
add_test(NAME SimpleAdvect.py
  COMMAND ${Python_EXECUTABLE} SimpleAdvect.py
)
set_tests_properties(SimpleAdvect.py
  PROPERTIES
    ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:${PYTHON_TARGET_NAME}>"
)

# WriteYAMLFile_py.py test
add_test(NAME WriteYAMLFile_py.py
  COMMAND ${Python_EXECUTABLE} WriteYAMLFile_py.py
)
set_tests_properties(WriteYAMLFile_py.py PROPERTIES ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:${PYTHON_TARGET_NAME}>")

# AdvectBMI_py.py test
add_test(NAME AdvectBMI_py.py
  COMMAND ${Python_EXECUTABLE} AdvectBMI_py.py
)
set_tests_properties(AdvectBMI_py.py
  PROPERTIES 
    ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:${PYTHON_TARGET_NAME}>"
    DEPENDS "WriteYAMLFile_py.py"
    # AdvectBMI_py.py requires PHREEQCRM_WITH_YAML_CPP=ON
    # DISABLED "$<$<NOT:$<BOOL:${PHREEQCRM_WITH_YAML_CPP}>>:TRUE>"
)

# TestAllMethods_py.py test
add_test(NAME TestAllMethods_py.py
  COMMAND ${Python_EXECUTABLE} TestAllMethods_py.py
)
set_tests_properties(TestAllMethods_py.py
  PROPERTIES 
    ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:${PYTHON_TARGET_NAME}>"
)


# Pytest
add_test(NAME run_pytest
  COMMAND ${Python_EXECUTABLE} -m pytest
)
set_tests_properties(run_pytest
  PROPERTIES ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:${PYTHON_TARGET_NAME}>"
    # Some Pytest tests require PHREEQCRM_WITH_YAML_CPP=ON
    # DISABLED "$<$<NOT:$<BOOL:${PHREEQCRM_WITH_YAML_CPP}>>:TRUE>"
)


if(MSVC)
  # setup python Debugging
  # Configuration Properties->Debugging

  # Get the base name
  get_filename_component(base_name ${Python_EXECUTABLE} NAME_WLE)

  # Get the extension
  get_filename_component(ext ${Python_EXECUTABLE} LAST_EXT)

  # Get the directory of the input file
  get_filename_component(dir_name ${Python_EXECUTABLE} DIRECTORY)

  # set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES VS_DEBUGGER_COMMAND "${dir_name}/${base_name}$<$<CONFIG:Debug>:_d>${ext}")
  set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES VS_DEBUGGER_COMMAND ${Python_EXECUTABLE})
  set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES VS_DEBUGGER_ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:${PYTHON_TARGET_NAME}>")
  set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES VS_DEBUGGER_COMMAND_ARGUMENTS "SimpleAdvect.py")
endif()

if(DEFINED SKBUILD)
  # install _phreeqcrm.${Python_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}
  # ie _phreeqcrm.cp38-win_amd64.pyd
  install(TARGETS ${PYTHON_TARGET_NAME} DESTINATION ${SKBUILD_PROJECT_NAME})

  # install phreeqcrm.py
  get_target_property(support_files ${PYTHON_TARGET_NAME} SWIG_SUPPORT_FILES)
  install(FILES ${support_files} DESTINATION ${SKBUILD_PROJECT_NAME})

  # install yamlphreeqcrm.py
  install(FILES yamlphreeqcrm.py DESTINATION ${SKBUILD_PROJECT_NAME})
endif()

# # Get all propreties that cmake supports
# if(NOT CMAKE_PROPERTY_LIST)
#     execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)
    
#     # Convert command output into a CMake list
#     string(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
#     string(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
#     list(REMOVE_DUPLICATES CMAKE_PROPERTY_LIST)
# endif()
    
# function(print_properties)
#     message("CMAKE_PROPERTY_LIST = ${CMAKE_PROPERTY_LIST}")
# endfunction()
    
# function(print_target_properties target)
#     if(NOT TARGET ${target})
#       message(STATUS "There is no target named '${target}'")
#       return()
#     endif()

#     foreach(property ${CMAKE_PROPERTY_LIST})
#         string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" property ${property})

#         # Fix https://stackoverflow.com/questions/32197663/how-can-i-remove-the-the-location-property-may-not-be-read-from-target-error-i
#         if(property STREQUAL "LOCATION" OR property MATCHES "^LOCATION_" OR property MATCHES "_LOCATION$")
#             continue()
#         endif()

#         get_property(was_set TARGET ${target} PROPERTY ${property} SET)
#         if(was_set)
#             get_target_property(value ${target} ${property})
#             message("${target} ${property} = ${value}")
#         endif()
#     endforeach()
# endfunction()


# include(CMakePrintHelpers)
# cmake_print_variables(is_yaml_not_set)
# cmake_print_variables(Python_EXECUTABLE)
# cmake_print_variables(Python_d_EXECUTABLE)
# cmake_print_variables(IS_MULTI_CONFIG)
# cmake_print_variables(Python_Development_FOUND)

# cmake_print_variables(NumPy_INCLUDE_DIRS)

# cmake_print_variables(Python_INCLUDE_DIRS)
# cmake_print_variables(Python_LIBRARIES)
# cmake_print_variables(Python_LIBRARY_RELEASE)
# cmake_print_variables(Python_LIBRARY_DEBUG)

# print_target_properties(Python::Python)
# print_target_properties(Python::NumPy)
# print_target_properties(Python::Module)

# cmake_print_variables(Python::Python)
