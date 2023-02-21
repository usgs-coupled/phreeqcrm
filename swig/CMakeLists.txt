cmake_minimum_required(VERSION 3.14) # SWIG_MODULE_NAME

# SWIG: use SWIG_MODULE_NAME property.
# for more run 'cmake --help-policy CMP0086'
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 NEW)
endif()

project(phreeqcrm-bmi)

find_package(SWIG REQUIRED)

find_package(Python REQUIRED COMPONENTS Interpreter Development)

# this defines swig_add_library and was set by FindSWIG aka find_package(SWIG ...)
# see https://cmake.org/cmake/help/latest/module/FindSWIG.html
# and https://cmake.org/cmake/help/latest/module/UseSWIG.html
include(${SWIG_USE_FILE})

# phreeqc.dat
configure_file(phreeqc.dat phreeqc.dat COPYONLY)

# units.pqi
configure_file(units.pqi units.pqi COPYONLY)

# advection_cpp.py
configure_file(advection_cpp.py advection_cpp.py COPYONLY)

# phreeqcrm-bmi.i
set_property(SOURCE phreeqcrm-bmi.i PROPERTY SWIG_MODULE_NAME phreeqcrm)
set_property(SOURCE phreeqcrm-bmi.i PROPERTY CPLUSPLUS ON)
swig_add_library(phreeqcrm LANGUAGE python SOURCES phreeqcrm-bmi.i)

target_link_libraries(phreeqcrm
  PRIVATE
    Python::Python
    PhreeqcRM::PhreeqcRM
)
