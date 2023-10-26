cmake_minimum_required(VERSION 3.9)

project(openmp Fortran)

find_package(OpenMP REQUIRED)
  
add_executable(openmp openmp.f90)
 
target_link_libraries(openmp PUBLIC OpenMP::OpenMP_Fortran)

enable_testing()

add_test(NAME openmp COMMAND openmp)
