cmake_minimum_required(VERSION 3.9)

project(hw Fortran)
  
add_executable(hw hw.f90)

enable_testing()

add_test(NAME hw COMMAND hw)
