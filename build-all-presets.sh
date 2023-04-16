#!/bin/bash

declare -a allowed=(
    ninja
    ninja-vcpkg
    vs2019
    vs2019-vcpkg
)

if [ $# -ne 1 ]; then
  echo "Usage: $0 <arg>"
  echo "where <arg> is one of: (${allowed[*]})"
  exit 1
fi

if [[ ! " ${allowed[*]} " =~ " ${1} " ]]; then
  echo "Error: First argument must be one of: (${allowed[*]})"
  exit 1
fi

declare -a presets=(
  static-debug
  static-release
  static-yaml-debug
  static-yaml-release
  static-yaml-python-debug
  static-yaml-python-release
  shared-debug
  shared-release
  shared-yaml-debug
  shared-yaml-release
  shared-yaml-python-debug
  shared-yaml-python-release
)

for preset in "${presets[@]}"
do 
  cmake --preset="$1-$preset"
  cmake --build --preset="$1-$preset"
done

for preset in "${presets[@]}"
do 
  ctest --preset="$1-$preset"
done
