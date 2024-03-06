USE_YAML = true

ifeq ($(OS),Windows_NT)
	DOXYGEN ?= "C:/Program Files/doxygen/bin/doxygen.exe"
else
	DOXYGEN ?= doxygen
endif

all:
	cp ../src/BMI_interface.F90 ../src/BMI_interface.F90.save
	cp ../src/RM_interface.F90 ../src/RM_interface.F90.save
	cp ../src/YAML_interface.F90 ../src/YAML_interface.F90.save
	./lower_case.sh
	$(DOXYGEN)
	cp ../src/BMI_interface.F90.save ../src/BMI_interface.F90
	cp ../src/RM_interface.F90.save ../src/RM_interface.F90
	cp ../src/YAML_interface.F90.save ../src/YAML_interface.F90
	#./upper_case.sh
	./upper_case2.sh
	#./upper_case2a.sh

win:
	lower_case.bat
	$(DOXYGEN)
	upper_case.bat
	#svn revert ..\src\RM_interface.F90
