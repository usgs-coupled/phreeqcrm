DOXYGEN ?= doxygen

all:
	cp ../src/RM_interface.F90 ../src/RM_interface.F90.save
	./lower_case.sh
	$(DOXYGEN)
	./upper_case.sh
	mv ../src/RM_interface.F90.save ../src/RM_interface.F90

win:
	lower_case.bat
	$(DOXYGEN)
	upper_case.bat
	svn revert ..\src\RM_interface.F90
