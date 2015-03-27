all:
	./lower_case.sh
	doxygen
	./upper_case.sh
	svn revert ../src/RM_interface.F90

win:
	lower_case.bat
	doxygen
	upper_case.bat
	svn revert ..\src\RM_interface.F90
