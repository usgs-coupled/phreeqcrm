EXTRA_DIST = CMakeLists.txt

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/IPhreeqcPhast -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp  -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords -I$(top_builddir)

TESTS = TestRM TestRMdtor TestBMIdtor
check_PROGRAMS = TestRM TestRMdtor TestBMIdtor

TestRM_SOURCES =\
	Advect_c.c\
	AdvectBMI_cpp.cpp\
	AdvectBMI_cpp_test.cpp\
	Advect_cpp.cpp\
	Gas_c.c\
	Gas_cpp.cpp\
	main.cpp\
	SimpleAdvect_c.c\
	SimpleAdvect_cpp.cpp\
	Species_c.c\
	Species_cpp.cpp\
	TestAllMethods_c.c\
	TestAllMethods_cpp.cpp\
	WriteYAMLFile_cpp.cpp\
	WriteYAMLFile_cpp_test.cpp

TestRM_LDADD = $(top_builddir)/libphreeqcrm.la
TestRM_LINK = $(CXXLINK)

# TestBMIdtor

TestBMIdtor_SOURCES = TestBMIdtor.cpp

TestBMIdtor_LDADD = $(top_builddir)/libphreeqcrm.la
TestBMIdtor_LINK = $(CXXLINK)

# TestRMdtor

TestRMdtor_SOURCES = TestRMdtor.cpp

TestRMdtor_LDADD = $(top_builddir)/libphreeqcrm.la
TestRMdtor_LINK = $(CXXLINK)


CLEANFILES =\
	Advect_c_utility.txt\
	Advect_c.chem.txt\
	Advect_c.dmp\
	Advect_c.log.txt\
	Advect_cpp_units_utility.txt\
	Advect_cpp_units_worker.chem.txt\
	Advect_cpp_units_worker.log.txt\
	Advect_cpp_units_worker.txt\
	Advect_cpp_utility.txt\
	Advect_cpp.chem.txt\
	Advect_cpp.dmp\
	Advect_cpp.log.txt\
	Advect_f90_utility.txt\
	Advect_f90.chem.txt\
	Advect_f90.dmp\
	Advect_f90.log.txt\
	AdvectBMI_cpp_test.chem.txt\
	AdvectBMI_cpp_test.log.txt\
	AdvectBMI_cpp_test.yaml\
	AdvectBMI_cpp_units_utility.txt\
	AdvectBMI_cpp_units_worker.txt\
	AdvectBMI_cpp_utility.txt\
	AdvectBMI_cpp.chem.txt\
	AdvectBMI_cpp.dmp\
	AdvectBMI_cpp.log.txt\
	AdvectBMI_cpp.yaml\
	AdvectBMI_f90_test.chem.txt\
	AdvectBMI_f90_test.log.txt\
	AdvectBMI_f90_test.yaml\
	AdvectBMI_f90.chem.txt\
	AdvectBMI_f90.dmp\
	AdvectBMI_f90.log.txt\
	AdvectBMI_f90.yaml\
	AdvectBMI_test_py.chem.txt\
	AdvectBMI_test_py.dump\
	AdvectBMI_test_py.log.txt\
	Advectcpp_utility.txt\
	Advectcpp.dmp\
	Gas_c_utility.txt\
	Gas_c.chem.txt\
	Gas_c.dmp\
	Gas_c.log.txt\
	Gas_cpp_utility.txt\
	Gas_cpp.chem.txt\
	Gas_cpp.dmp\
	Gas_cpp.log.txt\
	Gas_f90.chem.txt\
	Gas_f90.log.txt\
	SimpleAdvect_c.chem.txt\
	SimpleAdvect_c.log.txt\
	SimpleAdvect_cpp.chem.txt\
	SimpleAdvect_cpp.log.txt\
	SimpleAdvect_f90.chem.txt\
	SimpleAdvect_f90.log.txt\
	Species_c_utility.txt\
	Species_c.chem.txt\
	Species_c.dmp\
	Species_c.log.txt\
	Species_cpp_utility.txt\
	Species_cpp.chem.txt\
	Species_cpp.dmp\
	Species_cpp.log.txt\
	Species_f90_utility.txt\
	Species_f90.chem.txt\
	Species_f90.dmp\
	Species_f90.log.txt\
	TestAllMethods_c.chem.txt\
	TestAllMethods_c.log.txt\
	TestAllMethods_c.yaml\
	TestAllMethods_cpp.chem.txt\
	TestAllMethods_cpp.dump\
	TestAllMethods_cpp.log.txt\
	TestAllMethods_cpp.yaml\
	TestAllMethods_f90.chem.txt\
	TestAllMethods_f90.dmp\
	TestAllMethods_f90.log.txt\
	TestAllMethods_f90.yaml\
	TestAllMethods_py.chem.txt\
	TestAllMethods_py.dump\
	TestAllMethods_py.log.txt\
	testbmi_cpp.yaml\
	Units_Worker.chem.txt\
	Units_Worker.log.txt\
	Utility_c.txt\
	Utility_cpp.out

FCDEFS=

if PHREEQCRM_TEST_FORTRAN

DEFS+=-DTEST_FORTRAN
FCDEFS+=-DTEST_FORTRAN

TestRM_SOURCES+=\
	../src/BMI_interface.F90\
	../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.F90\
	../src/RM_interface.F90\
	../src/YAML_interface.F90\
	Advect_f90.F90\
	AdvectBMI_f90.F90\
	AdvectBMI_f90_test.F90\
	Gas_f90.F90\
	SimpleAdvect_f90.F90\
	Species_f90.F90\
	TestAllMethods_f90.F90\
	WriteYAMLFile_f90.F90\
	WriteYAMLFile_f90_test.F90

TestRM_LDADD+=$(FCLIBS)

mostlyclean-local: clean-modules
.PHONY: clean-modules
clean-modules:
	-test -z $(FC_MODEXT) || $(RM) *.$(FC_MODEXT)

AdvectBMI_f90.$(OBJEXT) : AdvectBMI_f90.F90 Advect_f90.$(OBJEXT) ../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.$(OBJEXT) ../src/BMI_interface.$(OBJEXT)

# mydata.mod
Advect_f90.$(OBJEXT) : Advect_f90.F90 ../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.$(OBJEXT) ../src/RM_interface.$(OBJEXT)

Gas_f90.$(OBJEXT) : Gas_f90.F90 Advect_f90.$(OBJEXT) ../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.$(OBJEXT) ../src/RM_interface.$(OBJEXT)

SimpleAdvect_f90.$(OBJEXT) : SimpleAdvect_f90.F90 Advect_f90.$(OBJEXT) ../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.$(OBJEXT) ../src/RM_interface.$(OBJEXT)

Species_f90.$(OBJEXT) : Species_f90.F90 ../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.$(OBJEXT) ../src/RM_interface.$(OBJEXT)

# BMIPhreeqcRM.mod
../src/BMI_interface.$(OBJEXT) : ../src/BMI_interface.F90 ../src/RM_interface.$(OBJEXT)

AdvectBMI_f90_test.$(OBJEXT) : AdvectBMI_f90_test.F90 ../src/BMI_interface.$(OBJEXT)

# IPhreeqc.mod
../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.$(OBJEXT) : ../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.F90

# PhreeqcRM.mod
../src/RM_interface.$(OBJEXT) : ../src/RM_interface.F90

# YAML_interface.mod
../src/YAML_interface.$(OBJEXT) : ../src/YAML_interface.F90

TestAllMethods_f90.$(OBJEXT) : TestAllMethods_f90.F90 ../src/YAML_interface.$(OBJEXT) ../src/BMI_interface.$(OBJEXT)

WriteYAMLFile_f90.$(OBJEXT) : WriteYAMLFile_f90.F90 ../src/YAML_interface.$(OBJEXT)

WriteYAMLFile_f90_test.$(OBJEXT) : WriteYAMLFile_f90_test.F90 ../src/YAML_interface.$(OBJEXT)

endif

# We must override the default PPFCCOMPILE since ifort can't
# handle the FC_FUNC define correctly so $(DEFS) is replaced
# by $(FCDEFS)
#
if PHREEQCRM_WITH_YAMLCPP
# Ubuntu 22.04.1
# sudo apt install libyaml-cpp-dev
FCDEFS+=-DUSE_YAML
endif
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
        $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)
LTPPFCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=compile $(FC) $(FCDEFS) \
        $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
        $(AM_FCFLAGS) $(FCFLAGS)
