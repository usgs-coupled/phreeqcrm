EXTRA_DIST = CMakeLists.txt

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/IPhreeqcPhast -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp  -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords -I$(top_builddir)

TESTS = test
check_PROGRAMS = test

test_SOURCES =\
	Advect_cpp.cpp\
	Advect_c.c\
	Advect_cpp.cpp\
	Gas_c.c\
	Gas_cpp.cpp\
	main.cpp\
	SimpleAdvect_cpp.cpp\
	Species_c.c\
	Species_cpp.cpp\
	WriteYAMLFile_cpp.cpp\
	YAMLPhreeqcRM.cpp\
	YAMLPhreeqcRM.h


# Ubuntu 22.04.1
# sudo apt install libyaml-cpp-dev
test_LDADD = $(top_builddir)/libphreeqcrm.la -lyaml-cpp
test_LINK = $(CXXLINK)

CLEANFILES =\
	AdvectBMI_cpp.chem.out\
	AdvectBMI_cpp.dmp\
	AdvectBMI_cpp.log.out\
	AdvectBMI_cpp_units_utility.out\
	AdvectBMI_cpp_units_worker.out\
	AdvectBMI_cpp.yaml\
	Advect_c.chem.out\
	Advect_c.dmp\
	Advect_c.log.out\
	Advect_cpp.chem.out\
	Advect_cpp.log.out\
	Advect_cpp.dmp\
	Advect_cpp_units_utility.out\
	Advect_cpp_units_worker.out\
	Gas_c.chem.out\
	Gas_c.dmp\
	Gas_c.log.out\
	Gas_c_utility.out\
	Gas_cpp.chem.out\
	Gas_cpp.dmp\
	Gas_cpp.log.out\
	Gas_cpp_utility.out\
	SimpleAdvect_cpp.chem.out\
	SimpleAdvect_cpp.log.out\
	Species_c.chem.out\
	Species_c.dmp\
	Species_c.log.out\
	Species_c_utility.out\
	Species_cpp.chem.out\
	Species_cpp.dmp\
	Species_cpp.log.out\
	Species_cpp_utility.out\
	Units_Worker.chem.out\
	Units_Worker.log.out\
	Utility_c.out\
	Utility_cpp.out

if PHREEQCRM_TEST_FORTRAN

test_SOURCES+=\
	Advect_f90.F90\
	Gas_f90.F90\
	SimpleAdvect_f90.F90\
	Species_f90.F90\
	../src/RM_interface.F90\
	../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.F90

test_LDADD+=$(FCLIBS)

mostlyclean-local: clean-modules
.PHONY: clean-modules
clean-modules:
	-test -z $(FC_MODEXT) || $(RM) *.$(FC_MODEXT)

advection_f90.$(OBJEXT) : ../src/RM_interface.$(OBJEXT) ../src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.$(OBJEXT)

species_f90.$(OBJEXT) : ../src/RM_interface.$(OBJEXT)

gas_f90.$(OBJEXT) : ../src/RM_interface.$(OBJEXT)


endif

# We must override the default PPFCCOMPILE since ifort can't
# handle the FC_FUNC define correctly so $(DEFS) is replaced
# by $(FCDEFS)
#
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
        $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)
LTPPFCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=compile $(FC) $(FCDEFS) \
        $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
        $(AM_FCFLAGS) $(FCFLAGS)
