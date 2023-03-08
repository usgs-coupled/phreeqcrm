EXTRA_DIST = CMakeLists.txt

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/IPhreeqcPhast -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp  -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords -I$(top_builddir)

TESTS = test
check_PROGRAMS = test

test_SOURCES =\
	Advect_c.c\
	AdvectBMI_cpp.cpp\
	Advect_cpp.cpp\
	Gas_c.c\
	Gas_cpp.cpp\
	main.cpp\
	SimpleAdvect_c.c\
	SimpleAdvect_cpp.cpp\
	Species_c.c\
	Species_cpp.cpp\
	WriteYAMLFile_cpp.cpp


# Ubuntu 22.04.1
# sudo apt install libyaml-cpp-dev
test_LDADD = $(top_builddir)/libphreeqcrm.la -lyaml-cpp
test_LINK = $(CXXLINK)

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
	AdvectBMI_cpp_units_utility.txt\
	AdvectBMI_cpp_units_worker.txt\
	AdvectBMI_cpp_utility.txt\
	AdvectBMI_cpp.chem.txt\
	AdvectBMI_cpp.dmp\
	AdvectBMI_cpp.log.txt\
	AdvectBMI_cpp.yaml\
	Gas_c_utility.txt\
	Gas_c.chem.txt\
	Gas_c.dmp\
	Gas_c.log.txt\
	Gas_cpp_utility.txt\
	Gas_cpp.chem.txt\
	Gas_cpp.dmp\
	Gas_cpp.log.txt\
	SimpleAdvect_c.chem.txt\
	SimpleAdvect_c.log.txt\
	SimpleAdvect_cpp.chem.txt\
	SimpleAdvect_cpp.log.txt\
	Species_c_utility.txt\
	Species_c.chem.txt\
	Species_c.dmp\
	Species_c.log.txt\
	Species_cpp_utility.txt\
	Species_cpp.chem.txt\
	Species_cpp.dmp\
	Species_cpp.log.txt\
	Units_Worker.chem.txt\
	Units_Worker.log.txt\
	Utility_c.txt\
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
