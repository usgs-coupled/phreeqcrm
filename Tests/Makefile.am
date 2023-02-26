EXTRA_DIST = CMakeLists.txt

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/IPhreeqcPhast -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp  -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords -I$(top_builddir)

TESTS = test
check_PROGRAMS = test

test_SOURCES =\
	advection_bmi_cpp.cpp\
	advection_c.c\
	advection_cpp.cpp\
	gas_c.c\
	gas_cpp.cpp\
	main.cpp\
	SimpleAdvection.cpp\
	species_c.c\
	species_cpp.cpp\
	WriteYAMLFile.cpp\
	YAMLPhreeqcRM.cpp\
	YAMLPhreeqcRM.h


# Ubuntu 22.04.1
# sudo apt install libyaml-cpp-dev
test_LDADD = $(top_builddir)/libphreeqcrm.la -lyaml-cpp
test_LINK = $(CXXLINK)

CLEANFILES =\
	Advect_bmi_cpp.chem.txt\
	Advect_bmi_cpp.log.txt\
	Advect_c.chem.txt\
	Advect_c.log.txt\
	Advect_cpp.chem.txt\
	Advect_cpp.log.txt\
	advect.yaml\
	advection_bmi_cpp.dmp\
	advection_c.dmp\
	advection_cpp.dmp\
	gas_c.chem.txt\
	gas_c.dmp\
	gas_c.log.txt\
	gas_cpp.chem.txt\
	gas_cpp.dmp\
	gas_cpp.log.txt\
	gas_utility_c.txt\
	gas_utility_cpp.txt\
	SimpleAdvection_cpp.chem.txt\
	SimpleAdvection_cpp.log.txt\
	Species_c.chem.txt\
	species_c.dmp\
	Species_c.log.txt\
	Species_cpp.chem.txt\
	species_cpp.dmp\
	Species_cpp.log.txt\
	species_utility_c.txt\
	species_utility_cpp.txt\
	Units_InitialPhreeqc_1.chem.txt\
	Units_InitialPhreeqc_1.log.txt\
	Units_InitialPhreeqc_2.chem.txt\
	Units_InitialPhreeqc_2.log.txt\
	Units_utility.out\
	Units_Worker.chem.txt\
	Units_Worker.log.txt\
	utility_c.txt\
	utility_cpp.txt

if PHREEQCRM_TEST_FORTRAN

test_SOURCES+=\
	advection_f90.F90\
	gas_f90.F90\
	species_f90.F90\
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
