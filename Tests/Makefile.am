EXTRA_DIST = CMakeLists.txt

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/IPhreeqcPhast -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords -I$(top_builddir)

TESTS = test
check_PROGRAMS = test

test_SOURCES = advection_c.c advection_cpp.cpp main.cpp species_c.c species_cpp.cpp
test_LDADD = $(top_builddir)/libphreeqcrm.la
test_LINK = $(CXXLINK)

CLEANFILES =\
	species_cpp.dmp\
	Species_cpp.log.txt\
	species_c.dmp\
	advection_c.dmp\
	Advect_c.chem.txt\
	Units_Worker.chem.txt\
	Species_c.chem.txt\
	Units_InitialPhreeqc_1.log.txt\
	species_utility_c.txt\
	Advect_cpp.chem.txt\
	Species_c.log.txt\
	utility_c.txt\
	advection_cpp.dmp\
	Units_InitialPhreeqc_2.log.txt\
	species_utility_cpp.txt\
	Units_Worker.log.txt\
	Units_utility.out\
	Advect_cpp.log.txt\
	Species_cpp.chem.txt\
	Units_InitialPhreeqc_2.chem.txt\
	Advect_c.log.txt\
	utility_cpp.txt\
	Units_InitialPhreeqc_1.chem.txt

if PHREEQCRM_TEST_FORTRAN

test_SOURCES+=\
	advection_f90.F90\
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
