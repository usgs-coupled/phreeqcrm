ACLOCAL_AMFLAGS=-I m4
EXTRA_DIST=CMakeLists.txt
AM_CPPFLAGS=-I$(top_srcdir)/src/IPhreeqcPhast -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords

# libtool libraries to build
lib_LTLIBRARIES=libphreeqcrm.la

# library sources for libphreeqcrm.la
libphreeqcrm_la_SOURCES=\
	src/PhreeqcRM.cpp\
	src/RM_interface_C.cpp\
	src/IPhreeqcPhast/IPhreeqcPhast.cxx\
	src/IPhreeqcPhast/IPhreeqcPhast.h\
	src/IPhreeqcPhast/IPhreeqcPhastLib.cpp\
	src/IPhreeqcPhast/IPhreeqcPhastLib.h\
	src/IPhreeqcPhast/IPhreeqc/thread.h\
	src/IPhreeqcPhast/IPhreeqc/CSelectedOutput.cpp\
	src/IPhreeqcPhast/IPhreeqc/CSelectedOutput.hxx\
	src/IPhreeqcPhast/IPhreeqc/CVar.hxx\
	src/IPhreeqcPhast/IPhreeqc/Debug.h\
	src/IPhreeqcPhast/IPhreeqc/ErrorReporter.hxx\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqc.cpp\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqc.h\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqcCallbacks.h\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface_F.cpp\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface_F.h\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqcLib.cpp\
	src/IPhreeqcPhast/IPhreeqc/Var.c\
	src/IPhreeqcPhast/IPhreeqc/Var.h\
	src/IPhreeqcPhast/IPhreeqc/Version.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/advection.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/basicsubs.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cl1.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cvdense.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cvdense.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cvode.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cvode.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxKinetics.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxKinetics.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxMix.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxMix.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dense.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dense.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dumper.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dumper.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dw.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/gases.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/global_structures.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/input.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/integrate.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/inverse.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/isotopes.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/kinetics.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/mainsubs.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/model.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector_serial.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector_serial.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/parse.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/phqalloc.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/phqalloc.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/phrqtype.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/pitzer.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/pitzer_structures.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/prep.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/print.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/read.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/readtr.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/runner.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/runner.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sit.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/smalldense.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/smalldense.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/spread.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/step.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/structures.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialsmath.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialsmath.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialstypes.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/tally.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/tidy.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/transport.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/utilities.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartHandler.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartHandler.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartObject.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartObject.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/CurveObject.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/CurveObject.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ExchComp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ExchComp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Exchange.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Exchange.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasComp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasComp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasPhase.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasPhase.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolution.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolution.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolutionComp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolutionComp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/KineticsComp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/KineticsComp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NA.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NumKeyword.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NumKeyword.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PBasic.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PBasic.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_base.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_base.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io_output.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblage.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblage.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblageComp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblageComp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Parser.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Parser.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Phreeqc.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Phreeqc.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Pressure.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Pressure.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Reaction.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Reaction.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ReadClass.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SS.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SS.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SSassemblage.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SSassemblage.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SScomp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SScomp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SelectedOutput.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SelectedOutput.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Solution.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Solution.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SolutionIsotope.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SolutionIsotope.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBinList.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBinList.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Surface.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Surface.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceCharge.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceCharge.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceComp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceComp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/System.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/System.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Temperature.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Temperature.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Use.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Use.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/UserPunch.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/UserPunch.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Utils.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Utils.h


include_HEADERS=\
	$(top_srcdir)/src/IrmResult.h\
	$(top_srcdir)/src/PhreeqcRM.h\
	$(top_srcdir)/src/RM_interface_C.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqcPhast.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/IPhreeqcCallbacks.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/Var.h


if BUILD_FORTRAN

DEFS+=\
	-DTEST_FORTRAN

libphreeqcrm_la_SOURCES+=\
	src/RM_interface_F.h\
	src/RM_interface_F.cpp\
	src/IPhreeqcPhast/IPhreeqc/fwrap.cpp\
	src/IPhreeqcPhast/IPhreeqc/fwrap.h

include_HEADERS+=\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.f90.inc\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.f.inc

##install-data-hook:
##	@if test -n "$(FC_MODEXT)"; then \
##	  $(INSTALL_DATA) $(top_builddir)/*.$(FC_MODEXT) $(DESTDIR)$(includedir)/. ; \
##	fi

##uninstall-hook:
##	@if test -n "$(FC_MODEXT)"; then \
##	  $(RM) $(DESTDIR)$(includedir)/[pP][hH][rR][eE][eE][qQ][cC][rR][mM].$(FC_MODEXT) ; \
##	fi

##mostlyclean-local: clean-modules
##.PHONY: clean-modules
##clean-modules:
##	-test -z $(FC_MODEXT) || $(RM) *.$(FC_MODEXT)

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

SUBDIRS = Doxygen Tests
