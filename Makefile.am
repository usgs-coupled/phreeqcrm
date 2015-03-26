ACLOCAL_AMFLAGS=-I m4
AM_CPPFLAGS=-I$(top_srcdir)/IPhreeqcPhast -I$(top_srcdir)/IPhreeqcPhast/IPhreeqc -I$(top_srcdir)/IPhreeqcPhast/IPhreeqc/phreeqcpp -I$(top_srcdir)/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords

# libtool libraries to build
lib_LTLIBRARIES=libphreeqcrm.la

# library sources for libphreeqcrm.la
libphreeqcrm_la_SOURCES=\
	PhreeqcRM.cpp\
	RM_interface_C.cpp\
	IPhreeqcPhast/IPhreeqcPhast.cxx\
	IPhreeqcPhast/IPhreeqcPhast.h\
	IPhreeqcPhast/IPhreeqcPhastLib.cpp\
	IPhreeqcPhast/IPhreeqcPhastLib.h\
	IPhreeqcPhast/IPhreeqc/thread.h\
	IPhreeqcPhast/IPhreeqc/CSelectedOutput.cpp\
	IPhreeqcPhast/IPhreeqc/CSelectedOutput.hxx\
	IPhreeqcPhast/IPhreeqc/CVar.hxx\
	IPhreeqcPhast/IPhreeqc/Debug.h\
	IPhreeqcPhast/IPhreeqc/ErrorReporter.hxx\
	IPhreeqcPhast/IPhreeqc/IPhreeqc.cpp\
	IPhreeqcPhast/IPhreeqc/IPhreeqc.h\
	IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp\
	IPhreeqcPhast/IPhreeqc/IPhreeqcCallbacks.h\
	IPhreeqcPhast/IPhreeqc/IPhreeqc_interface_F.cpp\
	IPhreeqcPhast/IPhreeqc/IPhreeqc_interface_F.h\
	IPhreeqcPhast/IPhreeqc/IPhreeqcLib.cpp\
	IPhreeqcPhast/IPhreeqc/Var.c\
	IPhreeqcPhast/IPhreeqc/Var.h\
	IPhreeqcPhast/IPhreeqc/Version.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/advection.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/basicsubs.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/cl1.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/cvdense.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/cvdense.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/cvode.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/cvode.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxKinetics.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxKinetics.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxMix.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/cxxMix.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/dense.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/dense.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/dumper.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/dumper.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/dw.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/gases.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/global_structures.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/input.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/integrate.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/inverse.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/isotopes.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/kinetics.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/mainsubs.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/model.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector_serial.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector_serial.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/parse.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/phqalloc.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/phqalloc.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/phrqtype.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/pitzer.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/pitzer_structures.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/prep.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/print.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/read.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/readtr.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/runner.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/runner.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/sit.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/smalldense.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/smalldense.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/spread.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/step.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/structures.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialsmath.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialsmath.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialstypes.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/tally.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/tidy.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/transport.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/utilities.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartHandler.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartHandler.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartObject.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/ChartObject.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/CurveObject.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/CurveObject.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/ExchComp.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/ExchComp.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Exchange.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Exchange.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/GasComp.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/GasComp.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/GasPhase.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/GasPhase.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolution.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolution.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolutionComp.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolutionComp.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/KineticsComp.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/KineticsComp.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/NA.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/NumKeyword.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/NumKeyword.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PBasic.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PBasic.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_base.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_base.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io_output.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblage.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblage.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblageComp.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblageComp.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Parser.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Parser.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Phreeqc.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Phreeqc.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Pressure.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Pressure.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Reaction.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Reaction.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/ReadClass.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SS.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SS.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SSassemblage.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SSassemblage.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SScomp.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SScomp.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SelectedOutput.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SelectedOutput.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Solution.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Solution.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SolutionIsotope.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SolutionIsotope.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBinList.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBinList.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Surface.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Surface.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceCharge.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceCharge.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceComp.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceComp.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/System.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/System.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Temperature.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Temperature.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Use.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Use.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/UserPunch.cpp\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/UserPunch.h\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Utils.cxx\
	IPhreeqcPhast/IPhreeqc/phreeqcpp/Utils.h


include_HEADERS=\
	$(top_srcdir)/IrmResult.h\
	$(top_srcdir)/PhreeqcRM.h\
	$(top_srcdir)/RM_interface_C.h\
	$(top_srcdir)/IPhreeqcPhast/IPhreeqcPhast.h\
	$(top_srcdir)/IPhreeqcPhast/IPhreeqc/IPhreeqc.h\
	$(top_srcdir)/IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp\
	$(top_srcdir)/IPhreeqcPhast/IPhreeqc/IPhreeqcCallbacks.h\
	$(top_srcdir)/IPhreeqcPhast/IPhreeqc/Var.h


if BUILD_FORTRAN

DEFS+=\
	-DTEST_FORTRAN

libphreeqcrm_la_SOURCES+=\
	RM_interface_F.h\
	RM_interface_F.cpp\
	IPhreeqcPhast/IPhreeqc/fwrap.cpp\
	IPhreeqcPhast/IPhreeqc/fwrap.h

include_HEADERS+=\
	$(top_srcdir)/IPhreeqcPhast/IPhreeqc/IPhreeqc.f90.inc\
	$(top_srcdir)/IPhreeqcPhast/IPhreeqc/IPhreeqc.f.inc

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

SUBDIRS = Tests
