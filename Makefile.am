EXTRA_DIST = fwrap2.cpp fwrap3.cpp

AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/src/phreeqcpp -I$(top_srcdir)/src/phreeqcpp/phreeqc

# libtool libraries to build
lib_LTLIBRARIES=libiphreeqc.la

# library sources for libiphreeqc.la
libiphreeqc_la_SOURCES=\
	CVar.hxx\
	Debug.h\
	ErrorReporter.hxx\
	IPhreeqc.cpp\
	IPhreeqcLib.cpp\
	SelectedOutput.cpp\
	SelectedOutput.hxx\
	Var.c\
	phreeqcpp/cxxKinetics.cxx\
	phreeqcpp/cxxKinetics.h\
	phreeqcpp/cxxMix.cxx\
	phreeqcpp/cxxMix.h\
	phreeqcpp/dumper.cpp\
	phreeqcpp/dumper.h\
	phreeqcpp/Exchange.cxx\
	phreeqcpp/Exchange.h\
	phreeqcpp/ExchComp.cxx\
	phreeqcpp/ExchComp.h\
	phreeqcpp/GasPhase.cxx\
	phreeqcpp/GasPhase.h\
	phreeqcpp/ISolution.cxx\
	phreeqcpp/ISolution.h\
	phreeqcpp/ISolutionComp.cxx\
	phreeqcpp/ISolutionComp.h\
	phreeqcpp/Keywords.cpp\
	phreeqcpp/Keywords.h\
	phreeqcpp/KineticsComp.cxx\
	phreeqcpp/KineticsComp.h\
	phreeqcpp/NameDouble.cxx\
	phreeqcpp/NameDouble.h\
	phreeqcpp/NumKeyword.cxx\
	phreeqcpp/NumKeyword.h\
	phreeqcpp/Parser.cxx\
	phreeqcpp/Parser.h\
	phreeqcpp/PBasic.cpp\
	phreeqcpp/PBasic.h\
	phreeqcpp/Phreeqc.cpp\
	phreeqcpp/Phreeqc.h\
	phreeqcpp/Phreeqc_class.h\
	phreeqcpp/PHRQ_base.cxx\
	phreeqcpp/PHRQ_base.h\
	phreeqcpp/PHRQ_io.cpp\
	phreeqcpp/PHRQ_io.h\
	phreeqcpp/PPassemblage.cxx\
	phreeqcpp/PPassemblage.h\
	phreeqcpp/PPassemblageComp.cxx\
	phreeqcpp/PPassemblageComp.h\
	phreeqcpp/Pressure.cxx\
	phreeqcpp/Pressure.h\
	phreeqcpp/Reaction.cxx\
	phreeqcpp/Reaction.h\
	phreeqcpp/ReadClass.cxx\
	phreeqcpp/runner.cpp\
	phreeqcpp/runner.h\
	phreeqcpp/Solution.cxx\
	phreeqcpp/Solution.h\
	phreeqcpp/SolutionIsotope.cxx\
	phreeqcpp/SolutionIsotope.h\
	phreeqcpp/SolutionIsotopeList.cxx\
	phreeqcpp/SolutionIsotopeList.h\
	phreeqcpp/SSassemblage.cxx\
	phreeqcpp/SSassemblage.h\
	phreeqcpp/SSassemblageSS.cxx\
	phreeqcpp/SSassemblageSS.h\
	phreeqcpp/StorageBin.cxx\
	phreeqcpp/StorageBin.h\
	phreeqcpp/StorageBinList.cpp\
	phreeqcpp/StorageBinList.h\
	phreeqcpp/Surface.cxx\
	phreeqcpp/Surface.h\
	phreeqcpp/SurfaceCharge.cxx\
	phreeqcpp/SurfaceCharge.h\
	phreeqcpp/SurfaceComp.cxx\
	phreeqcpp/SurfaceComp.h\
	phreeqcpp/System.cxx\
	phreeqcpp/System.h\
	phreeqcpp/Temperature.cxx\
	phreeqcpp/Temperature.h\
	phreeqcpp/Utils.cxx\
	phreeqcpp/Utils.h\
	phreeqcpp/phreeqc/advection.cpp\
	phreeqcpp/phreeqc/basicsubs.cpp\
	phreeqcpp/phreeqc/cl1.cpp\
	phreeqcpp/phreeqc/cvdense.cpp\
	phreeqcpp/phreeqc/cvdense.h\
	phreeqcpp/phreeqc/cvode.cpp\
	phreeqcpp/phreeqc/cvode.h\
	phreeqcpp/phreeqc/dense.cpp\
	phreeqcpp/phreeqc/dense.h\
	phreeqcpp/phreeqc/dw.cpp\
	phreeqcpp/phreeqc/global_structures.h\
	phreeqcpp/phreeqc/input.cpp\
	phreeqcpp/phreeqc/integrate.cpp\
	phreeqcpp/phreeqc/inverse.cpp\
	phreeqcpp/phreeqc/isotopes.cpp\
	phreeqcpp/phreeqc/kinetics.cpp\
	phreeqcpp/phreeqc/mainsubs.cpp\
	phreeqcpp/phreeqc/model.cpp\
	phreeqcpp/phreeqc/NA.h
	phreeqcpp/phreeqc/nvector.cpp\
	phreeqcpp/phreeqc/nvector.h\
	phreeqcpp/phreeqc/nvector_serial.cpp\
	phreeqcpp/phreeqc/nvector_serial.h\
	phreeqcpp/phreeqc/parse.cpp\
	phreeqcpp/phreeqc/phqalloc.cpp\
	phreeqcpp/phreeqc/phqalloc.h\
	phreeqcpp/phreeqc/phrqproto.h\
	phreeqcpp/phreeqc/phrqtype.h\
	phreeqcpp/phreeqc/pitzer.cpp\
	phreeqcpp/phreeqc/pitzer.h\
	phreeqcpp/phreeqc/pitzer_structures.cpp\
	phreeqcpp/phreeqc/prep.cpp\
	phreeqcpp/phreeqc/print.cpp\
	phreeqcpp/phreeqc/read.cpp\
	phreeqcpp/phreeqc/readtr.cpp\
	phreeqcpp/phreeqc/sit.cpp\
	phreeqcpp/phreeqc/smalldense.cpp\
	phreeqcpp/phreeqc/smalldense.h\
	phreeqcpp/phreeqc/spread.cpp\
	phreeqcpp/phreeqc/step.cpp\
	phreeqcpp/phreeqc/structures.cpp\
	phreeqcpp/phreeqc/sundialsmath.cpp\
	phreeqcpp/phreeqc/sundialsmath.h\
	phreeqcpp/phreeqc/sundialstypes.h\
	phreeqcpp/phreeqc/tally.cpp\
	phreeqcpp/phreeqc/tidy.cpp\
	phreeqcpp/phreeqc/transport.cpp\
	phreeqcpp/phreeqc/utilities.cpp

include_HEADERS=\
	$(top_srcdir)/include/IPhreeqc.h\
	$(top_srcdir)/include/IPhreeqc.hpp\
	$(top_srcdir)/include/IPhreeqcCallbacks.h\
	$(top_srcdir)/include/Var.h


if BUILD_FORTRAN

libiphreeqc_la_SOURCES+=\
	fwrap.h\
	fwrap.cpp\
	IPhreeqcF.f

include_HEADERS+=\
	$(top_srcdir)/include/IPhreeqc.f90.inc\
	$(top_srcdir)/include/IPhreeqc.f.inc

endif

COMPILE = $(CXXCOMPILE)
LTCOMPILE = $(LTCXXCOMPILE)
