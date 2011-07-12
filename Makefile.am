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
	phreeqcpp/Dictionary.h\
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
	phreeqcpp/KineticsComp.cxx\
	phreeqcpp/KineticsComp.h\
	phreeqcpp/NameDouble.cxx\
	phreeqcpp/NameDouble.h\
	phreeqcpp/NumKeyword.cxx\
	phreeqcpp/NumKeyword.h\
	phreeqcpp/Parser.cxx\
	phreeqcpp/Parser.h\
	phreeqcpp/Phreeqc.cpp\
	phreeqcpp/Phreeqc.h\
	phreeqcpp/Phreeqc_class.h\
	phreeqcpp/PPassemblage.cxx\
	phreeqcpp/PPassemblage.h\
	phreeqcpp/PPassemblageComp.cxx\
	phreeqcpp/PPassemblageComp.h\
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
	phreeqcpp/phreeqc/advection.c\
	phreeqcpp/phreeqc/basic.c\
	phreeqcpp/phreeqc/basic.h\
	phreeqcpp/phreeqc/basicsubs.c\
	phreeqcpp/phreeqc/cl1.c\
	phreeqcpp/phreeqc/cvdense.c\
	phreeqcpp/phreeqc/cvdense.h\
	phreeqcpp/phreeqc/cvode.c\
	phreeqcpp/phreeqc/cvode.h\
	phreeqcpp/phreeqc/dense.c\
	phreeqcpp/phreeqc/dense.h\
	phreeqcpp/phreeqc/dw.c\
	phreeqcpp/phreeqc/global.h\
	phreeqcpp/phreeqc/global_structures.h\
	phreeqcpp/phreeqc/input.c\
	phreeqcpp/phreeqc/input.h\
	phreeqcpp/phreeqc/integrate.c\
	phreeqcpp/phreeqc/inverse.c\
	phreeqcpp/phreeqc/isotopes.c\
	phreeqcpp/phreeqc/kinetics.c\
	phreeqcpp/phreeqc/kinetics.h\
	phreeqcpp/phreeqc/mainsubs.c\
	phreeqcpp/phreeqc/model.c\
	phreeqcpp/phreeqc/nvector.c\
	phreeqcpp/phreeqc/nvector.h\
	phreeqcpp/phreeqc/nvector_serial.c\
	phreeqcpp/phreeqc/nvector_serial.h\
	phreeqcpp/phreeqc/output.c\
	phreeqcpp/phreeqc/output.h\
	phreeqcpp/phreeqc/p2c.h\
	phreeqcpp/phreeqc/p2clib.c\
	phreeqcpp/phreeqc/parse.c\
	phreeqcpp/phreeqc/phqalloc.c\
	phreeqcpp/phreeqc/phqalloc.h\
	phreeqcpp/phreeqc/phreeqc_files.c\
	phreeqcpp/phreeqc/phrqproto.h\
	phreeqcpp/phreeqc/phrqtype.h\
	phreeqcpp/phreeqc/pitzer.c\
	phreeqcpp/phreeqc/pitzer.h\
	phreeqcpp/phreeqc/pitzer_structures.c\
	phreeqcpp/phreeqc/pitzer_structures.h\
	phreeqcpp/phreeqc/prep.c\
	phreeqcpp/phreeqc/print.c\
	phreeqcpp/phreeqc/read.c\
	phreeqcpp/phreeqc/readtr.c\
	phreeqcpp/phreeqc/sit.c\
	phreeqcpp/phreeqc/smalldense.c\
	phreeqcpp/phreeqc/smalldense.h\
	phreeqcpp/phreeqc/spread.c\
	phreeqcpp/phreeqc/step.c\
	phreeqcpp/phreeqc/structures.c\
	phreeqcpp/phreeqc/sundialsmath.c\
	phreeqcpp/phreeqc/sundialsmath.h\
	phreeqcpp/phreeqc/sundialstypes.h\
	phreeqcpp/phreeqc/tally.c\
	phreeqcpp/phreeqc/tidy.c\
	phreeqcpp/phreeqc/transport.c\
	phreeqcpp/phreeqc/utilities.c\
	phreeqcpp/phreeqc/NA.h

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
