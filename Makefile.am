ACLOCAL_AMFLAGS=-I m4

PHREEQCRM_EXTRA_DIST=\
	CMakeLists.txt\
	phreeqcrm.cmake.in\
	PhreeqcRMConfig.cmake.in\
	swig/CMakeLists.txt\
	swig/numpy.i\
	swig/PhreeqcRM.i\
	swig/pyfragments.swg.in\
	swig/python/AdvectBMI_py.py\
	swig/python/BMIPhreeqcRM_docstrings.swg\
	swig/python/CMakeLists.txt\
	swig/python/conftest.py\
	swig/python/constants.py\
	swig/python/PhreeqcRM_docstrings.swg\
	swig/python/SimpleAdvect.py\
	swig/python/test_basic.py\
	swig/python/test_get_value.py\
	swig/python/test_init.py\
	swig/python/test_irf.py\
	swig/python/test_set_value.py\
	swig/python/TestAllMethods_py.py\
	swig/python/WriteYAMLFile_py.py\
	swig/python/yamlphreeqcrm.py\
	swig/python/YAMLPhreeqcRM.swg

EXTRA_DIST=$(PHREEQCRM_EXTRA_DIST)

AM_CPPFLAGS=-I$(top_srcdir)/src/IPhreeqcPhast -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common -I$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords

AM_LDFLAGS=-release 3.0.0
# AM_LDFLAGS= -no-undefined -release 3.0.0

pkgconfigdir = $(datadir)/pkgconfig
pkgconfig_DATA = phreeqcrm.pc

# libtool libraries to build
lib_LTLIBRARIES=libphreeqcrm.la

# library sources for libphreeqcrm.la
libphreeqcrm_la_SOURCES=\
	src/BMI_interface_C.cpp\
	src/BMI_interface_C.h\
	src/BMI_interface_F.cpp\
	src/BMI_interface_F.h\
	src/bmi.hxx\
	src/BMIPhreeqcRM.cpp\
	src/BMIPhreeqcRM.h\
	src/BMIVariant.cpp\
	src/BMIVariant.h\
	src/IPhreeqcPhast/IPhreeqc/CSelectedOutput.cpp\
	src/IPhreeqcPhast/IPhreeqc/CSelectedOutput.hxx\
	src/IPhreeqcPhast/IPhreeqc/CVar.hxx\
	src/IPhreeqcPhast/IPhreeqc/Debug.h\
	src/IPhreeqcPhast/IPhreeqc/ErrorReporter.hxx\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface_F.cpp\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface_F.h\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqc.cpp\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqc.h\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqcCallbacks.h\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqcLib.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/advection.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/basicsubs.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/cl1.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/Parser.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/Parser.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_base.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_base.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_exports.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_io.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_io.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/phrqtype.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/Utils.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/Utils.h\
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
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Dictionary.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Dictionary.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dumper.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/dumper.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Exchange.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Exchange.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ExchComp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ExchComp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasComp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasComp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/gases.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasPhase.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/GasPhase.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/global_structures.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/input.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/integrate.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/inverse.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolution.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolution.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolutionComp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ISolutionComp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/isotopes.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/kinetics.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/KineticsComp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/KineticsComp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/mainsubs.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/model.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NA.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NumKeyword.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NumKeyword.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector_serial.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector_serial.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/nvector.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/parse.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PBasic.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PBasic.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/phqalloc.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/phqalloc.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Phreeqc.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Phreeqc.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PHRQ_io_output.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/pitzer_structures.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/pitzer.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblage.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblage.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblageComp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PPassemblageComp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/prep.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Pressure.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Pressure.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/print.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Reaction.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Reaction.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/read.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/ReadClass.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/readtr.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/runner.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/runner.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SelectedOutput.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SelectedOutput.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Serializer.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Serializer.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sit.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/smalldense.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/smalldense.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Solution.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Solution.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SolutionIsotope.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SolutionIsotope.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/spread.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SS.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SS.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SSassemblage.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SSassemblage.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SScomp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SScomp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/step.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBinList.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBinList.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/structures.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialsmath.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialsmath.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/sundialstypes.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Surface.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Surface.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceCharge.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceCharge.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceComp.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/SurfaceComp.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/System.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/System.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/tally.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Temperature.cxx\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Temperature.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/tidy.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/transport.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Use.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/Use.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/UserPunch.cpp\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/UserPunch.h\
	src/IPhreeqcPhast/IPhreeqc/phreeqcpp/utilities.cpp\
	src/IPhreeqcPhast/IPhreeqc/thread.h\
	src/IPhreeqcPhast/IPhreeqc/Var.c\
	src/IPhreeqcPhast/IPhreeqc/Var.h\
	src/IPhreeqcPhast/IPhreeqc/Version.h\
	src/IPhreeqcPhast/IPhreeqcPhast.cxx\
	src/IPhreeqcPhast/IPhreeqcPhast.h\
	src/IPhreeqcPhast/IPhreeqcPhastLib.cpp\
	src/IPhreeqcPhast/IPhreeqcPhastLib.h\
	src/PhreeqcRM.cpp\
	src/RM_interface_C.cpp\
	src/RM_interface_F.cpp\
	src/RM_interface_F.h\
	src/RMVARS.h\
	src/VarManager.cpp\
	src/VarManager.h\
	src/YAML_interface_C.cpp\
	src/YAML_interface_C.h\
	src/YAML_interface_F.cpp\
	src/YAML_interface_F.h\
	src/YAMLPhreeqcRM.cpp\
	src/YAMLPhreeqcRM.h


include_HEADERS=\
	$(top_builddir)/irm_dll_export.h\
	$(top_srcdir)/src/BMI_interface_C.h\
	$(top_srcdir)/src/BMI_interface_F.h\
	$(top_srcdir)/src/bmi.hxx\
	$(top_srcdir)/src/BMIPhreeqcRM.h\
	$(top_srcdir)/src/BMIVariant.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/IPhreeqc.hpp\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/IPhreeqcCallbacks.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/Parser.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_base.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_exports.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/PHRQ_io.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/common/phrqtype.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/NameDouble.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/PhreeqcKeywords/Keywords.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/StorageBin.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/phreeqcpp/System.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqc/Var.h\
	$(top_srcdir)/src/IPhreeqcPhast/IPhreeqcPhast.h\
	$(top_srcdir)/src/IrmResult.h\
	$(top_srcdir)/src/PhreeqcRM.h\
	$(top_srcdir)/src/RM_interface_C.h\
	$(top_srcdir)/src/RMVARS.h\
	$(top_srcdir)/src/VarManager.h\
	$(top_srcdir)/src/YAML_interface_C.h\
	$(top_srcdir)/src/YAML_interface_F.h\
	$(top_srcdir)/src/YAMLPhreeqcRM.h

FCDEFS=
if PHREEQCRM_TEST_FORTRAN
DEFS+=-DTEST_FORTRAN
FCDEFS+=-DTEST_FORTRAN
endif

# Docs
#
sourcedir=$(docdir)/src

dist_source_DATA=$(SRC)

SRC=\
	src/BMI_interface.F90\
	src/BMI_not_implemented.inc\
	src/IPhreeqcPhast/IPhreeqc/IPhreeqc_interface.F90\
	src/RM_interface.F90\
	src/YAML_interface.F90


# We must override the default PPFCCOMPILE since ifort can't
# handle the FC_FUNC define correctly so $(DEFS) is replaced
# by $(FCDEFS)
#
if PHREEQCRM_WITH_YAMLCPP
FCDEFS+=-DUSE_YAML
endif
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
        $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)
LTPPFCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=compile $(FC) $(FCDEFS) \
        $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
        $(AM_FCFLAGS) $(FCFLAGS)

SUBDIRS = database doc Doxygen Tests
