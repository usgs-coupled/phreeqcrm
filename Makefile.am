#AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/src -I$(top_srcdir)/src/phreeqcpp -I$(top_srcdir)/src/phreeqcpp/phreeqc
AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/src/phreeqcpp -I$(top_srcdir)/src/phreeqcpp/phreeqc

SUBDIRS=phreeqcpp

lib_LTLIBRARIES=libiphreeqc.la

##{{
##libiphreeqc_la_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/src/phreeqcpp -I$(top_srcdir)/src/phreeqcpp/phreeqc
##libiphreeqc_la_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/src -I$(top_srcdir)/src/phreeqcpp -I$(top_srcdir)/src/phreeqcpp/phreeqc
##}}

# library sources
libiphreeqc_la_SOURCES=\
	CVar.hxx\
	Debug.h\
	ErrorReporter.hxx\
	IPhreeqc.cpp\
	IPhreeqcLib.cpp\
	SelectedOutput.cpp\
	SelectedOutput.hxx\
	Var.c
	
libiphreeqc_la_LIBADD=\
	phreeqcpp/libphreeqcpp.la\
	phreeqcpp/phreeqc/libphreeqc.la	

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
