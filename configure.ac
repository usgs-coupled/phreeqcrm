<#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.61])
AC_INIT([PhreeqcRM], [0.1.0-9067], [dlpark@usgs.gov])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign no-define subdir-objects])
AC_CONFIG_SRCDIR([src/PhreeqcRM.cpp])

# enable silent rules when available (automake 1.11 or later).
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CANONICAL_HOST

AC_DEFINE([SWIG_SHARED_OBJ], , [Must be set for PhreeqcRM])
AC_DEFINE([USE_PHRQ_ALLOC], , [Must be set for PhreeqcRM])
AC_DEFINE([NDEBUG], , [Must be set for PhreeqcRM])

# check for --with-mpi
# If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
# If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
# If --with-mpi=no is used, use a standard C compiler instead.
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
    [compile with MPI (parallelization) support. If none is found,
    MPI is not used. Default: no])
],,[with_mpi=no])


# Checks for programs.
AX_PROG_CC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
  use_mpi=no
  if test x"$with_mpi" = xyes; then
    AC_MSG_ERROR([MPI compiler requested, but couldn't use MPI.])
  fi
])
AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
  use_mpi=no
  if test x"$with_mpi" = xyes; then
    AC_MSG_ERROR([MPI compiler requested, but couldn't use MPI.])
  fi
])
AC_PROG_AWK
AC_PROG_CPP


# Check if the fortran test should be included
AC_MSG_CHECKING([if Fortran test is added])
AC_ARG_ENABLE([fortran-test],
     [AS_HELP_STRING([--enable-fortran-test],[test Fortran interface @<:@default=no@:>@])],
     [IRM_FORTRAN=$enableval],
     [IRM_FORTRAN=no])

if test "X$IRM_FORTRAN" = "Xyes"; then
  AC_MSG_RESULT(yes)
  AX_PROG_FC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
    use_mpi=no
    if test x"$with_mpi" = xyes; then
      AC_MSG_ERROR([MPI compiler requested, but couldn't use MPI.])
    fi
  ])
  AC_FC_LIBRARY_LDFLAGS
  AC_FC_WRAPPERS 
  AC_FC_MAIN
  m4_ifdef([AC_FC_MODULE_EXTENSION], [AC_FC_MODULE_EXTENSION], [
    AX_F90_MODULE_EXTENSION
    FC_MODEXT=mod
    if test x$ax_cv_f90_modext != xunknown ; then
      FC_MODEXT=$ax_cv_f90_modext
    fi
    AC_SUBST(FC_MODEXT)
  ])
  AC_SUBST(F77, $FC)
  AC_SUBST(FFLAGS, $FCFLAGS) 
else
  AC_MSG_RESULT(no)
fi
AM_CONDITIONAL([PHREEQCRM_TEST_FORTRAN], [test "X$IRM_FORTRAN" = "Xyes"])

AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT


# Checks for libraries.
if test x"$with_mpi" = xyes; then
  enable_openmp=no
fi

if test x"$enable_openmp" != xno; then
  AC_LANG_PUSH(C++)
  AC_OPENMP
  AC_LANG_POP(C++)
  AC_SUBST(AM_CXXFLAGS, $OPENMP_CXXFLAGS)
fi

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h memory.h stddef.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD

# libc functions
AC_CHECK_FUNCS([memmove], , AC_ERROR([memmove not found in libc]))
AC_CHECK_FUNCS([memset],  , AC_ERROR([memset not found in libc]))
AC_CHECK_FUNCS([strchr],  , AC_ERROR([strchr not found in libc]))
AC_CHECK_FUNCS([strcspn], , AC_ERROR([strcspn not found in libc]))
AC_CHECK_FUNCS([strtol],  , AC_ERROR([strtol not found in libc]))

# libm functions
AC_CHECK_FUNCS([floor],   , AC_CHECK_LIB(m, floor, , AC_ERROR([cannot find floor])) )
AC_CHECK_FUNCS([pow],     , AC_CHECK_LIB(m, pow,   , AC_ERROR([cannot find pow])) )
AC_CHECK_FUNCS([sqrt],    , AC_CHECK_LIB(m, sqrt,  , AC_ERROR([cannot find sqrt])) )

case $host_os in
     solaris*)
	CXXFLAGS="$CXXFLAGS -library=stlport4"
	LDFLAGS="$LDFLAGS -library=stlport4"
	;;
esac


# isfinite
AC_LANG_PUSH(C++)
AC_CACHE_CHECK([for isfinite], ac_cv_isfinite,
 [AC_TRY_LINK([#include <math.h>],
 [double x; int y; y = isfinite(x);],
 ac_cv_isfinite=yes,
 ac_cv_isfinite=no
)])
if test x"$ac_cv_isfinite" = x"yes"; then
  AC_DEFINE(HAVE_ISFINITE, [], [Has isfinite])
fi
AC_LANG_POP(C++)

# finite
AC_LANG_PUSH(C++)
AC_CACHE_CHECK([for finite], ac_cv_finite,
 [AC_TRY_LINK([#include <math.h>],
 [double x; int y; y = finite(x);],
 ac_cv_finite=yes,
 ac_cv_finite=no
)])
if test x"$ac_cv_finite" = x"yes"; then
  AC_DEFINE(HAVE_FINITE, [], [Has finite])
fi
AC_LANG_POP(C++)

# isnan
AC_CHECK_FUNCS([isnan],    , AC_CHECK_LIB(m, isnan))

## Taken from hdf5-1.8.14
## ----------------------------------------------------------------------
## Fake --with-xxx option to allow us to create a help message for the
## following --with-xxx options which can take either a =DIR or =INC,LIB
## specifier.
##
AC_ARG_WITH([fnord],
  [
 For the following --with-xxx options, you can specify where the header
 files and libraries are in two different ways:

    --with-xxx=INC,LIB - Specify individually the include directory and
                         library directory separated by a comma
    --with-xxx=DIR     - Specify only the directory which contains the
                         include/ and lib/ subdirectories
  ])

#  ----------------------------------------------------------------------
#  Is the zlib present? It has a header file `zlib.h' and a library
#  `-lz' and their locations might be specified with the `--with-zlib'
#  command-line switch. The value is an include path and/or a library path.
#  If the library path is specified then it must be preceded by a comma.
#
AC_ARG_WITH([zlib],
            [AC_HELP_STRING([--with-zlib=DIR],
                            [Use the zlib library [default=no]])],,withval=no)

case $withval in
  no)
    ;;
  yes)
    HAVE_ZLIB="yes"
    AC_CHECK_HEADERS(zlib.h, , AC_MSG_ERROR([zlib.h not found]))
    AC_CHECK_LIB(z, gzwrite, , AC_MSG_ERROR([couldn't find zlib library]))
    ;;
  *)
    HAVE_ZLIB="yes"
    case "$withval" in
      *,*)
        zlib_inc="`echo $withval |cut -f1 -d,`"
        zlib_lib="`echo $withval |cut -f2 -d, -s`"
        ;;
      *)
        if test -n "$withval"; then
          zlib_inc="$withval/include"
          zlib_lib="$withval/lib"
        fi
        ;;
    esac

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"

    if test -n "$zlib_inc"; then
      CPPFLAGS="$CPPFLAGS -I$zlib_inc"
    fi

    AC_CHECK_HEADERS(zlib.h,, AC_MSG_ERROR([zlib.h not found]))

    if test -n "$zlib_lib"; then
      LDFLAGS="$LDFLAGS -L$zlib_lib"
    fi

    AC_CHECK_LIB(z, gzwrite,, LDFLAGS="$saved_LDFLAGS"; unset HAVE_ZLIB)

    if test -z "$HAVE_ZLIB"; then
      AC_MSG_ERROR(couldn't find zlib library)
    fi
    ;;
esac
if test -n "$HAVE_ZLIB"; then
  AC_DEFINE([USE_GZ], [], [Set if zlib is present])
fi

AM_CONDITIONAL([BUILD_OPENMP], [test "X$enable_openmp" != "Xno"])
AM_CONDITIONAL([BUILD_MPI], [test "X$with_mpi" = "Xyes"])

if test "X$with_mpi" != "Xno"; then
  AC_DEFINE([USE_MPI], [], [Set to use MPI])
  PHAST_SUFFIX=mpi
else
  if test "X$enable_openmp" != "Xno"; then
    AC_DEFINE([USE_OPENMP], [], [Set to to use OpenMP])
  fi
fi

AC_CONFIG_FILES([
  Makefile
  database/Makefile
  doc/Makefile
  Doxygen/Makefile
  Tests/Makefile
  Tests/advect.pqi
  Tests/phreeqc.dat
  Tests/units.pqi
])

AC_OUTPUT
echo "*************************************"
echo "Summary:"
echo "C++ Compiler                 = $CXX"
echo "CPPFLAGS                     = $CPPFLAGS"
echo "CXXFLAGS                     = $CXXFLAGS"
echo "Fortran Compiler             = $FC"
echo "Fortran Libs                 = $FCLIBS"
echo "LDFLAGS                      = $LDFLAGS"
if test "X$enable_openmp" != "Xno"; then
  echo "Enable OpenMP                = yes"
else
  echo "Enable OpenMP                = no"
fi
echo "Enable MPI                   = $with_mpi"
echo "*************************************"
