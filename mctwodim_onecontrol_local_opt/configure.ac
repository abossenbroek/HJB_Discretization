AC_PREREQ(2.59)
AC_INIT([mctwodim_onecontrol_local_opt], [0.01], [anton.bossenbroek@me.com])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/markovchain.cpp])
AC_CONFIG_HEADER([config.h])

AC_GNU_SOURCE

AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AM_PROG_AS
AM_PROG_CC_C_O

AC_LANG([C++])

ACLOCAL_AMFLAGS="-I m4 $ACLOCAL_AMFLAGS"
AC_SUBST(ACLOCAL_AMFLAGS)

# Checks for libraries.
AX_BOOST_BASE([1.37])


AC_CACHE_SAVE

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

AM_PATH_GSL(1.11, , [AC_MSG_ERROR(Please ensure that the proper GSL version is installed)])


if test "$ac_cv_c_inline" != no ; then
  AC_DEFINE(HAVE_INLINE,1, [The compiler supports inline])
  AC_SUBST(HAVE_INLINE)
fi


AC_CHECK_LIB(lis, lis_matrix_create,
		[AC_CHECK_HEADERS([/usr/local/include/lis.h],
					[CPPFLAGS="${CPPFLAGS} -I/usr/local/include"; LDFLAGS="$LDFLAGS -llis"],
					[AC_MSG_ERROR(Could not find required header lis)])],
					[AC_MSG_ERROR(Could not find required library lis)])

# Check if optimization flags must be used.
AC_ARG_WITH([optimize-flags],
	AS_HELP_STRING([--with-optimize-flags], [use compiler optimization flags (default is yes)]),
	[
    if test "$withval" = "no"; then
			want_optimizer="no"
    elif test "$withval" = "yes"; then
      want_optimizer="yes"
    else
	    echo "Unkown switch for --optimize-flags=(yes|no)!" ; exit -1
	fi
    ],
    [want_optimizer="yes"])

if test "x$want_optimizer" = "xyes"; then
  AX_CHECK_COMPILER_FLAGS([-fvect-cost-model], 
                        [CPPFLAGS="${CPPFLAGS} -fvect-cost-model"] ,)

  AX_CHECK_COMPILER_FLAGS([-ftree-vectorize], 
                        [CPPFLAGS="${CPPFLAGS} -ftree-vectorize"] ,)

  AX_CHECK_COMPILER_FLAGS([-mfpmath=sse -msse -msse2 -funroll-loops], 
                        [CPPFLAGS="${CPPFLAGS} -mfpmath=sse -msse -msse2 -funroll-loops"] ,)

  AX_CHECK_COMPILER_FLAGS([-msse3], [CPPFLAGS="${CPPFLAGS} -msse3"], )
  AX_CHECK_COMPILER_FLAGS([-O3], [CPPFLAGS="${CPPFLAGS} -O3"], )
else
  AX_CHECK_COMPILER_FLAGS([-O0], [CPPFLAGS="${CPPFLAGS} -O0 -ggdb"], )
fi

AX_CHECK_COMPILER_FLAGS([-fopenmp], [CPPFLAGS="${CPPFLAGS} -fopenmp"], AC_MSG_ERROR(Require -fopenmp))
# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile
				 ])
AC_OUTPUT

# vim: 
