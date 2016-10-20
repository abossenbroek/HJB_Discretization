AC_PREREQ(2.59)
AC_INIT([markovchain], [0.02], [anton.bossenbroek@me.com])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/markovchain.cpp])
AC_CONFIG_HEADER([config.h])

AC_GNU_SOURCE

AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_LIBTOOL
AM_PROG_AS
AM_PROG_CC_C_O

AC_LANG([C++])

ACLOCAL_AMFLAGS="-I m4 $ACLOCAL_AMFLAGS"
AC_SUBST(ACLOCAL_AMFLAGS)

# Checks for libraries.

AC_ARG_WITH([accelerate], 
	AS_HELP_STRING([--with-accelerate], [specify if accelerate is to be used for blas (default is no)]),
	[
	 if test $withval = "no"; then
		 want_accelerate="no"
	 elif test $withval = "yes"; then
		 want_accelerate="yes"
	 else
		 echo "Unkown value for --with-accelerate=(yes|no) !" ; exit -1
	 fi
	 ], [want_accelerate="no"])

if test "x$want_accelerate" = "xyes"; then
	AX_CHECK_COMPILER_FLAGS([-framework Accelerate], 
													 [AC_DEFINE(HAVE_ACCELERATE, 1, [The compiler supports the accelerate framework])
													 LDFLAGS="${LDFLAGS} -framework Accelerate"] , 
													 [-framework Accerlate is not supported by the compiler; exit -1])
else
	ACX_BLAS(,[AC_MSG_ERROR(BLAS could not be found)])
fi

AC_ARG_WITH([luajit],
	AS_HELP_STRING([--with-luajit], [specify if luajit should be used (default no)]),
		[
		 if test $withval = "yes"; then
			 want_luajit="yes"
		 elif test $withval = "no" ; then
			 want_luajit="no"
		 else
		 echo "Unkown value for --with-luajit=(yes|no) !" ; exit -1
	 fi
	 ], [want_luajit="no"])

if test "x$want_luajit" = "xyes"; then
	AC_SEARCH_LIBS([luaJIT_run], [lua], [AC_DEFINE(HAVE_LUAJIT, 1, [AC_MSG_ERROR(The luajit library was found)])],
		  [AC_MSG_ERROR(The luaJIT_run function cannot be found in a library)])
	AC_EGREP_HEADER("jit", /usr/local/include/lualib.h, [LUAHEADER="/usr/local/inlude"],
			[AC_EGREP_HEADER("jit", /usr/include/lua5.1/lualib.h, [LUAHEADER="/usr/inlude/lua5.1"],
				[AC_EGREP_HEADER("jit", /opt/local/include/lualib.h, [LUAHEADER="/opt/local/include"],
					[AC_MSG_ERROR(No jit function was found in a header)])])])
else
	# Check LUA stuff.
	AC_CHECK_LIB(lua5.1, luaL_openlibs,
    [AC_CHECK_HEADERS([/usr/include/lua5.1/lua.h],
        [CPPFLAGS="${CPPFLAGS} -I/usr/include/lua5.1"; LDFLAGS="${LDFLAGS} -llua5.1"],
		[AC_CHECK_HEADERS([/opt/local/include/lua.h],
        		[CPPFLAGS="${CPPFLAGS} -I/opt/local/include"; LDFLAGS="${LDFLAGS} -llua5.1"],
                  [AC_CHECK_HEADERS(lua.h, ,
                    [echo "We need lua headers files"; exit -1])])])],
    [AC_CHECK_LIB(lua, luaL_openlibs,
      [AC_CHECK_HEADERS([/usr/include/lua5.1/lua.h],
        [CPPFLAGS="${CPPFLAGS} -I/usr/include/lua5.1"; LDFLAGS="${LDFLAGS} -llua"],
		[AC_CHECK_HEADERS([/opt/local/include/lua.h],
        		[CPPFLAGS="${CPPFLAGS} -I/opt/local/include"; LDFLAGS="${LDFLAGS} -llua"],
                  [AC_CHECK_HEADERS(lua.h, ,
                    [echo "We need lua headers files"; exit -1])])])])])
fi

AX_BOOST_BASE([1.37])



AC_CACHE_SAVE

AM_PATH_GSL(1.11, , [AC_MSG_ERROR(Please ensure that the proper GSL version is installed)])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

if test "$ac_cv_c_inline" != no ; then
  AC_DEFINE(HAVE_INLINE,1, [The compiler supports inline])
  AC_SUBST(HAVE_INLINE)
fi


# Check if vector range checking is required.
AC_ARG_WITH([gsl-range-check],
	AS_HELP_STRING([--with-gsl-range-check], [turn on GSL vector and matrix range check (default is no)]),
	[
    if test "$withval" = "no"; then
			want_gsl_range_check="no"
    elif test "$withval" = "yes"; then
      want_gsl_range_check="yes"
    else
	    echo "Unkown switch for --gsl-range-check=(yes|no) !" ; exit -1
	fi
    ],
    [want_gsl_range_check="no"])

if test "x$want_gsl_range_check" = "xno"; then
  AC_DEFINE(GSL_RANGE_CHECK_OFF, 1, [Do not check for vector and matrix dimension matches]) 
fi

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
                        [CPPFLAGS="${CPPFLAGS} -fvect-cost-model"] , 
                        [echo "Try using gcc 4.2 or higher."; exit -1])

  AX_CHECK_COMPILER_FLAGS([-ftree-vectorize], 
                        [CPPFLAGS="${CPPFLAGS} -ftree-vectorize"] , 
                        [echo "Try using gcc 4.2 or higher."; exit -1])

  AX_CHECK_COMPILER_FLAGS([-mfpmath=sse -msse -msse2 -funroll-loops], 
                        [CPPFLAGS="${CPPFLAGS} -mfpmath=sse -msse -msse2 -funroll-loops"] , 
                        [echo "Some optimization flags are not supported."])

  AX_CHECK_COMPILER_FLAGS([-msse3], [CPPFLAGS="${CPPFLAGS} -O3 -msse3"], )
else
  AX_CHECK_COMPILER_FLAGS([-O0], [CPPFLAGS="${CPPFLAGS} -O0 -ggdb"], )
fi

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/markovsolverlib/Makefile
				 ])
AC_OUTPUT

# vim: 
