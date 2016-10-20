markovchain_src = markovchain.cpp \
					markovchain.h 
			
problemcalloption_src = problemcalloption.cpp \
												problemcalloption.h

problemlua_src = problemlua.h \
								 problemlua.cpp

solverexplicitgm_src = solverexplicitgm.h \
										 solverexplicitgm.cpp

solverimplicitbm_src = solverimplicitbm.h \
										 solverimplicitbm.cpp

problemmerton_src = problemmerton.h \
										 problemmerton.cpp

problemindiffpricing_src = problemindiffpricing.h \
													 problemindiffpricing.cpp


AM_CPPFLAGS = -g $(BOOST_CPPFLAGS) -Wall -Wextra $(GSL_CPPFLAGS) $(LUAHEADER)
AM_LDFLAGS = $(LDFLAGS) -lboost_program_options-mt $(BOOST_LDFLAGS) $(GSL_LDFLAGS) $(GSL_LIBS) $(BLAS_LIBS) 
markovchain_LDADD = markovsolverlib/libmarkovsolver.la 
markovchainlogval_LDADD = markovsolverlib/libmarkovsolver.la 

noinst_PROGRAMS = markovchain markovchainlogval

EXTRA_DIST = *.lua

markovchain_SOURCES = $(markovchain_src) \
											$(problemcalloption_src) \
											$(solverexplicitgm_src) \
											$(solverimplicitbm_src) \
											$(problemlua_src) \
											$(problemindiffpricing_src)

markovchainlogval_SOURCES = $(markovchain_src) \
											$(problemcalloption_src) \
											$(solverexplicitgm_src) \
											$(solverimplicitbm_src) \
 											$(problemlua_src) \
											$(problemindiffpricing_src)


markovchainlogval_CPPFLAGS = -DLOG_VAL


SUBDIRS = markovsolverlib 
