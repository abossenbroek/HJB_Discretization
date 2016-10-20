markovchain_src = markovchain.cpp \
					markovchain.h 
			
solver_src = solver.h \
										 solver.cpp

problem_src = problem.h

AM_CPPFLAGS = -g -Wall -Wextra $(BOOST_CPPFLAGS) $(GSL_CPPFLAGS)
AM_LDFLAGS = $(LDFLAGS) -lboost_program_options-mt  $(BOOST_LDFLAGS) $(GSL_LDFLAGS) $(GSL_LIBS) 

noinst_PROGRAMS = markovchain 

markovchain_SOURCES = $(markovchain_src) \
											$(solver_src) \
											$(problem_src)

