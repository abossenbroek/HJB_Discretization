abstract_classes_src = solver.h \
									 solver.cpp \
									 problem.h \
									 problem.cpp

noinst_LTLIBRARIES = libmarkovsolver.la

libmarkovsolver_la_SOURCES = $(abstract_classes_src) 
libmarkovsolver_la_CPPFLAGS = -g $(BOOST_CPPFLAGS) 
libmarkovsolver_la_LDFLAGS = -lboost_program_options-mt $(BOOST_LDFLAGS) 

