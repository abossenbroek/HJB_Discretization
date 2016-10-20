/*
 * =====================================================================================
 *
 *       Filename:  problemcalloption.cpp
 *
 *    Description:  Impelementation of a Call Option.
 *
 *        Version:  1.0
 *        Created:  22/01/2009 09:36:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Anton Bossenbroek (ab), anton.bossenbroek@me.com
 *        Company:  University of Amsterdam
 *
 * =====================================================================================
 */

#include <cmath>

#include <iostream>
using namespace std;

#include <gsl/gsl_vector.h>

#include "markovsolverlib/problem.h"
#include "problemcalloption.h"


ProblemCallOption::ProblemCallOption(const string filename) :
	Problem(filename)
{
	// Load additional parameters from the lua file.
	get_lua_num_param("volatility", volatility, PARAM_CHECK_POS);
	get_lua_num_param("interest", interest, PARAM_CHECK_POS);
	get_lua_num_param("strike", strike, PARAM_CHECK_STRICT_POS);
	get_lua_num_param("stockinit", stockinit, PARAM_CHECK_STRICT_POS);
	get_lua_num_param("maturity", maturity, PARAM_CHECK_STRICT_POS);
}

void 
ProblemCallOption::terminal_cost(gsl_vector* terminal_cost, gsl_vector* state)
{
  size_t state_idx;

  for (state_idx = 0; state_idx < state_size; state_idx++) 
      // terminal_cost = fmax(state - strike)
      gsl_vector_set(terminal_cost, state_idx, 
            fmax(*gsl_vector_ptr(state, state_idx) - strike , 0.0));
}

double 
ProblemCallOption::boundary_state_max(double& time, double& last_state,
    __attribute__ ((unused)) double& prev_state)
{
    return fmax(last_state - strike * exp(interest 
        * (maturity - time)), 0.0);
}

double 
ProblemCallOption::boundary_state_min(double& time, double& last_state,
    __attribute__ ((unused)) double& prev_state)
{
    return fmax(last_state - strike * exp(interest 
        * (maturity - time)), 0.0);
}

void 
ProblemCallOption::running_cost(
    gsl_vector* running_cost, 
    __attribute__ ((unused)) double& time, 
    __attribute__ ((unused)) gsl_vector* state, 
    __attribute__ ((unused)) double& control)
{
  gsl_vector_set_all(running_cost, 0.0);
}

void 
ProblemCallOption::drift(
    gsl_vector* drift, 
    __attribute__ ((unused)) double& time, 
    gsl_vector* state, 
    __attribute__ ((unused)) double& control)
{
  gsl_vector_memcpy(drift, state);
  gsl_vector_scale(drift, interest);
}

void 
ProblemCallOption::diffusion(gsl_vector* diffusion,
    __attribute__ ((unused)) double& time, 
    gsl_vector* state, __attribute__ ((unused)) double& control)
{
  gsl_vector_memcpy(diffusion, state);
  gsl_vector_scale(diffusion, volatility);
}

/* vim:set spell spelllang=en_gb cindent foldmethod=syntax tw=80 ts=2 et: */

