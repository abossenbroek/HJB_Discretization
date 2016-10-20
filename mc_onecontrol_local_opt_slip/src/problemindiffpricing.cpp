/*
 * =====================================================================================
 *
 *       Filename:  problem_indiffpricing.cpp
 *
 *    Description:  Execute functions from a lua file to define the problem.
 *
 *        Version:  1.0
 *        Created:  22/01/2009 09:09:28
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Anton Bossenbroek (ab), anton.bossenbroek@me.com
 *        Company:  University of Amsterdam
 *
 * =====================================================================================
 */

#include <gsl/gsl_vector.h>

#include <cmath>

#include "markovsolverlib/problem.h"
#include "problemindiffpricing.h"


extern "C" {
	#include <lua.h>
	#include <lualib.h>
	#include <lauxlib.h>
}

ProblemIndiffpricing::ProblemIndiffpricing(const string filename) :
	Problem(filename)
{
	// Load additional parameters from the lua file.
	get_lua_num_param("mu", mu, PARAM_CHECK_NONE);
	get_lua_num_param("r", r, PARAM_CHECK_NONE);
	get_lua_num_param("sigma", sigma, PARAM_CHECK_NONE);
	get_lua_num_param("gamma", gamma, PARAM_CHECK_NONE);
}


void 
ProblemIndiffpricing::terminal_cost(gsl_vector* terminal_cost, 
    gsl_vector* state)
{
  for (size_t i = 0; i < state->size; i++) 
    (*gsl_vector_ptr(terminal_cost, i)) = 
      -exp(-gamma * *gsl_vector_ptr(state, i));

}

double 
ProblemIndiffpricing::boundary_state_max(
    __attribute__ ((unused)) double& time, double& last_state, 
    double& prev_state)
{
  return 2 * last_state - prev_state;
}

double 
ProblemIndiffpricing::boundary_state_min(
    __attribute__ ((unused)) double& time, double& last_state,
    double& prev_state)
{
  return 2 * last_state - prev_state;
}

void 
ProblemIndiffpricing::running_cost(gsl_vector* running_cost, 
    __attribute__ ((unused)) double& time, 
    gsl_vector* state, 
    __attribute__ ((unused)) double& control)
{
  double *running_cost_elm = running_cost->data;
  double *running_cost_elm_end = running_cost->data + running_cost->size;
  double *state_elm = state->data;

  while (running_cost_elm != running_cost_elm_end) {
    (*running_cost_elm) = -exp(-gamma * *state_elm);
    state_elm++;
    running_cost_elm++;
  }
//  for (size_t i = 0; i < state->size; i++) 
//    (*gsl_vector_ptr(running_cost, i)) = 
//      -exp(gamma * *gsl_vector_ptr(state, i));
}

void 
ProblemIndiffpricing::drift(gsl_vector* drift, 
    __attribute__ ((unused)) double& time, 
    __attribute__ ((unused)) gsl_vector* state, 
    double& control)
{
  gsl_vector_set_all(drift, r + control * (mu - r));
}

void 
ProblemIndiffpricing::diffusion(gsl_vector* diffusion,
     __attribute__ ((unused)) double& time, 
     __attribute__ ((unused)) gsl_vector* state, 
     double& control)
{
  gsl_vector_set_all(diffusion, sigma * control);
}

/* vim:set spell spelllang=en_gb cindent foldmethod=syntax tw=80 ts=2 et: */

