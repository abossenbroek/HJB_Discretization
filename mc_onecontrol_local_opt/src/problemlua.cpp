/*
 * =====================================================================================
 *
 *       Filename:  problemlua.cpp
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

#include "markovsolverlib/problem.h"
#include "problemlua.h"


extern "C" {
	#include <lua.h>
	#include <lualib.h>
	#include <lauxlib.h>
}

void 
ProblemLua::terminal_cost(gsl_vector* terminal_cost, 
    gsl_vector* state)
{
  double *terminal_cost_elm = terminal_cost->data;
  double *terminal_cost_end = terminal_cost->data + terminal_cost->size;
  double *state_elm = state->data;

  while (terminal_cost_elm != terminal_cost_end) {
    lua_getglobal(L, "terminal_cost");
    lua_pushnumber(L, *state_elm);
	  lua_call(L, 1, 1);
    (*terminal_cost_elm) = static_cast<double>(lua_tonumber(L, -1));
    lua_pop(L, 1);
    terminal_cost_elm++;
    state_elm++;
  }
}

double 
ProblemLua::boundary_state_max(double& time, double& last_state, 
    double& prev_state)
{
	double lua_val;
	
	lua_getglobal(L, "boundary_state_max");
	lua_pushnumber(L, time);
	lua_pushnumber(L, last_state);
	lua_pushnumber(L, prev_state);
	lua_call(L, 3, 1);
	
	lua_val = (double)lua_tonumber(L, -1);
	lua_pop(L, 1);

	return lua_val;
}

double 
ProblemLua::boundary_state_min(double& time, double& last_state,
    double& prev_state)
{
	double lua_val;
	
	lua_getglobal(L, "boundary_state_min");
	lua_pushnumber(L, time);
	lua_pushnumber(L, last_state);
	lua_pushnumber(L, prev_state);
	lua_call(L, 3, 1);
	
	lua_val = (double)lua_tonumber(L, -1);
	lua_pop(L, 1);

	return lua_val;
}

void 
ProblemLua::running_cost(gsl_vector* running_cost, double& time, 
     gsl_vector* state, double& control)
{
  for (size_t i = 0; i < state->size; i++) {
    lua_getglobal(L, "running_cost");
    lua_pushnumber(L, time);
    lua_pushnumber(L, control);
    lua_pushnumber(L, *gsl_vector_ptr(state, i));
	  lua_call(L, 3, 1);
    (*gsl_vector_ptr(running_cost, i)) 
      = static_cast<double>(lua_tonumber(L, -1));
    lua_pop(L, 1);
  }
}

void 
ProblemLua::drift(gsl_vector* drift, double& time, 
    gsl_vector* state, double& control)
{
  for (size_t i = 0; i < state->size; i++) {
    lua_getglobal(L, "drift");
    lua_pushnumber(L, time);
    lua_pushnumber(L, control);
    lua_pushnumber(L, *gsl_vector_ptr(state, i));
	  lua_call(L, 3, 1);
    (*gsl_vector_ptr(drift, i)) 
      = static_cast<double>(lua_tonumber(L, -1));
    lua_pop(L, 1);
  }
}

void 
ProblemLua::diffusion(gsl_vector* diffusion,
     double& time, __attribute__ ((unused)) gsl_vector* state, 
     double& control)
{
  for (size_t i = 0; i < state->size; i++) {
    lua_getglobal(L, "diffusion");
    lua_pushnumber(L, time);
    lua_pushnumber(L, control);
    lua_pushnumber(L, *gsl_vector_ptr(state, i));
	  lua_call(L, 3, 1);
    (*gsl_vector_ptr(diffusion, i)) 
      = static_cast<double>(lua_tonumber(L, -1));
    lua_pop(L, 1);
  }
}

/* vim:set spell spelllang=en_gb cindent foldmethod=syntax tw=80 ts=2 et: */

