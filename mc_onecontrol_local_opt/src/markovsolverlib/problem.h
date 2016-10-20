/*
 * =====================================================================================
 *
 *       Filename:  problem.h
 *
 *    Description:  Defines the type of problem under investigation.
 *
 *        Version:  1.0
 *        Created:  14/01/2009 15:07:10
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Anton Bossenbroek (ab), anton.bossenbroek@me.com
 *        Company:  University of Amsterdam
 *
 * =====================================================================================
 */

#ifndef _PROBLEM_H_
#define _PROBLEM_H_

#include <cstring>
#include <iostream>
#include <fstream>
using namespace std;

#include <gsl/gsl_vector.h>


// Files required for lua
extern "C" {
	#include <lua.h>
	#include <lualib.h>
	#include <lauxlib.h>
}

/*
 * =====================================================================================
 *        Class:  Problem
 *  Description:  Defines the problem which has to be solved.
 * =====================================================================================
 */
class Problem
{
protected:
  
  ///
  /// Macros which permit to specify which check should be performed by
  /// get_lua_num_param.
  /// 
  enum {
    PARAM_CHECK_STRICT_POS = 0, /// Parameter should be strictly positive.
    PARAM_CHECK_NONE = 1,       /// No check has to be performed.
    PARAM_CHECK_POS = 2         /// Parameters should be positive.
  } param_val_req;

	/* the Lua interpreter */
	lua_State* L;

  ///
  /// Get a parameter from the lua state.
  /// 
  /// @param param_name   The name of the parameter.
  /// @param param_val    The reference where the value of the parameter should
  ///                     be stored.
  /// @param r            The check which should be performed by the function.
  void get_lua_num_param(const char *param_name, double& param_val, int r);

public:
  double state_step;
  double state_space_min;
  double state_space_max;
  double policy_step;
  double policy_min;
  double policy_max;
  double time_step;
  double time_max;
  size_t time_size;
  size_t state_size;
  int policy_size;
  string problem_file_name;
  ofstream *trc_file;
  ofstream *optval_file;
  ofstream *ctl_file;
  bool is_valid;
  double dh;
  double dh2;

  /// 
  /// A problem should always be configured with a lua file.
  /// 
  /// @param filename    The filename of the lua file.
  /// 
  Problem(const string filename);
  virtual ~Problem();

  /// 
  /// Compute the terminal cost for the given problem.
  ///
  /// @param terminal_cost  The vector holding the computed costs.
  /// @param state          The vector holding the states in the grid.
  ///
  virtual void terminal_cost(gsl_vector *terminal_cost, gsl_vector* state) = 0;

  /// 
  /// Give the boundary condition for the maximum of the grid.
  ///
  /// @param time    The current time.
  /// @param state   The maximum state in the grid.
  /// 
  /// @returns   The optimal value at the boundary.
  ///
  virtual double boundary_state_max(double& time, double& last_state, 
      double& prev_state) = 0;

  /// 
  /// Give the boundary condition for the minimum of the grid.
  ///
  /// @param time    The current time.
  /// @param state   The minimum state in the grid.
  /// 
  /// @returns   The optimal value at the boundary.
  ///
	virtual double boundary_state_min(double& time, double& last_state,
      double& prev_state) = 0;

  /// 
  /// Compute the running cost for the given problem.
  /// 
  /// @param running_cost   The vector holding the computed running costs.
  /// @param time           The time in the grid.
  /// @param state          The vector holding the states in the grid.
  /// @param control        The current control / policy in the problem.
  ///
	virtual void running_cost(gsl_vector *running_cost, double& time, 
      gsl_vector* state, double& control) = 0;

  /// 
  /// Compute the drift for the given problem.
  /// 
  /// @param drift          The vector holding the computed drift.
  /// @param time           The time in the grid.
  /// @param state          The vector holding the states in the grid.
  /// @param control        The current control / policy in the problem.
  ///
  virtual void drift(gsl_vector *drift, double& time, gsl_vector* state, 
      double& control) = 0;
  /// 
  /// Compute the diffusion for the given problem.
  /// 
  /// @param diffusion      The vector holding the computed drift.
  /// @param time           The time in the grid.
  /// @param state          The vector holding the states in the grid.
  /// @param control        The current control / policy in the problem.
  ///
	virtual void diffusion(gsl_vector *diffusion, double& time, gsl_vector* state,
      double& control) = 0;
	
  ///
  /// Print the characteristics of the problem.
  /// 
  /// @param out            The out stream to which should be printed.
  /// @param s              The problem under consideration.
  /// 
  /// @returns              the out stream plus the problem description.
  friend ostream& operator<<(ostream& out, const Problem& s);
};

#endif

/* vim:set spell spelllang=en_gb cindent foldmethod=syntax tw=80 ts=2 et: */

