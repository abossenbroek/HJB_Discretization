/*
 * =====================================================================================
 *
 *       Filename:  problem.cpp
 *
 *    Description:  Reads a problem from a file where the problem is defined by a lua
 *									script.
 *
 *        Version:  1.0
 *        Created:  14/01/2009 18:11:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Anton Bossenbroek (ab), anton.bossenbroek@me.com
 *        Company:  University of Amsterdam
 *
 * =====================================================================================
 */

#include <string>
#include <fstream>
using namespace std;

#include <config.h>

#include <cmath>

// Files required for lua
extern "C" {
	//#include <luaxlib.h>
	#include <lua.h>
	#include <lualib.h>
}


#include "problem.h"


Problem::Problem(const string filename) :
	problem_file_name(filename)
{

	// initialize Lua 
	L = luaL_newstate();

	// load various Lua libraries 
	luaL_openlibs(L);

	luaL_dofile(L, filename.c_str());

	get_lua_num_param("statestep", state_step, PARAM_CHECK_STRICT_POS);
	get_lua_num_param("statemin", state_space_min, PARAM_CHECK_NONE);
	get_lua_num_param("statemax", state_space_max, PARAM_CHECK_NONE);
	get_lua_num_param("policystep", policy_step, PARAM_CHECK_STRICT_POS);
	get_lua_num_param("policymin", policy_min, PARAM_CHECK_NONE);
	get_lua_num_param("policymax", policy_max, PARAM_CHECK_NONE);
	get_lua_num_param("timestep", time_step, PARAM_CHECK_STRICT_POS);
	get_lua_num_param("timemax", time_max, PARAM_CHECK_STRICT_POS);

  if (state_space_min > state_space_max) {
    cerr << "State minimum should be smaller than maximum" << endl;
    throw new exception();
  }
  if (policy_min > policy_max) {
    cerr << "Policy minimum should be smaller than maximum" << endl;
    throw new exception();
  }

	// Read the trace file parameter from lua.	
	lua_getglobal(L, "tracefile");
	if (!lua_isstring(L, -1)) {
		cerr << "tracefile should be a string" << endl;
		throw new exception();
	}
	string tracefn = lua_tostring(L, -1);
	// Open the trace file.
	trc_file = new ofstream(tracefn.c_str());

	// Read the optimal value file parameter from lua.	
	lua_getglobal(L, "optimalvaluefile");
	if (!lua_isstring(L, -1)) {
		cerr << "optimalvaluefile should be a string" << endl;
		throw new exception();
	}
	string optvalfn = lua_tostring(L, -1);
	// Open the optimal value file.
	optval_file = new ofstream(optvalfn.c_str());

  // Read the optimal value file parameter from lua.	
	lua_getglobal(L, "controlfile");
	if (!lua_isstring(L, -1)) {
		cerr << "controlfile should be a string" << endl;
		throw new exception();
	}
	string ctlfn = lua_tostring(L, -1);
	// Open the optimal value file.
	ctl_file = new ofstream(ctlfn.c_str());


	// Compute the size of the grids to use.
	time_size = static_cast<int>(lround(time_max / time_step));
	cerr << "Time size " << time_size << endl;
	state_size = static_cast<int>(lround((state_space_max - state_space_min) /
			state_step));
	cerr << "State size " << state_size << endl;
	policy_size = static_cast<int>(lround((policy_max - policy_min) /
			policy_step));
	cerr << "policy size " << policy_size << endl;
	
	// Compute these values because they will be used a lot.
	this->dh = time_step / state_step;
	this->dh2 = time_step / (state_step * state_step);
	cerr << "dh: " << dh << " dh2: " << dh2 << endl;
}

Problem::~Problem()
{
	trc_file->close();
	optval_file->close();
  ctl_file->close();
	delete trc_file;
	delete optval_file;
	/* cleanup Lua */
	lua_close(L);
}

ostream& operator<<(ostream& out, const Problem& p)
{
	return out << "Policy [" << p.policy_min << ", " << p.policy_max 
		<< "] step " << p.policy_step << endl
		<< "Time [0, " << p.time_max << "] step " << p.time_step << endl
		<< "State [" << p.state_space_min << ", " << p.state_space_max 
		<< "] step " << p.state_step << endl;

	return out;
}

void
Problem::get_lua_num_param(const char *param_name, double& param_val, int r)
{
	lua_getglobal(L, param_name);

	if (!lua_isnumber(L, -1)) {
		cerr << param_name << " must be a number" << endl;
		throw new exception();
	}

	param_val = static_cast<double>(lua_tonumber(L, -1));

	switch (r) {
		case PARAM_CHECK_STRICT_POS:
			if (!(param_val > 0)) {
				cerr << param_name << " is " << param_val << endl;
				cerr << param_name << " must by strictly positive." << endl;
				throw new exception();
			}
			break;

		case PARAM_CHECK_POS:
			if (param_val < 0) {
				cerr << param_name << " is " << param_val << endl;
				cerr << param_name << " must by strictly positive." << endl;
				throw new exception();
			}
			break;

	//	case PARAM_CHECK_NONE:
	//	default:
	}
}

/* vim:set spell spelllang=en_gb cindent foldmethod=syntax tw=80 ts=2 et: */

