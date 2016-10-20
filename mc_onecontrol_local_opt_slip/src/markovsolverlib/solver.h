/*
 * =====================================================================================
 *
 *       Filename:  solver.h
 *
 *    Description:  The abstract class of a solver.
 *
 *        Version:  1.0
 *        Created:  26/01/2009 15:24:03
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Anton Bossenbroek (ab), anton.bossenbroek@me.com
 *        Company:  University of Amsterdam
 *
 * =====================================================================================
 */

#ifndef _SOLVER_H_
#define _SOLVER_H_

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

#include <gsl/gsl_vector.h>

#include "problem.h"

class Solver
{
protected:
  Problem*  p;

  // Use a static logger.
  void log_opt_val(gsl_vector* state);
  void log_opt_ctl(double& ctl);
  void log_opt_ctl(gsl_vector* ctl);

public:
  Solver(Problem* p);
  virtual ~Solver()
  {
  }

  virtual void solve() = 0;
};

#endif 

/* vim:set spell spelllang=en_gb cindent foldmethod=syntax tw=80 ts=2 et: */

