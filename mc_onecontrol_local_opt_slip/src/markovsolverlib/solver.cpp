/*
 * =====================================================================================
 *
 *       Filename:  solver.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  26/01/2009 15:34:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Anton Bossenbroek (ab), anton.bossenbroek@me.com
 *        Company:  University of Amsterdam
 *
 * =====================================================================================
 */

#include "solver.h"
#include "problem.h"

#include <gsl/gsl_vector.h>

Solver::Solver(Problem* _p) :
	p(_p)
{
}

void 
Solver::log_opt_val(gsl_vector* state) 
{
  size_t state_idx;

  // Write the old state to file.
  for (state_idx = 0; state_idx < p->state_size; state_idx++) 
    (*(p->trc_file)) << *gsl_vector_ptr(state, state_idx) << " ";

  // Add a end of line at the end of the line ;)
  (*(p->trc_file)) << endl;
}

void
Solver::log_opt_ctl(double& ctl)
{
  (*(p->ctl_file)) << ctl << endl;
}

void
Solver::log_opt_ctl(gsl_vector* ctl)
{
  double *ctl_elm = ctl->data;
  double *ctl_end = ctl->data + ctl->size;

  while (ctl_elm != ctl_end) {
    (*(p->ctl_file)) << *ctl_elm << " ";
    ctl_elm++;
  }
  (*(p->ctl_file)) << endl;
}

/* vim:set spell spelllang=en_gb cindent foldmethod=syntax tw=80 ts=2 et: */

