/*
 * =====================================================================================
 *
 *       Filename:  solver.cpp
 *
 *    Description:  The actual solver of the utility maximization defined in a
 *                  problem.
 *
 *        Version:  1.0
 *        Created:  08/01/2009 13:22:16
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Anton Bossenbroek (ab), anton.bossenbroek@me.com
 *        Company:  University of Amsterdam
 *
 * =====================================================================================
 */

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

#ifdef HAVE_ACCELERATE
  #include <Accelerate/Accelerate.h>
#elif defined(HAVE_BLAS)
  #include <cblas.h>
#else
  #include <gsl/gsl_blas.h>
  #include <gsl/gsl_cblas.h>
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
using namespace std;

#include <boost/progress.hpp>
using namespace boost;

#include "solverexplicitgm.h"
#include "markovsolverlib/solver.h"
#include "markovsolverlib/problem.h"

#define VEC_LAST_ELM(a) (*gsl_vector_ptr(a, a->size - 1))
#define VEC_FRST_ELM(a) (*gsl_vector_ptr(a, 0))

#define NUM_BANDS 3

// Defines which percentage of the cost-to-go should be used to evaluate the
// maximum. It counts from the centre.
#define WINDOW_SIZE 0.7

// Change this with care!!!
#define MINMAX 1

SolverExplicitGM::SolverExplicitGM(Problem* _p) :
  Solver(_p)
{
  val_old = gsl_vector_alloc(p->state_size);
  val_new = gsl_vector_alloc(p->state_size);
}

SolverExplicitGM::~SolverExplicitGM()
{
  gsl_vector_free(val_old);
  gsl_vector_free(val_new);
}

void
SolverExplicitGM::solve()
{
	int time_idx;
  size_t i = 0;
	int policy_idx;

  double time;

  double control;
  double opt_ctl;
  double state;

  // Required to compute the boundary value.
  double prob_bnd_up;
  double prob_bnd_down;

  gsl_vector* prob_diag_vec = gsl_vector_calloc(p->state_size);
  gsl_vector* prob_su_diag_vec = gsl_vector_calloc(p->state_size - 1);
  gsl_vector* drift = gsl_vector_calloc(p->state_size);
  gsl_vector* diffusion = gsl_vector_calloc(p->state_size);
  gsl_vector* state_vec = gsl_vector_calloc(p->state_size);

  // Use vector views to construction the probability matrix.
  gsl_vector_const_view diffusion_view_sup = 
    gsl_vector_const_subvector(diffusion, 0, p->state_size - 1);
  gsl_vector_const_view diffusion_view_sub = 
    gsl_vector_const_subvector(diffusion, 1, p->state_size - 1);

  // Allocate vectors.  
  gsl_vector* cost_to_go = gsl_vector_alloc(p->state_size);

  // Use a general matrix to store the probabilities.
  gsl_matrix* prob_matrix = gsl_matrix_calloc(p->state_size, p->state_size);

  // Declare the three diagonal views.
  gsl_vector_view prob_sub_diag_view = 
    gsl_matrix_subdiagonal(prob_matrix, 1);
  gsl_vector_view prob_sup_diag_view = 
    gsl_matrix_superdiagonal(prob_matrix, 1);
  gsl_vector_view prob_diag_view = 
    gsl_matrix_diagonal(prob_matrix);

  int subvector_minimum = lround((cost_to_go->size - WINDOW_SIZE 
        * (cost_to_go->size - 1)) / 2.0);
  int subvector_size = lround(WINDOW_SIZE * (cost_to_go->size - 1));

  cerr << "Subvector minimum " << subvector_minimum << endl;
  cerr << "Subvector size " << subvector_size << endl;
      
  gsl_vector_view cost_to_go_subvector = 
    gsl_vector_subvector(cost_to_go, subvector_minimum, subvector_size);
  gsl_vector_view val_new_subvector = 
    gsl_vector_subvector(val_new, subvector_minimum, subvector_size);
  // Initialize the state vector and save the values to file.
  for (i = 0; i < p->state_size; i++) {
    // Create a temporary variable.
    state = p->state_step * i + p->state_space_min;
    // Save the scalar to the vector.
    gsl_vector_set(state_vec, i, state);
  }

  // Initialize the old value function with the terminal cost.
  p->terminal_cost(val_old, state_vec);
	
  // Print states to file.
  for (i = 0; i < state_vec->size; i++) 
    (*(p->optval_file)) << gsl_vector_get(state_vec, i) << " ";
  // Add new line to file.
  (*(p->optval_file)) << endl;
  // Save the terminal cost to file.
	for (i = 0; i < val_old->size; i++) 
    (*(p->optval_file)) << gsl_vector_get(val_old, i) << " ";
  // Add new line to file.
  (*(p->optval_file)) << endl;
    
  // Initialize a progress bar.
#ifndef DEBUG_PRINT
  progress_display show_progress((p->time_size - 1) * p->policy_size);
#endif // DEBUG_PRINT

#ifdef DEBUG_PRINT
  cerr << setiosflags(ios::fixed) << setprecision(4);
#endif // DEBUG_PRINT

#ifdef DEBUG_PRINT
  cerr << "**** DEBUG the state vector: " << endl;
  for (i = 0; i < state_vec->size; i++) 
    cerr << *gsl_vector_ptr(state_vec, i) << endl;
#endif // DEBUG_PRINT

	// Work backwards in time.
	for (time_idx = (p->time_size - 2); time_idx >= 0; time_idx--) {
		// Create a temporary variable for time.
		time = time_idx * p->time_step;

    // Initialize the value function to a minimum.
    gsl_vector_set_all(val_new, -HUGE_VAL);

    // The optimal control is not yet known.
    opt_ctl = NAN;
 
    // Iterate over policy space.
    for (policy_idx = 0; policy_idx < p->policy_size; policy_idx++) {
      // Determine the current value of the control.
      control = p->policy_min + policy_idx * p->policy_step;

      //
      // Create the probability matrix.
      //
      // Compute the drift and diffusion for this control and multiply by the
      // appropriate variables.
      p->diffusion(diffusion, time, state_vec, control);
      // diffusion = diffusion^2
      gsl_vector_mul(diffusion, diffusion);
      p->drift(drift, time, state_vec, control);
      // diffusion =  p->dh2 * diffusion;
      gsl_vector_scale(diffusion, p->dh2);
      // drift = p->dh * drift
      gsl_vector_scale(drift, p->dh);

      // First compute the straight probabilities.
      // Initialize and perform operation on the prob_diag_vec.
      for (i = 0; i < prob_diag_vec->size; i++) 
          // prob_diag_vec[i] = 1 - fabs(drift[i])
          gsl_vector_set(prob_diag_vec, i, 1.0 
              - fabs(*gsl_vector_ptr(drift, i)));
      
      // prob_diag_vec -= diffusion
      gsl_vector_sub(prob_diag_vec, diffusion);
      // Copy the vector to the matrix.
      gsl_vector_memcpy(&prob_diag_view.vector, prob_diag_vec);

      // For the downwards and upwards probabilities the diffusion vector must
      // be scaled.
      gsl_vector_scale(diffusion, 0.5);

#ifdef DEBUG_PRINT
      cerr << "*** DEBUG diffusion vector: ";
      for (i = 0; i < diffusion->size; i++)
        cerr << *gsl_vector_ptr(diffusion, i) << " ";
      cerr << endl;
      cerr << "*** DEBUG drift vector    : ";
      for (i = 0; i < drift->size; i++)
        cerr << *gsl_vector_ptr(drift, i) << " ";
      cerr << endl;
#endif // DEBUG_PRINT

      // Compute the downwards probability.
      // Initialize and perform operation on the prob_diag_vec.
#ifdef MINMAX
      for (i = 0; i < prob_su_diag_vec->size; i++) 
          // prob_diag_vec[i] = GSL_MIN(drift[i], 0.0)
          gsl_vector_set(prob_su_diag_vec, i, 
            GSL_MIN(*gsl_vector_ptr(drift, i + 1), 0.0));
#else
      for (i = 0; i < prob_su_diag_vec->size; i++) 
          // prob_diag_vec[i] = GSL_MIN(drift[i], 0.0)
          gsl_vector_set(prob_su_diag_vec, i, 
            GSL_MAX(*gsl_vector_ptr(drift, i + 1), 0.0));
#endif  // MINMAX
      // Add the lower elements of the diffusion vector to the probability
      // vector.
      gsl_vector_add(prob_su_diag_vec, &diffusion_view_sub.vector);
      // Copy the vector to the matrix.
      gsl_vector_memcpy(&prob_sub_diag_view.vector, prob_su_diag_vec);
      
      // Compute the upwards probability.
      // Initialize and perform operation on the prob_diag_vec.
#ifdef MINMAX
      for (i = 0; i < prob_su_diag_vec->size; i++) 
          // prob_diag_vec[i] = GSL_MAX(drift[i], 0.0)
          gsl_vector_set(prob_su_diag_vec, i, 
            GSL_MAX(*gsl_vector_ptr(drift, i), 0.0));
#else
      for (i = 0; i < prob_su_diag_vec->size; i++) 
          // prob_diag_vec[i] = GSL_MAX(drift[i], 0.0)
          gsl_vector_set(prob_su_diag_vec, i, 
            GSL_MIN(*gsl_vector_ptr(drift, i), 0.0));
#endif  // MINMAX
      // Add the lower elements of the diffusion vector to the probability
      // vector.
      gsl_vector_add(prob_su_diag_vec, &diffusion_view_sup.vector);
      // Copy the vector to the matrix.
      gsl_vector_memcpy(&prob_sup_diag_view.vector, prob_su_diag_vec);

#ifdef MINMAX
      // Compute the remaining probability required for the boundary.
      prob_bnd_down = VEC_FRST_ELM(diffusion) + GSL_MIN(VEC_FRST_ELM(drift), 0.0);
      // Compute the remaining probability.
      prob_bnd_up = VEC_LAST_ELM(diffusion) + GSL_MAX(VEC_LAST_ELM(drift), 0.0);
      // Finished building the probability matrix.
#else 
      // Compute the remaining probability required for the boundary.
      prob_bnd_down = VEC_FRST_ELM(diffusion) + GSL_MIN(VEC_FRST_ELM(drift), 0.0);
      // Compute the remaining probability.
      prob_bnd_up = VEC_LAST_ELM(diffusion) + GSL_MAX(VEC_LAST_ELM(drift), 0.0);
      // Finished building the probability matrix.
#endif  // MINMAX
#ifdef DEBUG_PRINT
      cerr << "**** DEBUG the probability matrix: " << endl;
      for (i = 0; i < prob_matrix->size1; i++) {
        for (size_t j = 0; j < prob_matrix->size2; j++) {
          cerr << *gsl_matrix_ptr(prob_matrix, i, j) << " ";
        }
        cerr << endl;
      }
#endif // DEBUG_PRINT

      // Because of BLAS the running_cost has to be stored in the cost_to_go.
      p->running_cost(cost_to_go, time, state_vec, control);

#ifdef DEBUG_PRINT
      cerr << "**** DEBUG running_cost: ";
      for (i = 0; i < cost_to_go->size; i++) 
          cerr << *gsl_vector_ptr(cost_to_go, i) << " ";
      cerr << endl;
      cerr << "**** DEBUG old value   : ";
      for (i = 0; i < val_old->size; i++) 
          cerr << *gsl_vector_ptr(val_old, i) << " ";
      cerr << endl;
      cerr << "**** DEBUG scaled diffusion: ";
      for (i = 0; i < diffusion->size; i++) 
          cerr << *gsl_vector_ptr(diffusion, i) << " ";
      cerr << endl;
      cerr << "**** DEBUG scaled drift : ";
      for (i = 0; i < drift->size; i++) 
          cerr << *gsl_vector_ptr(drift, i) << " ";
      cerr << endl;
#endif // DEBUG_PRINT
      // Compute cost_to_go = prop_matrix * val_old + time_step * cost_to_go.
      cblas_dgemv(CblasRowMajor, CblasNoTrans, prob_matrix->size1, 
          prob_matrix->size2, 1.0, prob_matrix->data, prob_matrix->tda, 
          val_old->data, val_old->stride, p->time_step, cost_to_go->data, 
          cost_to_go->stride);
      
      //gsl_blas_dgemv(CblasTrans, 1.0, prob_matrix, val_old, p->time_step, cost_to_go);
#ifdef DEBUG_PRINT 
      cerr << "*** DEBUG prob boundary down: " << prob_bnd_down << " * "
        << p->boundary_state_min(time, 
            VEC_FRST_ELM(val_old), 
            *gsl_vector_ptr(val_old, val_old->size - 2)) << " + " 
        << VEC_FRST_ELM(cost_to_go) << endl;
      cerr << "*** DEBUG prob boundary up: " << prob_bnd_up << " * "
        << p->boundary_state_max(time, 
            VEC_LAST_ELM(val_old), *gsl_vector_ptr(val_old, 1)) << " + " 
        << VEC_LAST_ELM(cost_to_go) << endl;
#endif // DEBUG_PRINT
      // Boundary conditions.
      (VEC_LAST_ELM(cost_to_go)) += p->boundary_state_max(time, 
          VEC_LAST_ELM(val_old), 
          *gsl_vector_ptr(val_old, val_old->size - 2)) * prob_bnd_up;
      
      (VEC_FRST_ELM(cost_to_go)) += p->boundary_state_min(time, 
          VEC_FRST_ELM(val_old), 
          *gsl_vector_ptr(val_old, 1)) * prob_bnd_down;
     
      // Copy the cost to go the optimal value if the maximum of the cost to go
      // sub vector is larger than the optimal value found until now. 
      // The sub vectors consist of the vector minus the boundaries.
      if (gsl_vector_max(&cost_to_go_subvector.vector) 
          > gsl_vector_max(&val_new_subvector.vector)) {
        opt_ctl = control;
        gsl_vector_memcpy(val_new, cost_to_go);
      }
#ifdef DEBUG_PRINT
      cerr << "**** DEBUG the cost_to_go value for control " << control << " : " << endl;
      for (i = 0; i < cost_to_go->size; i++) 
          cerr << *gsl_vector_ptr(cost_to_go, i) << " ";
      cerr << endl;
      cerr << "**** DEBUG winning control @ time " << time << " : " << opt_ctl << endl;
#endif // DEBUG_PRINT
      // Increase the progress bar.
#ifndef DEBUG_PRINT
      ++show_progress;
#endif // DEBUG_PRINT
    }
    gsl_vector_memcpy(val_old, val_new);
#ifdef LOG_VAL
    log_opt_val(val_new);
#endif // LOG_VAL
    // Log the optimal control.
    log_opt_ctl(opt_ctl);
  }

  for (i = 0; i < val_new->size; i++) 
    (*(p->optval_file)) << gsl_vector_get(val_new, i) << " ";

  // Add a end of line at the end of the line ;)
  (*(p->optval_file)) << endl;

  // Free all the data structures used in the solver.
  gsl_vector_free(cost_to_go);
  gsl_matrix_free(prob_matrix);
  gsl_vector_free(drift);
  gsl_vector_free(diffusion);
  gsl_vector_free(state_vec);
  gsl_vector_free(prob_diag_vec);
  gsl_vector_free(prob_su_diag_vec);
}
  
/* vim:set spell spelllang=en_gb cindent foldmethod=syntax tw=80 ts=2 et: */

