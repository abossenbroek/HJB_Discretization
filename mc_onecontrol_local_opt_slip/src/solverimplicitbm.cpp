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
#else
  #include <gsl/gsl_cblas.h>
  #include <gsl/gsl_blas.h>
#endif

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>

#include <lis.h>

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

#include <emmintrin.h>

#include "solverimplicitbm.h"
#include "markovsolverlib/solver.h"
#include "markovsolverlib/problem.h"

#define VEC_LAST_ELM(a) (*gsl_vector_ptr(a, a->size - 1))
#define VEC_FRST_ELM(a) (*gsl_vector_ptr(a, 0))

#define NUM_BANDS 3

SolverImplicitBM::SolverImplicitBM(Problem* _p) :
  Solver(_p)
{
  val_old = gsl_vector_alloc(p->state_size);
  val_new = gsl_vector_alloc(p->state_size);
  opt_ctl = gsl_vector_alloc(p->state_size);
}

SolverImplicitBM::~SolverImplicitBM()
{
  gsl_vector_free(val_old);
  gsl_vector_free(val_new);
  gsl_vector_free(opt_ctl);
}

void
SolverImplicitBM::solve()
{
	int time_idx;
  size_t i = 0;
	int policy_idx;

  double time;

  double control;
  double state;

  // Required to compute the boundary value.
  double prob_bnd_up;
  double prob_bnd_down;

  gsl_vector* prob_curr_state_vec = gsl_vector_calloc(p->state_size);
  gsl_vector* drift = gsl_vector_calloc(p->state_size);
  gsl_vector* diffusion = gsl_vector_calloc(p->state_size);
  gsl_vector* state_vec = gsl_vector_calloc(p->state_size);

  // Use vector views to construction the probability matrix.
  gsl_vector_view diffusion_view_sup = 
    gsl_vector_subvector(diffusion, 1, p->state_size - 1);
  gsl_vector_view diffusion_view_sub = 
    gsl_vector_subvector(diffusion, 0, p->state_size - 1);

  // Allocate vectors.  
  gsl_vector* b = gsl_vector_alloc(p->state_size);
  gsl_vector* running_cost = gsl_vector_alloc(p->state_size);
  gsl_vector* cost_to_go = gsl_vector_calloc(p->state_size);

  // Allocate the vector that will contain the sparse matrix data, column
  // information and row information.
  gsl_vector* sparse_matrix_data = 
    gsl_vector_calloc(p->state_size * 3);
  int sparse_matrix_idx[] = {-1, 0, 1};

  // Define the views that make the code compatible with previous code.
  gsl_vector_view prob_sub_diag_view = 
    gsl_vector_subvector(sparse_matrix_data, 
        1, p->state_size - 1);
  gsl_vector_view prob_sup_diag_view = 
    gsl_vector_subvector(sparse_matrix_data, 
        2 * p->state_size, p->state_size - 1);
  gsl_vector_view prob_diag_view = 
    gsl_vector_subvector(sparse_matrix_data, p->state_size, p->state_size);

  LIS_MATRIX prob_matrix;
  LIS_VECTOR b_lis;
  LIS_VECTOR cost_to_go_lis;
  LIS_SOLVER solver;
  int *lis_vec_ptrs = new int[p->state_size];
  for (i = 0; i < p->state_size; i++)
    lis_vec_ptrs[i] = i;

  lis_matrix_create(0, &prob_matrix);
  lis_matrix_set_size(prob_matrix, 0, p->state_size);
  lis_matrix_set_dia(3, sparse_matrix_idx, sparse_matrix_data->data, 
      prob_matrix);
  lis_matrix_assemble(prob_matrix);

  lis_vector_create(0, &b_lis);
  lis_vector_create(0, &cost_to_go_lis);

  lis_solver_create(&solver);
  lis_solver_set_option("-print mem", solver);
  lis_solver_set_optionC(solver);

  // Declare the three diagonal vectors.
  gsl_vector* prob_sub_diag = &prob_sub_diag_view.vector;
  gsl_vector* prob_sup_diag = &prob_sup_diag_view.vector;
  gsl_vector* prob_diag = &prob_diag_view.vector;

  double two_zeros[] = {0.0f, 0.0f};
      
  gsl_vector_view prob_curr_state_vec_sup_view =
    gsl_vector_subvector(prob_curr_state_vec, 1, prob_diag->size - 1);
  gsl_vector_view prob_curr_state_vec_sub_view =
    gsl_vector_subvector(prob_curr_state_vec, 0, prob_diag->size - 1);

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
  progress_display show_progress((p->time_size - 1) * p->policy_size);

  gsl_vector_set_all(b, 1.0);
	// Work backwards in time.
	for (time_idx = (p->time_size - 1); time_idx >= 0; time_idx--) {
		// Create a temporary variable for time.
		time = time_idx * p->time_step;

    // Initialize the value function to a minimum.
    gsl_vector_set_all(val_new, -HUGE_VAL);
    // The optimal control is not yet known.
    gsl_vector_set_all(opt_ctl, NAN);

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

      // First compute the probability that the current which should be used to
      // normalize the probabilities.
      //
      // Compute
      // prob_curr_state_vec = 1 + drift + diffusion
      // 
      // which is:
      // prob_curr_state_vec = 1 + delta / h * b() + delta / h^2 sigma()
      //
      gsl_vector_set_all(prob_curr_state_vec, 1.0);
      gsl_vector_add(prob_curr_state_vec, drift);
      gsl_vector_add(prob_curr_state_vec, diffusion);
      
      // The main diagonal should hold 1 since we have an equation of the form:
      // (I - P) * x = b
      gsl_vector_set_all(prob_diag, 1.0);

      // For the downwards and upwards probabilities the diffusion vector must
      // be scaled.
      //
      // The diffusion is now: 0.5 * delta / h^2 * sigma
      gsl_vector_scale(diffusion, 0.5);

      double *prob_sub_diag_end = prob_sub_diag->data 
        + prob_sub_diag->size - (prob_sub_diag->size % 2);
      double *prob_sub_diag_elm = prob_sub_diag->data;
      double *prob_sup_diag_elm = prob_sup_diag->data;
      double *drift_elm = drift->data;
      __m128d packed_zeros = _mm_loadu_pd(two_zeros);

      // 
      // SSE2 code used to decrease the amount of time required to compute the
      // probability matrix.
      //
      while (prob_sub_diag_elm != prob_sub_diag_end) {
        // Prefetch
        _mm_prefetch(drift_elm + 2, _MM_HINT_T0);
        _mm_prefetch(drift_elm + 3, _MM_HINT_T0);

        // Load the drift elements.
        __m128d packed_drift0 = _mm_loadu_pd(drift_elm);
        __m128d packed_drift1 = _mm_loadu_pd(drift_elm + 1);

        // 
        // First step for the downwards probability.
        //
        __m128d packed_sub_res = _mm_min_pd(packed_drift0, packed_zeros);
        // 
        // First step for the upwards probability.
        //
        __m128d packed_sup_res = _mm_max_pd(packed_drift1, packed_zeros);

        // Store the results in the array.
        _mm_storeu_pd(prob_sub_diag_elm, packed_sub_res);
        _mm_storeu_pd(prob_sup_diag_elm, packed_sup_res);

        // Increase the pointers.
        prob_sub_diag_elm += 2;
        prob_sup_diag_elm += 2;
        drift_elm += 2;
      }

      // For some cases it is the case that not the entire array is computed in
      // the loop above. This loop takes care of the rest.
      while (prob_sub_diag_elm != prob_sub_diag->data + prob_sub_diag->size) {
        // 
        // First step for the downwards probability.
        //
        (*prob_sub_diag_elm) = GSL_MIN(*drift_elm, 0.0);
        (*prob_sup_diag_elm) = GSL_MAX(*(drift_elm + 1), 0.0);
        prob_sub_diag_elm++;
        prob_sup_diag_elm++;
        drift_elm++;
      }
      
      // 
      // Compute the downwards probability.
      //
      // Add the lower elements of the diffusion vector to the probability
      // vector.
      gsl_vector_add(prob_sub_diag, &diffusion_view_sub.vector);
      // Divide by the probability to remain in the current state.
      gsl_vector_div(prob_sub_diag, &prob_curr_state_vec_sub_view.vector);
      // Compensate for the fact that the matrix is pulled to the left hand side
      // of the equation for the cost-to-go.
      gsl_vector_scale(prob_sub_diag, -1.0);

      // 
      // Compute the upwards probability.
      //
      // Add the lower elements of the diffusion vector to the probability
      // vector.
      gsl_vector_add(prob_sup_diag, &diffusion_view_sup.vector);
      // Divide by the probability to remain in the current state.
      gsl_vector_div(prob_sup_diag, &prob_curr_state_vec_sup_view.vector);
      gsl_vector_scale(prob_sup_diag, -1.0);

      // Compute the remaining probability required for the boundary.
      prob_bnd_up = (VEC_LAST_ELM(diffusion) + GSL_MAX(VEC_LAST_ELM(drift), 0.0))
        / VEC_LAST_ELM(prob_curr_state_vec);
      prob_bnd_down = (VEC_FRST_ELM(diffusion) + GSL_MIN(VEC_FRST_ELM(drift), 0.0))
        / VEC_FRST_ELM(prob_curr_state_vec);
      // Finished building the probability matrix.

      gsl_vector_memcpy(b, val_old);
      // Compute the left hand side of the equation.
      p->running_cost(running_cost, time, state_vec, control);
      // Compute an intermediate result in place to save memory.
      gsl_vector_scale(running_cost, p->time_step);

      gsl_vector_add(b, running_cost);
      gsl_vector_div(b, prob_curr_state_vec);
    
      // Boundary conditions.
      VEC_FRST_ELM(b) += prob_bnd_down 
        * p->boundary_state_min(time, VEC_FRST_ELM(val_old),
            *gsl_vector_ptr(val_old, 1));

      VEC_LAST_ELM(b) += prob_bnd_up 
        * p->boundary_state_max(time, VEC_LAST_ELM(val_old), 
            *gsl_vector_ptr(val_old, val_old->size - 2));

      // Most probably this is incredibly inefficient but we just use it as a
      // test.
      lis_vector_set_values(LIS_INS_VALUE, p->state_size, lis_vec_ptrs, 
          b->data, b_lis);
      
      lis_solve(prob_matrix, b_lis, cost_to_go_lis, solver);

      lis_vector_get_values(cost_to_go_lis, 0, p->state_size, cost_to_go->data);

      // Solve the system of linear equations.
//      gsl_linalg_solve_tridiag(prob_diag, prob_sup_diag, prob_sub_diag, b, 
//          cost_to_go);

      // Copy the cost to go the optimal value if the maximum of the cost to go
      // sub vector is larger than the optimal value found until now. 
      // The sub vectors consist of the vector minus the boundaries.
      
      double *cost_to_go_elm = cost_to_go->data;
      double *cost_to_go_end = cost_to_go->data + cost_to_go->size;
      double *val_new_elm = val_new->data;
      double *opt_ctl_elm = opt_ctl->data;

      while (cost_to_go_elm != cost_to_go_end) {
        if ((*val_new_elm = GSL_MAX(*val_new_elm, *cost_to_go_elm)) 
            == *cost_to_go_elm)
          *opt_ctl_elm = control;

        cost_to_go_elm++;
        val_new_elm++;
        opt_ctl_elm++;
      }

      // Increase the progress bar.
      ++show_progress;
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
  gsl_vector_free(b);
  gsl_vector_free(running_cost);
  gsl_vector_free(cost_to_go);
  gsl_vector_free(drift);
  gsl_vector_free(diffusion);
  gsl_vector_free(state_vec);
  gsl_vector_free(prob_curr_state_vec);
  gsl_vector_free(sparse_matrix_data);

  lis_vector_destroy(cost_to_go_lis);
  lis_vector_destroy(b_lis);
  lis_matrix_destroy(prob_matrix);
  lis_solver_destroy(solver);
}
  
/* vim:set spell spelllang=en_gb cindent foldmethod=syntax tw=80 ts=2 et: */

