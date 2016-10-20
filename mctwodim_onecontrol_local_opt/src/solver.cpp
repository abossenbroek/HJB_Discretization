

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

// Undefine if the progress bar is not required.
//#undef HAVE_BOOST
//
#include <gsl/gsl_cdf.h>

#include <lis.h>

#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
using namespace std;

#ifdef HAVE_BOOST
# include <boost/progress.hpp>
using namespace boost;
#endif // HAVE_BOOST

#include "solver.h"
#include "problem.h"

#define NUM_BANDS 9

#define COORDINATES (STEPS * STEPS)

#define UPPER_BOUNDARY_FACTOR 2
#define LEFT_BOUNDARY_FACTOR(time)  (1.2 * (time + DELTA) + 0.8)
#define RIGHT_BOUNDARY_FACTOR(time) (1.2 * (time + DELTA) + 0.8)

// Utility function:  U(blc) = blc^\gamma
// Running cost:      U(y(t)) - u_2(t)
#define RUNNING_COST(x, y)  \
  (pow(fmax((1.0 - DRIFT_BLC) * (sqrt(DELTA) * x), 0.0), GAMMA) - EFFORT_COST * y * y)

void
Solver::solve()
{
  int i, j, idx; 

  double time;

  double control1;
  double control2;

  cerr << setiosflags(ios::fixed) << setprecision(6);

  LIS_SCALAR* prob_matrix_data = new LIS_SCALAR[COORDINATES * NUM_BANDS];

  double dh = DELTA / H;
  double dh2 = DELTA / (H * H);
  double *opt_control1 = new double[COORDINATES];
  double *opt_control2 = new double[COORDINATES];
  double C1minus;

  LIS_VECTOR cost_to_go;
  LIS_VECTOR b;
  LIS_VECTOR val_new;
  LIS_VECTOR val_old;
  LIS_MATRIX prob_matrix;
  LIS_SOLVER solver;

  lis_vector_create(0, &cost_to_go);
  lis_vector_create(0, &b);
  lis_vector_create(0, &val_new);
  lis_vector_create(0, &val_old);

  lis_vector_set_size(cost_to_go, 0, COORDINATES);
  lis_vector_set_size(b, 0, COORDINATES);
  lis_vector_set_size(val_new, 0, COORDINATES);
  lis_vector_set_size(val_old, 0, COORDINATES);


  int sparse_matrix_diags[] = {
    // The first block.
    -STEPS - 1 , -STEPS, -STEPS + 1, 
    // The second block.
    -1, 0, 1, 
    // The third block.
    STEPS - 1, STEPS, STEPS + 1};
  
  lis_matrix_create(0, &prob_matrix);
  lis_matrix_set_size(prob_matrix, 0, COORDINATES);
  lis_matrix_set_dia(NUM_BANDS, sparse_matrix_diags, prob_matrix_data, 
      prob_matrix);
  lis_matrix_assemble(prob_matrix);

  // Initialize the solver.
  lis_solver_create(&solver);
  lis_solver_set_option((char*)"-p ilut", solver);
  //lis_solver_set_optionC(solver);

  for (idx = 0; idx < NUM_BANDS * COORDINATES; idx++)
    prob_matrix_data[idx] = 0;
 
  // Compute the terminal cost.
  for (idx = 0; idx < COORDINATES; idx++) {
    // Retrieve the coordinate of the current state.
    i = idx % STEPS;
    j = idx / STEPS;
    // Compute the state.
    double blc_elm = H * j + BLC_MIN;
    double terminal_cost = RUNNING_COST(blc_elm, 0);

    // Compute the terminal cost.
    lis_vector_set_value(LIS_INS_VALUE, idx, terminal_cost, val_old);
  }
  // Dump optimal controls to file.
  for (idx = 0; idx < COORDINATES; idx++) {
    double val_old_elm;
    lis_vector_get_value(val_old, idx, &val_old_elm);
    *opt_val_end_file << val_old_elm << " ";
  }
  *opt_val_end_file << endl;
 
#ifdef HAVE_BOOST
  // Initialize a progress bar.
  progress_display show_progress(TIME_STEPS);
#endif // HAVE_BOOST

	// Work backwards in time.
	for (time = (1 - DELTA); time >= 0; time -= DELTA) {

    // Initialize the value function to a minimum.
    lis_vector_set_all(-HUGE_VAL, val_new);
    // The optimal control1 is not yet known.
    for (idx = 0; idx < COORDINATES; idx++) {
      opt_control1[idx] = NAN;
      opt_control2[idx] = NAN;
    }

    // Iterate over policy space.
    for (control1 = CONTROL1_MIN; control1 < CONTROL1_MAX; 
        control1 += (CONTROL1_MAX - CONTROL1_MIN) 
        / static_cast<double>(CONTROL1_STEPS)) {

      for (control2 = CONTROL2_MIN; control2 < CONTROL2_MAX; 
          control2 += (CONTROL2_MAX - CONTROL2_MIN) 
          / static_cast<double>(CONTROL2_STEPS)) {
        cerr << "About to compute @ time " << time << " control1 " 
          << control1 << " control2 " << control2 << endl;
        //======================================================================
        //
        // NOTE: The code has been specifically written for the bonus malus
        // problem. The probability matrix is built with 9 diagonals.
        //
        //======================================================================

        // Perform the first step to construct the vector b.
        // Copy  b <- val_old
        lis_vector_copy(val_old, b);

        for (idx = 0; idx < COORDINATES; idx++) {
          i = idx % STEPS;
          j = idx / STEPS;
          double pnl_elm = H * i + PNL_MIN;
          double blc_elm = H * j + BLC_MIN;
          // Step up and down of the current balance value. This is required for
          // the boundary conditions.
          double blc_elm_up = H * (j + 1) + BLC_MIN;
          double blc_elm_down = H * (j - 1) + BLC_MIN;

          double pnl_b = (DRIFT_PNL + IMPACT * control1) * pnl_elm;
          double blc_b = fmax(0.0, DRIFT_BLC * blc_elm) + PAYOUT * pnl_b;

          double pnl_b_p = fmax(0.0, pnl_b);
          double pnl_b_m = -fmin(0.0, pnl_b);
          double blc_b_p = fmax(0.0, blc_b);
          double blc_b_m = -fmin(0.0, blc_b);

          double axx = (DIFFUSION_PNL + IMPACT * control2) 
            * (DIFFUSION_PNL + IMPACT * control2);
          double axy = PAYOUT * PAYOUT * axx;
          double axyp = fmax(axy, 0.0);
          double axym = fmax(-axy, 0.0);
          double axya = fabs(axy);
          double ayy = PAYOUT * PAYOUT * axy;
          C1minus = -1.0 / (1.0 + dh * fabs(pnl_b) + dh * fabs(blc_b) 
              + dh2 * axx + dh2 * ayy - dh2 * axya);

          // Determine the relative index of the current coordinate.
          int diag_idx = 4 * COORDINATES + idx;
          // The diagonal will only hold 1's.
          prob_matrix_data[diag_idx] = 1;

          double running_cost = RUNNING_COST(blc_elm, control1);
          // First sub diagonal is the probability to jump to (-1, -1).
          if (((idx % STEPS) > 0) && ((idx - STEPS) >= 0)) {
            prob_matrix_data[idx] =  
              C1minus * (0.5 * dh2 * axyp);
          } else {
            double bnd_utility;
            bnd_utility = RUNNING_COST(blc_elm_down, control1);
            bnd_utility *= 0.5 * dh2 * axyp;
            lis_vector_set_value(LIS_ADD_VALUE, idx, bnd_utility, b);
          }

          // Last super diagonal is the probability to jump to (1, 1).
          if (((idx % STEPS) < (STEPS - 1)) && ((idx + STEPS) < COORDINATES)) {
            prob_matrix_data[idx + 8 * COORDINATES] = 
              C1minus * (0.5 * dh2 * axyp);
          } else {
            double bnd_utility = UPPER_BOUNDARY_FACTOR * RUNNING_COST(blc_elm_up, control1);
            bnd_utility *= 0.5 * dh2 * axyp;
            lis_vector_set_value(LIS_ADD_VALUE, idx, bnd_utility, b);
          }


          // Second sub diagonal is the probability to jump to (0, -1).
          if ((idx - STEPS) >= 0) {
            prob_matrix_data[idx + COORDINATES] = 
              C1minus * (dh * blc_b_m + 0.5 * dh2 * ayy - 0.5 * dh2 * axya);
          } else {
            double bnd_utility = RUNNING_COST(blc_elm_down, control1);
            bnd_utility *= dh * blc_b_m + 0.5 * dh2 * ayy - 0.5 * dh2 * axya;
            lis_vector_set_value(LIS_ADD_VALUE, idx, bnd_utility, b);
          }

          // Second super diagonal is the probability to jump to (0, 1).
          if ((idx + STEPS) < COORDINATES) {
            prob_matrix_data[idx + 7 * COORDINATES] = 
              C1minus * (dh * blc_b_p + 0.5 * dh2 * ayy - 0.5 * dh2 * axya);
          } else {
            double bnd_utility = UPPER_BOUNDARY_FACTOR * RUNNING_COST(blc_elm_up, control1);
            bnd_utility *= dh * blc_b_p + 0.5 * dh2 * ayy - 0.5 * dh2 * axya;
            lis_vector_set_value(LIS_ADD_VALUE, idx, bnd_utility, b);
          }


          // Third sub diagonal is the probability to jump to (1, -1)
          if (((idx % STEPS) < (STEPS - 1)) && ((idx - STEPS) >= 0)) {
            prob_matrix_data[idx + 2 * COORDINATES] =  
              C1minus * (0.5 * dh2 * axym);
          } else {
            double bnd_utility = RUNNING_COST(blc_elm_down, control1);
            bnd_utility *= 0.5 * dh2 * axym;
            lis_vector_set_value(LIS_ADD_VALUE, idx, bnd_utility, b);
          }


          // Third super diagonal is the probability to jump to (-1, 1)
          if (((idx % STEPS) > 0) && ((idx + STEPS) < COORDINATES)) {
            prob_matrix_data[idx + 6 * COORDINATES] = 
              C1minus * (0.5 * dh2 * axym);
          } else {
            double bnd_utility = UPPER_BOUNDARY_FACTOR * RUNNING_COST(blc_elm_up, control1);
            bnd_utility *= 0.5 * dh2 * axym;
            lis_vector_set_value(LIS_ADD_VALUE, idx, bnd_utility, b);
          }


          // First sub diagonal is the probability to jump to (-1, 0).
          if ((idx % STEPS) > 0) {
            prob_matrix_data[idx + 3 * COORDINATES] =
              C1minus * (dh * pnl_b_m + 0.5 * dh2 * axx - 0.5 * dh2 * axya);
          } else {
            double bnd_utility = LEFT_BOUNDARY_FACTOR(time) * running_cost;
            bnd_utility *= dh * pnl_b_m + 0.5 * dh2 * axx - 0.5 * dh2 * axya;
            lis_vector_set_value(LIS_ADD_VALUE, idx, bnd_utility, b);
          }


          // First sup diagonal is the probability to jump to (1, 0).
          if ((idx % STEPS) < (STEPS - 1)) {
            prob_matrix_data[idx + 5 * COORDINATES] =
              C1minus * (dh * pnl_b_p + 0.5 * dh2 * axx - 0.5 * dh2 * axya);
          } else {
            double bnd_utility = RIGHT_BOUNDARY_FACTOR(time) * running_cost;
            bnd_utility *= dh * pnl_b_p + 0.5 * dh2 * axx - 0.5 * dh2 * axya;
            lis_vector_set_value(LIS_ADD_VALUE, idx, bnd_utility, b);
          }

          // Finished building the probability matrix. 

          // Add the running cost to the vector b.

          lis_vector_set_value(LIS_ADD_VALUE, idx, running_cost, b);
              //exp(-GAMMA_TIME * time) * pow(fmax(blc_elm, 0.0), GAMMA) - control1, b);

          // Finished constructing the first part the vector b.
        } 

//      for (idx = 0; idx < COORDINATES; idx++) 
//        cerr << idx << ": "
//             << prob_matrix_data[idx] << " "  
//             << prob_matrix_data[idx + COORDINATES] << " " 
//             << prob_matrix_data[idx + 2 * COORDINATES] << " " 
//             << prob_matrix_data[idx + 3 * COORDINATES] << " " 
//             << prob_matrix_data[idx + 4 * COORDINATES] << " " 
//             << prob_matrix_data[idx + 5 * COORDINATES] << " " 
//             << prob_matrix_data[idx + 6 * COORDINATES] << " "
//             << prob_matrix_data[idx + 7 * COORDINATES] << " "
//             << prob_matrix_data[idx + 8 * COORDINATES] << "      "
//             << prob_matrix_data[idx] 
//                + prob_matrix_data[idx + COORDINATES]  
//                + prob_matrix_data[idx + 2 * COORDINATES]  
//                + prob_matrix_data[idx + 3 * COORDINATES]  
//                + prob_matrix_data[idx + 4 * COORDINATES]  
//                + prob_matrix_data[idx + 5 * COORDINATES]  
//                + prob_matrix_data[idx + 6 * COORDINATES]  
//                + prob_matrix_data[idx + 7 * COORDINATES]  
//                + prob_matrix_data[idx + 8 * COORDINATES]  << endl;

        // TODO: Multiply the vector b with C1minus and the time scale.
        lis_vector_scale(-C1minus * DELTA, b);
        // Dump the vector b screen.
//        for (idx = 0; idx < COORDINATES; idx++) {
//          double b_elm;
//          lis_vector_get_value(b, idx, &b_elm);
//          cerr << b_elm << " ";
//        }
//        cerr << endl;

        // Solve the system of linear equations.
        int err = lis_solve(prob_matrix, b, cost_to_go, solver);
        
        if (err != 0) {
          cerr << "Got error " << err << " from lis solver." << endl;
        }
        
//        // Dump the vector cost_to_go screen.
//        for (idx = 0; idx < COORDINATES; idx++) {
//          double cost_to_go_elm;
//          lis_vector_get_value(cost_to_go, idx, &cost_to_go_elm);
//          cerr << cost_to_go_elm << " ";
//        }
//        cerr << endl;

        // Copy the cost to go the optimal value if the maximum of the cost to go
        // sub vector is larger than the optimal value found until now. 
        // The sub vectors consist of the vector minus the boundaries.
        //cerr << "@ time " << time << " control1 " << control1 << " control2 " << control2 << endl;
        for (idx = 0; idx < COORDINATES; idx++) {
          double val_new_elm;
          double cost_to_go_elm;
          lis_vector_get_value(cost_to_go, idx, &cost_to_go_elm);
          lis_vector_get_value(val_new, idx, &val_new_elm);
          // cerr << idx << " " << val_new_elm << " " << cost_to_go_elm << endl;

          if (val_new_elm < cost_to_go_elm) {
            lis_vector_set_value(LIS_INS_VALUE, idx, cost_to_go_elm, val_new);
            opt_control1[idx] = control1;
            opt_control2[idx] = control2;
          }
        }
      }
    }

    // Dump optimal controls to file.
    for (idx = 0; idx < COORDINATES; idx++) {
      *ctl1_file << opt_control1[idx] << " ";
      *ctl2_file << opt_control2[idx] << " ";

      double val_new_elm;
      lis_vector_get_value(val_new, idx, &val_new_elm);
      *opt_val_end_file << val_new_elm << " ";
    }

    *ctl1_file << endl;
    *ctl2_file << endl;
    *opt_val_end_file << endl;

#ifdef HAVE_BOOST
    // Increase the progress bar.
    ++show_progress;
#endif // HAVE_BOOST

    lis_vector_copy(val_new, val_old);
  }

  // Free all the data structures used in the solver.
  lis_vector_destroy(cost_to_go);
  lis_vector_destroy(b);
  lis_vector_destroy(val_new);
  lis_vector_destroy(val_old);
  lis_matrix_destroy(prob_matrix);

  delete opt_control1;
  delete opt_control2;
}

  
/* vim:set spell spelllang=en_gb cindent foldmethod=syntax tw=80 ts=2 et: */

