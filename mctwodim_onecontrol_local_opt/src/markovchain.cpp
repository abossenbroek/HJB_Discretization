/*
 * =====================================================================================
 *
 *       Filename:  markovchain.cpp
 *
 *    Description:  Compute the strategy which is required to attain the
 *									maximal utility considering Merton's problem.
 *
 *        Version:  0.01
 *        Created:  07/01/2009 16:15:39
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

#include <lis.h>

#include <sysexits.h>

#include <string>
#include <iostream>
#include <fstream>
using namespace std;

#include "markovchain.h"
#include "solver.h"

int main(int argc, char *argv[])
{
    // Initialize lis with a dummy problem.
    lis_initialize(&argc, &argv);

    Solver s("control1", "control2", "optimal_val_end");

    s.solve();

    lis_finalize();
};

/* vim:set spell spelllang=en_gb cindent foldmethod=syntax tw=80 ts=2 et: */
