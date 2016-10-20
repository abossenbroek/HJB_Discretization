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

#include <sysexits.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <string>
#include <iostream>
#include <fstream>
using namespace std;

#include "markovchain.h"
#include "markovsolverlib/problem.h"
#include "markovsolverlib/solver.h"
#include "problemlua.h"
#include "problemcalloption.h"
#include "problemindiffpricing.h"
#include "solverexplicitgm.h"
#include "solverimplicitbm.h"

int main(int argc, char *argv[])
{
    string  cfg_file("");
    string  prblm_file("");
    bool    use_explicit;
    Solver  *s;
    Problem *p = 0;
    int     problem_type = PROBLEMLUA;

    po::options_description cmdline("Command line");
    cmdline.add_options()
        ("help", "produce help message")
        ("problemfile", po::value<string>(), "lua script to use as a problem")
        ("useexplicit", po::value<bool>(), "use an explicit solver")
        ("prblm", po::value<string>(), "use an explicit solver (indiffpricing|lua)")
        ;
    
    po::variables_map vm;

    po::store(po::parse_command_line(argc, argv, cmdline), vm);

    po::notify(vm);    


    if (vm.count("help")) {
      cout << cmdline << "\n";
      return EX_USAGE;
    }

    if (vm.count("problemfile")) {
        prblm_file = vm["problemfile"].as<string>();
    } else {
      cerr << "A lua file is required to read the problem" << endl;
      return  EX_USAGE;
    }
    if (vm.count("useexplicit")) {
        use_explicit = vm["useexplicit"].as<bool>();
    } else {
      use_explicit = false;
    }
    if (vm.count("prblm")) {
      if ((vm["prblm"].as<string>()).compare("indiffpricing") == 0)
        problem_type = PROBLEMINDIFFPRICING;
      else if ((vm["prblm"].as<string>()).compare("lua") == 0)
        problem_type = PROBLEMLUA;
      else {
        cerr << "Unkown problem " << vm["problemfile"].as<string>() << endl;
        throw new exception();
      }
    }

    switch (problem_type) {
      case PROBLEMLUA:
        cout << "Loading LUA problem" << endl;
        p = new ProblemLua(prblm_file);
        break;

      case PROBLEMINDIFFPRICING:
        cout << "Loading Indiffpricing problem" << endl;
        p = new ProblemIndiffpricing(prblm_file);
        break;

      default:
        cerr << "Problem type could not be found" << endl;
        throw new exception();
    }

    //ProblemCallOption p(prblm_file);
    cout << *p << endl;
    if (use_explicit)
      s = new SolverExplicitGM(p);
    else
      s = new SolverImplicitBM(p);
    
    s->solve();

    delete s;
    delete p;
};

/* vim:set spell spelllang=en_gb cindent foldmethod=syntax tw=80 ts=2 et: */
