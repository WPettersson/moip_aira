#include <atomic>
#include <condition_variable>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <mutex>
#include <thread>

#include <ilcplex/cplex.h>

#include <boost/program_options.hpp>

#include "cluster.h"
#include "problem.h"
#include "env.h"
#include "lockingvars.h"
#include "symgroup.h"
#include "result.h"
#include "solutions.h"
#include "errors.h"
#include "thread.h"

//#define DEBUG
//#define FINETIMING

#ifdef DEBUG
std::mutex debug_mutex;
#endif

int num_threads;
int cplex_threads;

/**
 * Are we splitting up the range of values an objective can take, such that
 * individual threads only find solutions in their own pre-determined range?
 */
bool split;

/**
 * Has a thread completed its search. If one thread is finished, it will have
 * synchronised along the way, and therefore it (and possibly its
 * cluster-friends) will have the complete list of solutions. That is, if
 * completed is true, then the search may exit.
 */
bool completed;

/**
 * Track how many individual IPs we solve.
 */
std::atomic<int> ipcount;


/**
 * Solves a CLMOIP and returns the status and result.
 * Variables:
 *
 * \param e An Env object holding an environment suitable for solving IPs
 * \param p The Problem object to be solved
 * \param result An array of ints of size p.objcnt where the result will be
 * stored
 * \param rhs An array of doubles holding the right hand side of the problem
 * \param perm_id The index of the permutation which denotes the hierarchy or
 * ordering of the objectives, from most significant to least.
 */
int solve(Env & e, Problem & p, int * result, double * rhs, Thread * t);
int solve(Env & e, Problem & p, int * result, double * rhs, const int * perm);


/**
 * Finds out the limit (either maximum or minimum) of one objective in the
 * problem. This is only used to set up the splitting algorithm correctly.
 */
int get_limit(Env & e, Problem & p, int obj, double * rhs, const Sense sense);

/**
 * Optimise!
 * This function runs through the selected algorithm of Pettersson and
 * Ozlen[1]. It is designed to be run multi-threaded, hence using locking
 * mechanisms.
 *
 * \pram pFilename The file name of the file holding the problem description
 * \param all A Solutions object into which all solutions will be placed
 * \param t A pointer to a Thread object that describes how this particular
 * thread should approach the problem.
 **/

/* Note that this function is templated. There's no evidence that this does
 * increase running time, but it's not like it decreases it either. */
template<Sense sense>
void optimise(const char * pFilename, Solutions & all, Thread *t);

namespace po = boost::program_options;

int main (int argc, char *argv[])
{

  int status = 0; /* Operation status */

  Env e;

  std::string pFilename, outputFilename;

  /* Timing */
  clock_t starttime, endtime;
  double cpu_time_used, elapsedtime, startelapsed;
  int solcount;

  bool spread = false;
  split = false;

  po::variables_map v;
  po::options_description opt("Options for aira");
  opt.add_options()
    ("help,h", "Show this help.")
    ("lp,p",
      po::value<std::string>(&pFilename),
     "The LP file to solve. Required.")
    ("output,o",
      po::value<std::string>(&outputFilename),
     "The output file. Optional.")
    ("split,",
     po::bool_switch(&split),
     "Split the range of the first objective into one strip per thread\n"
     "Optional, defaults to False.")
    ("spread,s",
     po::bool_switch(&spread),
     "Spread threads out over various subgroups of the symmetries, or cluster inside subgroups.\n"
     "Option, defaults to False")
    ("threads,t",
      po::value<int>(&num_threads)->default_value(1),
     "Number of threads to use internally. Optional, default to 1.")
    ("cplex_threads,c",
        po::value<int>(&cplex_threads)->default_value(1),
     "Number of threads to allocate to CPLEX.\n"
     "Note that each internal thread calls CPLEX, so the total number of "
     "threads used is threads*cplex_threads.\n"
     "Optional, defaults to 1.")
  ;

  po::store(po::parse_command_line(argc, argv, opt), v);
  po::notify(v);

  if (v.count("help")) {
    // usage();
    std::cout << opt << std::endl;
    return(1);
  }


  if (v.count("lp") == 0) {
    // usage();
    std::cerr << "Error: You must pass in a problem file. All other parameters are optional." << std::endl;
    std::cout << opt << std::endl;
    return(1);
  }




  std::ofstream outFile;

  if (v.count("output")) {
    outFile.open(outputFilename);
  } else {
    /* Output file: <filename>.out */
    size_t last = pFilename.find_last_of(".");
    outFile.open(pFilename.substr(0, last).append(".out"));
  }



  /* Initialize the CPLEX environment */
  e.env = CPXopenCPLEX (&status);

  Problem p(pFilename.c_str(), e);

  /* Set to deterministic parallel mode */
  status=CPXsetintparam(e.env, CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC);

  /* Set to only one thread */
  CPXsetintparam(e.env, CPXPARAM_Threads, cplex_threads);

  if (e.env == NULL) {
    std::cerr << "Could not open CPLEX environment." << std::endl;
    return -ERR_CPLEX;
  }

  status = CPXsetintparam(e.env, CPX_PARAM_SCRIND, CPX_OFF);
  if (status) {
    std::cerr << "Failure to turn off screen indicator, error." << std::endl ;
      return -ERR_CPLEX;
  }

  outFile << std::endl << "Using improved algorithm" << std::endl;

  /* Start the timer */
  starttime = clock();
  timespec start;
  clock_gettime(CLOCK_MONOTONIC, &start);
  startelapsed = start.tv_sec + start.tv_nsec/1e9;


  if (num_threads > S[p.objcnt].size())
    num_threads = S[p.objcnt].size();

  Solutions all(p.objcnt);

  completed = false;

  std::list<Locking_Vars*> locking_var_list;
  double start_point, stop_point;
  if (split) {
    if (p.objsen == MIN) {
      start_point = get_limit(e, p, p.objcnt-1, p.rhs, MAX);
    } else {
      start_point = get_limit(e, p, p.objcnt-1, p.rhs, MIN);
    }
    stop_point = get_limit(e, p, p.objcnt-1, p.rhs, p.objsen);
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread m found start point " << start_point;
          std::cout << " and stop point " << stop_point << std::endl;
          debug_mutex.unlock();
#endif
  }
  /* Free up memory as all we needed was a count of the number of objectives. */
  if ( e.lp != NULL ) {
     status = CPXfreeprob (e.env, &e.lp);
     if ( status ) {
       std::cerr << "CPXfreeprob failed." << std::endl;
     }
  }
  if ( e.env != NULL ) {
     status = CPXcloseCPLEX (&e.env);

     if ( status ) {
       std::cerr << "Could not close CPLEX environment." << std::endl;
     }
  }

  std::list<Thread*> threads;
  if (split) {
    double split_start = start_point;
    double split_stop;
    double step_size = (stop_point - start_point)/num_threads;
    for (int t = 0; t < num_threads; ++t) {
      split_stop = split_start + step_size;
      threads.push_back(new Thread(t, p.objcnt, split_start, split_stop));
      split_start = split_stop;
    }
  } else {
    // Not splitting.
    int * ordering = new int[p.objcnt];
    int ** share_from = new int*[p.objcnt] {nullptr};
    int ** share_to = new int*[p.objcnt] {nullptr};
    int ** share_bounds = new int*[p.objcnt] {nullptr};
    int ** share_limit = new int*[p.objcnt] {nullptr};
    Locking_Vars ** lvs = new Locking_Vars*[p.objcnt] {nullptr};
    Cluster c(num_threads, p.objcnt, p.objsen, spread, p.objcnt, ordering,
        share_from, share_to, share_bounds, share_limit, threads, lvs);
    delete[] ordering;
    delete[] share_to;
    delete[] share_from;
    delete[] share_bounds;
    delete[] share_limit;
  }
  std::list<std::thread> threadList;
  if (p.objsen == MIN) {
    for(Thread* thread: threads) {
      threadList.emplace_back(optimise<MIN>, pFilename.c_str(),
          std::ref(all), thread);
    }
  } else {
    for(Thread* & thread: threads) {
      threadList.emplace_back(optimise<MAX>, pFilename.c_str(),
          std::ref(all), thread);
    }
  }
  for (auto& thread: threadList)
    thread.join();

  for (auto lv: locking_var_list)
    delete lv;

  /* Stop the clock. Sort and print results.*/
  endtime = clock();
  cpu_time_used=((double) (endtime - starttime)) / CLOCKS_PER_SEC;
  clock_gettime(CLOCK_MONOTONIC, &start);
  elapsedtime = (start.tv_sec + start.tv_nsec/1e9 - startelapsed);

  all.sort();
  solcount = 0;

  /* List of things we've printed */
  std::list<const Result *> printed;
  for (const Result* soln: all) {
    if (soln->infeasible)
      continue;
    bool printMe = true;
    for (auto& done: printed) {
      if (*done == *soln) {
        printMe = false;
        break;
      }
    }
    if (printMe) {
      for (int i = 0; i < p.objcnt; i++) {
        outFile << soln->result[i] << "\t";
      }
      outFile << std::endl;
      printed.push_back(soln);
      solcount++;
    }
  }

  constexpr int width = 8;
  constexpr int precision = 3;
  outFile << std::endl << "---" << std::endl;
  outFile << std::setw(width) << std::setprecision(precision) << std::fixed;
  outFile << cpu_time_used << " CPU seconds" << std::endl;
  outFile << std::setw(width) << std::setprecision(precision) << std::fixed;
  outFile << elapsedtime << " elapsed seconds" << std::endl;
  outFile << std::setw(width) << std::setprecision(precision) << std::fixed;
  outFile << ipcount << " IPs solved" << std::endl;
  outFile << std::setw(width) << std::setprecision(precision) << std::fixed;
  outFile << solcount << " Solutions found" << std::endl;


  outFile.close();

  return (status);
}

/* Solve CLMOIP and return solution status */
int solve(Env & e, Problem & p, int * result, double * rhs, const int * perm) {

  int cur_numcols, status, solnstat;
  double objval;
  double * srhs;
  srhs = new double[p.objcnt];

  memcpy(srhs, rhs, p.objcnt * sizeof(double));

  cur_numcols = CPXgetnumcols(e.env, e.lp);

  for (int j_preimage = 0; j_preimage < p.objcnt; j_preimage++) {
    int j = perm[j_preimage];
    status = CPXchgobj(e.env, e.lp, cur_numcols, p.objind[j], p.objcoef[j]);
    if (status) {
      std::cerr << "Failed to set objective." << std::endl;
    }

    status = CPXchgrhs (e.env, e.lp, p.objcnt, p.conind, srhs);
    if (status) {
      std::cerr << "Failed to change constraint srhs" << std::endl;
    }

    /* solve for current objective*/
    status = CPXmipopt (e.env, e.lp);
    if (status) {
      std::cerr << "Failed to optimize LP." << std::endl;
    }

    // This is shared across threads, but it's an atomic integer (so read/writes
    // are atomic) and thus we don't need a lock/mutex
    ipcount++;

    solnstat = CPXgetstat (e.env, e.lp);
    if ((solnstat == CPXMIP_INFEASIBLE) || (solnstat == CPXMIP_INForUNBD)) {
       break;
    }
    status = CPXgetobjval (e.env, e.lp, &objval);
    if ( status ) {
      std::cerr << "Failed to obtain objective value." << std::endl;
      exit(0);
    }
    if ( objval > 1/p.mip_tolerance ) {
      while (objval > 1/p.mip_tolerance) {
        p.mip_tolerance /= 10;
      }
      CPXsetdblparam(e.env, CPXPARAM_MIP_Tolerances_MIPGap, p.mip_tolerance);
      status = CPXmipopt (e.env, e.lp);
      ipcount++;
      solnstat = CPXgetstat (e.env, e.lp);
      if ((solnstat == CPXMIP_INFEASIBLE) || (solnstat == CPXMIP_INForUNBD)) {
        break;
      }
      status = CPXgetobjval (e.env, e.lp, &objval);
      if ( status ) {
        std::cerr << "Failed to obtain objective value." << std::endl;
        exit(0);
      }
    }

    //p.result[j] = srhs[j] = round(objval);
    result[j] = srhs[j] = round(objval);
  }

  delete[] srhs;

  return solnstat;
}

int get_limit(Env & e, Problem & p, int obj, double * rhs, const Sense sense) {

  int cur_numcols, status, solnstat;
  double objval;
  double * srhs;
  srhs = new double[p.objcnt];

  memcpy(srhs, rhs, p.objcnt * sizeof(double));

  cur_numcols = CPXgetnumcols(e.env, e.lp);

  status = CPXchgobj(e.env, e.lp, cur_numcols, p.objind[obj], p.objcoef[obj]);
  if (status) {
    std::cerr << "Failed to set objective." << std::endl;
  }

  status = CPXchgrhs (e.env, e.lp, p.objcnt, p.conind, srhs);
  if (status) {
    std::cerr << "Failed to change constraint srhs" << std::endl;
  }

  int i_sense;
  if (sense == MIN) {
    i_sense = 1;
  } else {
    i_sense = -1;
  }
  status = CPXchgobjsen (e.env, e.lp, i_sense);
  if (status) {
    std::cerr << "Failed to change constraint srhs" << std::endl;
  }

  /* solve for current objective*/
  status = CPXmipopt (e.env, e.lp);
  if (status) {
    std::cerr << "Failed to optimize LP." << std::endl;
  }

  // This is shared across threads, but it's an atomic integer (so read/writes
  // are atomic) and thus we don't need a lock/mutex
  ipcount++;

  solnstat = CPXgetstat (e.env, e.lp);
  if ((solnstat == CPXMIP_INFEASIBLE) || (solnstat == CPXMIP_INForUNBD)) {
    if (sense == MIN) {
      return (int)CPX_INFBOUND;
    } else {
      return (int)-CPX_INFBOUND;
    }
  }
  status = CPXgetobjval (e.env, e.lp, &objval);
  if ( status ) {
    std::cerr << "Failed to obtain objective value." << std::endl;
    exit(0);
  }
  if ( objval > 1/p.mip_tolerance ) {
    while (objval > 1/p.mip_tolerance) {
      p.mip_tolerance /= 10;
    }
    CPXsetdblparam(e.env, CPXPARAM_MIP_Tolerances_MIPGap, p.mip_tolerance);
    status = CPXmipopt (e.env, e.lp);
    ipcount++;
    solnstat = CPXgetstat (e.env, e.lp);
    if ((solnstat == CPXMIP_INFEASIBLE) || (solnstat == CPXMIP_INForUNBD)) {
      if (sense == MIN) {
        return (int)CPX_INFBOUND;
      } else {
        return (int)-CPX_INFBOUND;
      }
    }
    status = CPXgetobjval (e.env, e.lp, &objval);
    if ( status ) {
      std::cerr << "Failed to obtain objective value." << std::endl;
      exit(0);
    }
  }
  delete[] srhs;

  return (int)objval;
}

int solve(Env & e, Problem & p, int * result, double * rhs, Thread * t) {

  int cur_numcols, status, solnstat;
  double objval;
  double * srhs;
  srhs = new double[p.objcnt];

  //for(int i = 0; i < p.objcnt; ++i)
  //  srhs[i] = rhs[perm[i]];

  memcpy(srhs, rhs, p.objcnt * sizeof(double));

  cur_numcols = CPXgetnumcols(e.env, e.lp);

  // TODO Permutation applies here.
  for (int j_preimage = 0; j_preimage < p.objcnt; j_preimage++) {
    int j = t->perm(j_preimage);
    status = CPXchgobj(e.env, e.lp, cur_numcols, p.objind[j], p.objcoef[j]);
    if (status) {
      std::cerr << "Failed to set objective." << std::endl;
    }

    status = CPXchgrhs (e.env, e.lp, p.objcnt, p.conind, srhs);
    if (status) {
      std::cerr << "Failed to change constraint srhs" << std::endl;
    }

    /* solve for current objective*/
    status = CPXmipopt (e.env, e.lp);
    if (status) {
      std::cerr << "Failed to optimize LP." << std::endl;
    }

    // This is shared across threads, but it's an atomic integer (so read/writes
    // are atomic) and thus we don't need a lock/mutex
    ipcount++;

    solnstat = CPXgetstat (e.env, e.lp);
    if ((solnstat == CPXMIP_INFEASIBLE) || (solnstat == CPXMIP_INForUNBD)) {
       break;
    }
    status = CPXgetobjval (e.env, e.lp, &objval);
    if ( status ) {
      std::cerr << "Failed to obtain objective value." << std::endl;
      exit(0);
    }
    if ( objval > 1/p.mip_tolerance ) {
      while (objval > 1/p.mip_tolerance) {
        p.mip_tolerance /= 10;
      }
      CPXsetdblparam(e.env, CPXPARAM_MIP_Tolerances_MIPGap, p.mip_tolerance);
      status = CPXmipopt (e.env, e.lp);
      ipcount++;
      solnstat = CPXgetstat (e.env, e.lp);
      if ((solnstat == CPXMIP_INFEASIBLE) || (solnstat == CPXMIP_INForUNBD)) {
        break;
      }
      status = CPXgetobjval (e.env, e.lp, &objval);
      if ( status ) {
        std::cerr << "Failed to obtain objective value." << std::endl;
        exit(0);
      }
    }

    //p.result[j] = srhs[j] = round(objval);
    result[j] = srhs[j] = round(objval);
  }

  delete[] srhs;

  return solnstat;
}

template<Sense sense>
void optimise(const char * pFilename, Solutions & all, Thread *t) {
  Env e;
  const bool sharing = (t->share_to != nullptr);
#ifdef DEBUG
  debug_mutex.lock();
  if (split) {
    std::cout << "Thread " << t->id << " has start ";
    std::cout << t->split_start << " and stop " << t->split_stop << std::endl;
  }
  if (sharing) {
    std::cout << "Thread " << t->id << " is sharing" << std::endl;
  }
  debug_mutex.unlock();
#endif
#ifdef FINETIMING
  double cplex_time = 0;
  double wait_time = 0;
  timespec start;
  clock_gettime(CLOCK_MONOTONIC, &start);
  double total_time = start.tv_sec + start.tv_nsec/1e9;
#endif
  {
    int status; // Only need status here.
    /* Initialize the CPLEX environment */
    e.env = CPXopenCPLEX (&status);
    if (status != 0) {
      std::cerr << "Could not open CPLEX environment." << std::endl;
    }

    /* Set to deterministic parallel mode */
    status=CPXsetintparam(e.env, CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC);

    /* Set to only one thread */
    CPXsetintparam(e.env, CPXPARAM_Threads, cplex_threads);

    if (e.env == NULL) {
      std::cerr << "Could not open CPLEX environment." << std::endl;
    }

    status = CPXsetintparam(e.env, CPX_PARAM_SCRIND, CPX_OFF);
    if (status) {
      std::cerr << "Failure to turn off screen indicator." << std::endl;
    }
  }

  Problem p(pFilename, e);

  Solutions s(p.objcnt);

  int infcnt;
  bool inflast;
  bool infeasible;
  int * max, * min,  * result, *resultStore;
  double * rhs;

  result = resultStore = new int[p.objcnt];
  rhs = new double[p.objcnt];
  for(int i = 0; i < p.objcnt; ++i) {
    rhs[i] = p.rhs[i];
  }
  if (split) {
    rhs[t->perm(p.objcnt-1)] = t->split_start;
  }

#ifdef FINETIMING
  clock_gettime(CLOCK_MONOTONIC, &start);
  double starttime = (start.tv_sec + start.tv_nsec/1e9);
#endif
  int solnstat = solve(e, p, result, rhs, t);
#ifdef FINETIMING
  clock_gettime(CLOCK_MONOTONIC, &start);
  cplex_time += (start.tv_sec + start.tv_nsec/1e9) - starttime;
#endif
#ifdef DEBUG
  debug_mutex.lock();
  std::cout << "Thread " << t->id << " with constraints ";
  for(int i = 0; i < p.objcnt; ++i) {
    if (rhs[i] > 1e09)
      std::cout << "∞";
    else if (rhs[i] < -1e09)
      std::cout << "-∞";
    else
      std::cout << rhs[i];
    std::cout << ",";
  }
  std::cout << " found ";
  for(int i = 0; i < p.objcnt; ++i) {
    std::cout << result[i] << ",";
  }
  std::cout << std::endl;
  debug_mutex.unlock();
#endif

  /* Need to add a result to the list here*/
  s.insert(rhs, result, solnstat == CPXMIP_INFEASIBLE);
  // Note that if we are splitting, we aren't sharing.
  if (!split) {
    if (sharing && (solnstat != CPXMIP_INFEASIBLE)) {
      int *objectives = new int[p.objcnt];
      for (int i = 0; i < p.objcnt; ++i) {
        objectives[i] = result[i];
      }
    }
  }
  min = new int[p.objcnt];
  max = new int[p.objcnt];

#ifdef DEBUG
  debug_mutex.lock();
  std::cout << "Thread " << t->id << " using permutation ";
  for (int i = 0; i < p.objcnt; ++i) {
    std::cout << t->perm(i) <<  ", ";
  }
  std::cout << std::endl;
  debug_mutex.unlock();
#endif

  if (sharing && (solnstat != CPXMIP_INFEASIBLE) && (p.objcnt > 1)) {
    int i = t->perm(1);
    if (t->share_to[i] != nullptr) {
      if (sense == MIN) {
        if (*t->share_to[i] < result[i]) {
          *t->share_to[i] = result[i];
        }
      } else {
        if (*t->share_to[i] > result[i]) {
          *t->share_to[i] = result[i];
        }
      }
    }
  }
  if (solnstat != CPXMIP_INFEASIBLE) {
    for (int j = 0; j < p.objcnt; j++) {
      min[j] = max[j] = result[j];
    }
  }


  /**
   * objective_counter marks which objectives are currently active. Note that
   * the first objective (highest priority) is 0, which we never make active.
   *
   * That is, if the current permutation is 1,2,3,4,5 and objective_counter is
   * 1, then we are trying to optimise objective 2. In other words, objectives
   * 3,4,5 all have their limits set to +- infinity, and these won't change
   * until objective_counter increases.
   *
   * If objective_counter is 3, then we are changing objectives 2,3 and 4.
   */

  /**
   * depth tracks which of the active objectives is the one being currently
   * modified. If objective_counter is 3, then the active objectives are 2,3
   * and 4. If depth is 1, then objective 2 is the one for which we are
   * changing the RHS.
   */

  /**
   * onwalk is true if the current active depth (as indicated by the depth
   * variable) was just increased.
   */
  for (int objective_counter = 1; objective_counter < p.objcnt; objective_counter++) {
    int objective = t->perm(objective_counter);
    int depth_level = 1; /* Track current "recursion" depth */
    int depth = t->perm(depth_level); /* Track depth objective */
    bool onwalk = false; /* Are we on the move? */
    //bool * found_any = new bool[p.objcnt] {false};
    infcnt = 0; /* Infeasible count*/
    inflast = false; /* Last iteration infeasible?*/

    /* Set all constraints back to infinity*/
    for (int j_pre = 1; j_pre < p.objcnt; j_pre++) {
      int j = t->perm(j_pre);
      if ((!sharing) || (t->share_from[j] == nullptr)) {
        if (sense == MIN) {
          rhs[j] = CPX_INFBOUND;
        } else {
          rhs[j] = -CPX_INFBOUND;
        }
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << j << "] to ±inf" << std::endl;
          debug_mutex.unlock();
#endif
      } else {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << j << "] to " << *t->share_from[j] << std::endl;
          debug_mutex.unlock();
#endif
        rhs[j] = *t->share_from[j];
      }
    }
    if (split) {
        rhs[t->perm(p.objcnt-1)] = t->split_start;
    }
    /* Set rhs of current depth */
    if (sense == MIN) {
      rhs[objective] = max[objective]-1;
    } else {
      rhs[objective] = min[objective]+1;
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << objective << "] to " << (min[objective]+1) << std::endl;
          debug_mutex.unlock();
#endif
    }
    max[objective] = (int) -CPX_INFBOUND;
    min[objective] = (int) CPX_INFBOUND;
    while (infcnt < objective_counter && !completed) {
      bool relaxed;
      int solnstat;
      /* Look for possible relaxations to the current problem*/
      const Result *relaxation;

#ifdef DEBUG_SOLUTION_SEARCH
      debug_mutex.lock();
      std::cout << "Thread " << t->id;
      debug_mutex.unlock();
#endif
      relaxation = s.find(rhs, p.objsen);
      relaxed = (relaxation != nullptr);
      if (relaxed) {
        infeasible = relaxation->infeasible;
        result = relaxation->result;
      } else {
        /* Solve in the absence of a relaxation*/
        result = resultStore;
#ifdef FINETIMING
        clock_gettime(CLOCK_MONOTONIC, &start);
        double starttime = (start.tv_sec + start.tv_nsec/1e9);
#endif
        solnstat = solve(e, p, result, rhs, t);
#ifdef FINETIMING
        clock_gettime(CLOCK_MONOTONIC, &start);
        cplex_time += (start.tv_sec + start.tv_nsec/1e9) - starttime;
#endif
        infeasible = ((solnstat == CPXMIP_INFEASIBLE) || (solnstat == CPXMIP_INForUNBD));
        /* Store result */
        s.insert(rhs, result, infeasible);
      }
#ifdef DEBUG
      debug_mutex.lock();
      std::cout << "Thread " << t->id << " with constraints ";
      for(int i = 0; i < p.objcnt; ++i) {
        if (rhs[i] > 1e09)
          std::cout << "∞";
        else if (rhs[i] < -1e09)
          std::cout << "-∞";
        else
          std::cout << rhs[i];
        std::cout << ",";
      }
      std::cout << " found ";
      if (relaxed)
        std::cout << "relaxation ";
      if (infeasible) {
        std::cout << "infeasible." << std::endl;
      } else {
        for(int i = 0; i < p.objcnt; ++i) {
          std::cout << result[i] << ",";
        }
        std::cout << std::endl;
      }
      debug_mutex.unlock();
#endif
      if (split) {
        if (!infeasible) {
          // check if we cross midpoint if we are thread 1
          if (sense == MIN) {
            if (result[p.objcnt-1] < t->split_stop) {
              infeasible = true;
#ifdef DEBUG
              debug_mutex.lock();
              std::cout << "Thread " << t->id << " reached split_stop";
              std::cout << " " << t->split_stop << std::endl;
              debug_mutex.unlock();
#endif
            }
          } else {
            if (result[p.objcnt-1] > t->split_stop) {
              infeasible = true;
#ifdef DEBUG
              debug_mutex.lock();
              std::cout << "Thread " << t->id << " reached split_stop";
              std::cout << " " << t->split_stop << std::endl;
              debug_mutex.unlock();
#endif
            }
          }
          /* Update maxima */
          for (int j = 0; j < p.objcnt; j++) {
            if (result[j] > max[j]) {
              max[j] = result[j];
            }
          }
          /* Update minima */
          for (int j = 0; j < p.objcnt; j++) {
            if (result[j] < min[j]) {
              min[j] = result[j];
            }
          }
        }
        if (infeasible) {
          infcnt++;
          inflast = true;
        } else {
          infcnt = 0;
          inflast = false;
        }
      } else if (sharing && t->partnered && t->locks && t->locks[t->perm(infcnt+1)]) {
        std::unique_lock<std::mutex> lk(t->locks[t->perm(infcnt+1)]->status_mutex);
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << " got lock on " << t->perm(infcnt+1) << std::endl;
          debug_mutex.unlock();
#endif
        // Note that perm[1] is only shared with one other thread, that only
        // ever reads it, and that this thread always improves the value of
        // this shared limit, so we can use this faster update.
        if (!infeasible && (p.objcnt > 1)) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " updating first bound on " << t->perm(1);
          std::cout << " to " << result[t->perm(1)] << std::endl;
          debug_mutex.unlock();
#endif
          if (t->share_to[t->perm(1)] != nullptr)
            *t->share_to[t->perm(1)] = result[t->perm(1)];
        }

        if (!infeasible && (t->share_from[t->perm(0)] != nullptr)) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "result on " << t->perm(0) << " is " << result[t->perm(0)] << " and ";
          std::cout << "first bound on " << t->perm(0) << " is " << *t->share_from[t->perm(0)] <<std::endl;
          debug_mutex.unlock();
#endif
          if (sense == MIN) {
            if (result[t->perm(0)] >= *t->share_from[t->perm(0)]) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " smaller result found by partner, bailing." << std::endl;
          debug_mutex.unlock();
#endif
              // Check if our partner found something.
              if (t->locks) {
                int obj = t->perm(infcnt+1);
                Locking_Vars *lv = t->locks[obj];
                if (lv != nullptr) {
                  if (lv->found_any) {
#ifdef DEBUG
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " ";
            std::cout << "partner found something." << std::endl;
            debug_mutex.unlock();
#endif
                    // Partner found something
                    infcnt = 0;
                    inflast = true;
                    depth_level = 1;
                    depth = t->perm(depth_level);
                  }
                }
              }
              // Pretend infeasible to backtrack properly
              infeasible = true;
            }
          } else {
            if (result[t->perm(0)] <= *t->share_from[t->perm(0)]) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " bigger result found by partner, bailing." << std::endl;
          debug_mutex.unlock();
#endif
              // Check if our partner found something.
              if (t->locks) {
                int obj = t->perm(infcnt+1);
                Locking_Vars *lv = t->locks[obj];
                if (lv != nullptr) {
                  if (lv->found_any) {
#ifdef DEBUG
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " ";
            std::cout << "partner found something." << std::endl;
            debug_mutex.unlock();
#endif
                    // Partner found something
                    infcnt = 0;
                    inflast = true;
                    depth_level = 1;
                    depth = t->perm(depth_level);
                  }
                }
              }
              // Pretend infeasible to backtrack properly
              infeasible = true;
            }
          }
          // Duplicate code as we are marking this result infeasible
          /* Update maxima */
          for (int j = 0; j < p.objcnt; j++) {
            if (result[j] > max[j]) {
              max[j] = result[j];
            }
          }
          /* Update minima */
          for (int j = 0; j < p.objcnt; j++) {
            if (result[j] < min[j]) {
              min[j] = result[j];
            }
          }
        }
        if (!infeasible) {
          if (t->locks) {
            int obj = t->perm(infcnt+1);
            Locking_Vars *lv = t->locks[obj];
            if (lv != nullptr) {
              lv->found_any = true;
#ifdef DEBUG
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " ";
            std::cout << "found something at infcnt=" << infcnt << std::endl;
            debug_mutex.unlock();
#endif
            }
          }
          infcnt = 0;
          inflast = false;
          /* Update maxima */
          for (int j = 0; j < p.objcnt; j++) {
            if (result[j] > max[j]) {
              max[j] = result[j];
            }
          }
          /* Update minima */
          for (int j = 0; j < p.objcnt; j++) {
            if (result[j] < min[j]) {
              min[j] = result[j];
            }
          }
        }
        if (infeasible) {
          if (t->locks) {
            int obj = t->perm(infcnt+1);
            Locking_Vars *lv = t->locks[obj];
            if (lv != nullptr) {
              if (lv->found_any) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "apparently partner found something" << std::endl;
          debug_mutex.unlock();
#endif
                infcnt = 0;
              }
            }
          }
          infcnt++;
          inflast = true;
        } else {
          infcnt = 0;
          inflast = false;
        }
      } else {
        if (infeasible) {
          infcnt++;
          inflast = true;
        } else {
          infcnt = 0;
          inflast = false;
          /* Update maxima */
          for (int j = 0; j < p.objcnt; j++) {
            if (result[j] > max[j]) {
              max[j] = result[j];
            }
          }
          /* Update minima */
          for (int j = 0; j < p.objcnt; j++) {
            if (result[j] < min[j]) {
              min[j] = result[j];
            }
          }
        }
      }

      // If we are sharing, and there are other threads, and this result
      // is infeasible, and it's a 2-objective problem
      if (sharing && infeasible && (infcnt+1) < p.objcnt) {
        int updated_objective = t->perm(infcnt+1);
#ifdef DEBUG
        debug_mutex.lock();
        std::cout << "Thread " << t->id <<  " done, ";
        std::cout << "infcnt is " << infcnt;
        std::cout << ", objective_counter is " << objective_counter;
        std::cout << ", updated_objective is " << updated_objective;
        if (sense == MIN) {
          std::cout << " and max[" << updated_objective << "] is " << max[updated_objective];
        } else {
          std::cout << " and min[" << updated_objective << "] is " << min[updated_objective];
        }
        std::cout << std::endl;
        debug_mutex.unlock();
#endif
        // If max/min have not been changed, then we have found zero new
        // solutions at this level. Normally this would only happen if
        // infcnt+1 == p.objcnt, but it can also happen if bounds on later
        // objectives are changed from other threads.
        if (sense == MIN) {
          if (max[updated_objective] == (int)-CPX_INFBOUND) {
            completed = true;
          }
        } else {
          if (min[updated_objective] == (int)CPX_INFBOUND) {
            completed = true;
          }
        }
        bool wait = false;
        Locking_Vars *lv = nullptr;
        if (t->locks && !completed)
          lv = t->locks[updated_objective];
        if ( lv != nullptr) {
          std::unique_lock<std::mutex> lk(lv->status_mutex);
          if ( lv->num_running_threads > 1) {
            lv->num_running_threads--; // This thread is no longer running.
#ifdef DEBUG_SYNC
            debug_mutex.lock();
            std::cout << "Thread " << t->id <<  " done, at " << __LINE__ << " w/ ";
            std::cout << lv->num_running_threads << " threads still going." << std::endl;
            debug_mutex.unlock();
#endif
            // Share all bounds!
            for(int pre_i=0; pre_i < p.objcnt; ++pre_i) {
              int i = t->perm(pre_i);
              if (t->share_bounds[i] != nullptr) {
                if (sense == MIN) {
                  if (*t->share_bounds[i] < max[i]) {
                    *t->share_bounds[i] = max[i];
                  }
                } else {
                  if (*t->share_bounds[i] > min[i]) {
                    *t->share_bounds[i] = min[i];
                  }
                }
              }
            }
            wait = true;
#ifdef FINETIMING
            clock_gettime(CLOCK_MONOTONIC, &start);
            double start_wait = start.tv_sec + start.tv_nsec/1e9;
#endif
            if (!completed) {
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "waiting on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
              lv->cv.wait(lk); // Wait for all other threads to update limits.
            }
#ifdef FINETIMING
            clock_gettime(CLOCK_MONOTONIC, &start);
            wait_time += (start.tv_sec + start.tv_nsec/1e9) - start_wait;
#endif
            // Share min/max values
            if (t->share_to[updated_objective] != nullptr) {
              if (sense == MIN) {
                if (*t->share_to[updated_objective] > max[updated_objective]) {
                  max[updated_objective] = *t->share_to[updated_objective];
                }
              } else {
                if (*t->share_to[updated_objective] < min[updated_objective]) {
                  min[updated_objective] = *t->share_to[updated_objective];
                }
              }
            }
            // Share bounds on all objectives
            for(int pre_i=0; pre_i < p.objcnt; ++pre_i) {
              int i = t->perm(pre_i);
              if (t->share_bounds[i] != nullptr) {
                if (sense == MIN) {
                  if (*t->share_bounds[i] > max[i]) {
                    max[i] = *t->share_bounds[i];
                  }
                } else {
                  if (*t->share_bounds[i] < min[i]) {
                    min[i] = *t->share_bounds[i];
                  }
                }
              }
            }
          } else {
            // This thread is last.
#ifdef DEBUG_SYNC
            debug_mutex.lock();
            std::cout << "Thread " << t->id <<  " done at " << __LINE__ << ", ";
            std::cout << "last in." << std::endl;
            debug_mutex.unlock();
#endif
            // Sharing bounds on all objectives
            for(int pre_i=0; pre_i < p.objcnt; ++pre_i) {
              int i = t->perm(pre_i);
              if (t->share_bounds[i] != nullptr) {
                if (sense == MIN) {
                  if (*t->share_bounds[i] < max[i]) {
                    *t->share_bounds[i] = max[i];
                  } else {
                    max[i] = *t->share_bounds[i];
                  }
                } else {
                  if (*t->share_bounds[i] > min[i]) {
                    *t->share_bounds[i] = min[i];
                  } else {
                    min[i] = *t->share_bounds[i];
                  }
                }
              }
            }
            // Sharing min/max
            if (t->share_to[updated_objective] != nullptr) {
              if (sense == MIN) {
                if (max[updated_objective] != (int)-CPX_INFBOUND) {
                  *t->share_to[updated_objective] = max[updated_objective];
                }
              } else {
                if (min[updated_objective] != (int)CPX_INFBOUND) {
                  *t->share_to[updated_objective] = min[updated_objective];
                }
              }
            }
            // Reset value of found_any for next time we visit this level.
            lv->found_any = false;
            // Update shared_limits for updated_objective
            int * limit = t->share_limit[updated_objective];
            if ((limit != nullptr) && (t->share_from[updated_objective] != nullptr)) {
              *limit = *t->share_from[updated_objective];
#ifdef DEBUG
              debug_mutex.lock();
              std::cout << "Thread " << t->id << " ";
              std::cout << "set limit[" << updated_objective << "] to ";
              std::cout << *limit << std::endl;
              debug_mutex.unlock();
#endif
            }
            lv->reset_num_running_threads(); // These threads are starting up again.
#ifdef DEBUG
            debug_mutex.lock();
            std::cout << "Thread " << t->id << ", ";
            std::cout << "max is " << max[updated_objective] << ", ";
            if (t->share_bounds[updated_objective] != nullptr)
              std::cout << "share_to is " << *t->share_to[updated_objective] << ", ";
            std::cout << lv->num_running_threads << " threads going." << std::endl;
            debug_mutex.unlock();
#endif
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "releasing on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
            lk.unlock();
            lv->cv.notify_all();
          }
        }
        if (!completed && lv != nullptr) {
          for(int i = 0; i <= infcnt; ++i) {
            if (sense == MIN) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "setting max[" << t->perm(i) << "] to -inf" << std::endl;
          debug_mutex.unlock();
#endif
              max[t->perm(i)] = (int)-CPX_INFBOUND;
              if (t->share_to[t->perm(i)] != nullptr) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "setting share_to[" << t->perm(i) << "] to inf" << std::endl;
          debug_mutex.unlock();
#endif
                *t->share_to[t->perm(i)] = (int)CPX_INFBOUND;
              }
              if (t->share_limit[t->perm(i)] != nullptr) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "setting limit_to[" << t->perm(i) << "] to inf" << std::endl;
          debug_mutex.unlock();
#endif
                *t->share_limit[t->perm(i)] = (int)CPX_INFBOUND;
              }
            } else {
              min[t->perm(i)] = (int)CPX_INFBOUND;
              if (t->share_to[t->perm(i)] != nullptr) {
                *t->share_to[t->perm(i)] = (int)-CPX_INFBOUND;
              }
              if (t->share_limit[t->perm(i)] != nullptr) {
                *t->share_limit[t->perm(i)] = (int)-CPX_INFBOUND;
              }
            }
          }

          // Reset share_bounds to their limit. First, reset to ±infinity.
            {
              std::unique_lock<std::mutex> lk(lv->status_mutex);
              lv->num_running_threads--;
              if (lv->num_running_threads > 0) {
                if (completed)
                  break;
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "waiting on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
                lv->cv.wait(lk);
              } else {
                // Last in. Reset the shared_bounds here.
                for(int pre_i = 0; pre_i <= infcnt+1; ++pre_i) {
                  int i = t->perm(pre_i);
                  if (t->share_bounds[i] != nullptr) {
                    if (sense == MIN) {
                      *t->share_bounds[i] = (int)-CPX_INFBOUND;
                    } else {
                      *t->share_bounds[i] = (int)CPX_INFBOUND;
                    }
                  }
                }
                lv->reset_num_running_threads();
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "releasing on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
                lk.unlock();
                lv->cv.notify_all();
              }
            }
          do {
            {
              std::unique_lock<std::mutex> lk(lv->status_mutex);
              lv->num_running_threads--;
              if (lv->num_running_threads > 0) {
                if (completed)
                  break;
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "waiting on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
                lv->cv.wait(lk);
              } else {
                // Last in
                lv->changed = false;
                lv->reset_num_running_threads();
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "releasing on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
                lk.unlock();
                lv->cv.notify_all();
              }
            }
            for (int i = 0; i <= infcnt; ++i) {
              int obj = t->perm(i);
              if (sense == MIN) {
                if (t->share_from[obj] != nullptr) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << " share_from[" << obj << "] is " << *t->share_from[obj] << std::endl;
          debug_mutex.unlock();
#endif
                  if (t->share_limit[obj] != nullptr) {
                    if (*t->share_limit[obj] > *t->share_from[obj]) {
                      lv->changed = true;
                      *t->share_limit[obj] = *t->share_from[obj];
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "setting share_limit[" << obj << "] to " << max[obj] << std::endl;
          debug_mutex.unlock();
#endif
                    }
                  }
                  if (t->share_to[obj] != nullptr) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << " share_to[" << obj << "] is " << *t->share_to[obj] << std::endl;
          debug_mutex.unlock();
#endif
                    if (*t->share_to[obj] > *t->share_from[obj]) {
                      lv->changed = true;
                      *t->share_to[obj] = *t->share_from[obj];
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "setting share_to[" << obj << "] to " << max[obj] << std::endl;
          debug_mutex.unlock();
#endif
                    }
                  }
                }
              } else {
                if (t->share_from[obj] != nullptr) {
                  if (t->share_limit[obj] != nullptr) {
                    if (*t->share_limit[obj] < *t->share_from[obj]) {
                      lv->changed = true;
                      *t->share_limit[obj] = *t->share_from[obj];
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "setting share_limit[" << obj << "] to " << max[obj] << std::endl;
          debug_mutex.unlock();
#endif
                    }
                  }
                  if (t->share_to[obj] != nullptr) {
                    if (*t->share_to[obj] < *t->share_from[obj]) {
                      lv->changed = true;
                      *t->share_to[obj] = *t->share_from[obj];
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "setting share_to[" << obj << "] to " << max[obj] << std::endl;
          debug_mutex.unlock();
#endif
                    }
                  }
                }
              }
            }
            {
              std::unique_lock<std::mutex> lk(lv->status_mutex);
              lv->num_running_threads--;
              if (lv->num_running_threads > 0) {
                if (completed)
                  break;
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "waiting on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
                lv->cv.wait(lk);
              } else {
                // Last in
                lv->reset_num_running_threads();
                lk.unlock();
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "releasing on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
                lv->cv.notify_all();
              }
            }
          } while(lv->changed);
#ifdef DEBUG
          debug_mutex.lock();
          for (int i=0; i < p.objcnt; ++i) {
            std::cout << "Thread " << t->id << " ";
            std::cout << "max[" << i << "] is " << max[i] << std::endl;
          }
          for (int i=0; i < p.objcnt; ++i) {
            std::cout << "Thread " << t->id << " ";
            std::cout << "min[" << i << "] is " << min[i] << std::endl;
          }
          for (int i=0; i < p.objcnt; ++i) {
            if (t->share_bounds[i]) {
              std::cout << "Thread " << t->id << " ";
              std::cout << "share_bounds[" << i << "] is " << *t->share_bounds[i] << std::endl;
            }
          }
          for (int i=0; i < p.objcnt; ++i) {
            if (t->share_from[i]) {
              std::cout << "Thread " << t->id << " ";
              std::cout << "share_from[" << i << "] is " << *t->share_from[i] << std::endl;
            }
          }
          for (int i=0; i < p.objcnt; ++i) {
            if (t->share_limit[i]) {
              std::cout << "Thread " << t->id << " ";
              std::cout << "share_limit[" << i << "] is " << *t->share_limit[i] << std::endl;
            }
          }
          for (int i=0; i < p.objcnt; ++i) {
            std::cout << "Thread " << t->id << " ";
            std::cout << "rhs[" << i << "] is " << rhs[i] << std::endl;
          }
          if (lv->found_any)
            std::cout << "Thread " << t->id << " found_any is true" << std::endl;
          std::cout << "Thread " << t->id << " ";
          std::cout << "synced." << std::endl;
          debug_mutex.unlock();
#endif
        }
      }

      if (sharing && (p.objcnt > 2) && (infcnt == objective_counter) && (infcnt == p.objcnt - 2)) {
        if (t->share_to[t->perm(p.objcnt-1)] == nullptr)
          continue;
#ifdef DEBUG
        debug_mutex.lock();
        std::cout << "Thread " << t->id << " updating bound on " << t->perm(p.objcnt-1);
        std::cout << " to ";
        if (sense == MIN) {
          std::cout << max[t->perm(p.objcnt-1)];
        } else {
          std::cout << min[t->perm(p.objcnt-1)];
        }
        std::cout << std::endl;
        debug_mutex.unlock();
#endif
        if (sense == MIN) {
          *t->share_to[t->perm(p.objcnt-1)] = max[t->perm(p.objcnt-1)];
        } else {
          *t->share_to[t->perm(p.objcnt-1)] = min[t->perm(p.objcnt-1)];
        }
      }
      if (infeasible && (infcnt == objective_counter-1)) {
        if ((p.objcnt > 2) && (objective_counter == p.objcnt - 1)) {
          if (t->share_to[t->perm(p.objcnt-1)] != nullptr) {
            if (sense == MIN) {
              *t->share_to[objective] = max[objective];
            } else {
              *t->share_to[objective] = min[objective];
            }
          }
        }
        /* Set all constraints back to infinity */
        for (int pre_j = 0; pre_j < p.objcnt; pre_j++) {
          int j  = t->perm(pre_j);
          if ((pre_j < infcnt) || (!sharing) || (t->share_limit[j] == nullptr && t->share_from[j] == nullptr)) {
            if (sense == MIN) {
              rhs[j] = CPX_INFBOUND;
            } else {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << j << "] to " << "-inf" << std::endl;
          debug_mutex.unlock();
#endif
              rhs[j] = -CPX_INFBOUND;
            }
          } else {
            int * share_from;
            if (t->share_limit[j] == nullptr)
              share_from = t->share_from[j];
            else
              share_from = t->share_limit[j];
            if (sense == MIN) {
#ifdef DEBUG
              debug_mutex.lock();
              std::cout << "Thread " << t->id << " setting rhs[" << j;
              std::cout << "] to " << (*share_from - 1) << std::endl;
              debug_mutex.unlock();
#endif
              rhs[j] = *share_from - 1;
            } else {
#ifdef DEBUG
              debug_mutex.lock();
              std::cout << "Thread " << t->id << " setting rhs[" << j;
              std::cout << "] to " << (*share_from + 1) << std::endl;
              debug_mutex.unlock();
#endif
              rhs[j] = *share_from + 1;
            }
            if (t->share_to[j] != nullptr) {
              if (sense == MIN) {
                if (*t->share_to[j] > *share_from) {
                  *t->share_to[j] = *share_from;
#ifdef DEBUG
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " setting share_to[" << j;
            std::cout << "] to " << *t->share_to[j] << std::endl;
            debug_mutex.unlock();
#endif
                }
              } else {
                if (*t->share_to[j] < *share_from) {
#ifdef DEBUG
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " changing share_to[" << j;
            std::cout << "] from " << *share_from << " to " << *t->share_to[j] << std::endl;
            debug_mutex.unlock();
#endif
                  *t->share_to[j] = *share_from;
                }
              }
            }
          }
        }
        /* In the case of a minimisation problem
         * set current level to max objective function value  -1 else set
         * current level to min objective function value  +1 */
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "setting rhs[" << objective << "] to " << max[objective]-1 <<std::endl;
          debug_mutex.unlock();
#endif
        if (sense == MIN) {
          rhs[objective] = max[objective]-1;
          max[objective] = (int) -CPX_INFBOUND;
        } else {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << objective << "] to " << (min[objective]+1) << std::endl;
          debug_mutex.unlock();
#endif
          rhs[objective] = min[objective]+1;
          min[objective] = (int) CPX_INFBOUND;
        }

        /* Reset depth */
        depth_level = 1;
        depth = t->perm(depth_level);
        onwalk = false;
      } else if (inflast && infcnt != objective_counter) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "increasing depth from " << depth_level << std::endl;
          debug_mutex.unlock();
#endif
        if (sense == MIN) {
          int * share_from = nullptr;
          if (t->share_limit[depth] != nullptr)
            share_from = t->share_limit[depth];
          else if (t->share_from[depth] != nullptr)
            share_from = t->share_from[depth];
          if (share_from != nullptr)
            rhs[depth] = *share_from - 1;
          else
            rhs[depth] = CPX_INFBOUND;
        } else {
          int * share_from = nullptr;
          if (t->share_limit[depth] != nullptr)
            share_from = t->share_limit[depth];
          else if (t->share_from[depth] != nullptr)
            share_from = t->share_from[depth];
          if (share_from != nullptr) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << depth << "] to " << (*share_from + 1) << std::endl;
          debug_mutex.unlock();
#endif
            rhs[depth] = *share_from + 1;
          } else {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << depth << "] to " << -CPX_INFBOUND << std::endl;
          debug_mutex.unlock();
#endif
            rhs[depth] = -CPX_INFBOUND;
          }
        }
        depth_level++;
        depth = t->perm(depth_level);
        if (sense == MIN) {
          if (t->share_limit[depth] && (*t->share_limit[depth] < max[depth])) {
#ifdef DEBUG
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
            std::cout << "setting rhs[" << depth << "] to " << (*t->share_limit[depth]+1) << std::endl;
            debug_mutex.unlock();
#endif
            rhs[depth] = *t->share_limit[depth] - 1;
          } else {
#ifdef DEBUG
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
            std::cout << "setting rhs[" << depth << "] to " << (max[depth]+1) << std::endl;
            debug_mutex.unlock();
#endif
            rhs[depth] = max[depth]-1;
          }
          max[depth] = (int) -CPX_INFBOUND;
        } else {
          if (t->share_limit[depth] && ((*t->share_limit[depth] > min[depth]) || (min[depth] == (int) CPX_INFBOUND))) {
#ifdef DEBUG
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
            std::cout << "setting rhs[" << depth << "] to " << (*t->share_limit[depth]+1) << std::endl;
            debug_mutex.unlock();
#endif
            rhs[depth] = *t->share_limit[depth] + 1;
          } else {
#ifdef DEBUG
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
            std::cout << "setting rhs[" << depth << "] to " << (min[depth]+1) << std::endl;
            debug_mutex.unlock();
#endif
            rhs[depth] = min[depth]+1;
          }
          min[depth] = (int) CPX_INFBOUND;
        }
        onwalk = true;
      } else if (!onwalk && infcnt != 1) {
        if (sense == MIN) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << depth << "] to " << (max[depth]-1) << std::endl;
          debug_mutex.unlock();
#endif
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        } else {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << depth << "] to " << (min[depth]+1) << std::endl;
          debug_mutex.unlock();
#endif
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
      } else if (onwalk && infcnt != 1)  {
        depth_level = 1;
        depth = t->perm(depth_level);
        if (sense == MIN) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << depth << "] to " << -CPX_INFBOUND << std::endl;
          debug_mutex.unlock();
#endif
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        } else {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << depth << "] to " << -CPX_INFBOUND << std::endl;
          debug_mutex.unlock();
#endif
          rhs[depth] = min[depth]+1;
          min[depth] = (int) CPX_INFBOUND;
        }
        onwalk = false;
      }
    }
  }
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "done!" << std::endl;
          debug_mutex.unlock();
#endif
  if (sharing) {
    completed = true;
    for(int i = 0; i < p.objcnt; ++i) {
      Locking_Vars *lv = nullptr;
      if (t->locks)
        lv = t->locks[i];
      if (lv != nullptr) {
        {
          std::unique_lock<std::mutex> lk(lv->status_mutex);
        }
        lv->cv.notify_all();
      }
    }
  }
#ifdef FINETIMING
  clock_gettime(CLOCK_MONOTONIC, &start);
  total_time = start.tv_sec + start.tv_nsec/1e9 - total_time;
  std::cout << "Thread " << t->id << " used " << cplex_time << "s in cplex";
  std::cout << ", waited for " << wait_time << "s";
  std::cout << " and " << total_time << "s overall." << std::endl;
#endif
  all.merge(s);
  delete[] resultStore;
  delete[] rhs;
  delete[] min;
  delete[] max;
}
