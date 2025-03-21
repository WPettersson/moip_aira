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

#if defined(DEBUG) || defined(DEBUG_SYNC) || defined(DEBUG_SHARES)
std::mutex debug_mutex;
#endif

extern std::string HASH;

int num_threads;
int cplex_threads;

/**
 * Are we splitting up the range of values an objective can take, such that
 * individual threads only find solutions in their own pre-determined range?
 */
bool split;

/**
 * If we are splitting, are we giving each thread the same sized range, or are
 * we assuming that the last objective is distributed normally, with a mean of
 * (highest+lowest)/2 and a std-dev of (highest-lowest)/6?
 */
bool split_normal;

/**
 * If we are assuming a normal distribution, here are the values we use. Only
 * have numbers for <= 12 threads, feel free to generate more.
 *
 * thread i, from X threads all together, has boundss normal_values[X][i] and
 * normal_values[X][i+1]
 */

double normal_values[13][13] {
  {0}, // 0 threads, to allow direct indexing
  {0, 1},
  {0, 0.5, 1},
  {0, 0.356, 0.644, 1},
  {0, 0.275, 0.5, 0.725, 1},
  {0, 0.219, 0.416, 0.584, 0.781, 1},
  {0, 0.178, 0.256, 0.5, 0.644, 0.822, 1},
  {0, 0.144, 0.311, 0.44, 0.56, 0.689, 0.856, 1},
  {0, 0.117, 0.275, 0.394, 0.5, 0.606, 0.725, 0.883, 1},
  {0, 0.093, 0.245, 0.356, 0.453, 0.547, 0.644, 0.755, 0.907, 1},
  {0, 0.073, 0.219, 0.325, 0.416, 0.5, 0.584, 0.675, 0.781, 0.927, 1},
  {0, 0.055, 0.197, 0.298, 0.384, 0.462, 0.538, 0.616, 0.702, 0.803, 0.945, 1},
  {0, 0.039, 0.178, 0.275, 0.356, 0.430, 0.5, 0.570, 0.644, 0.725, 0.822, 0.961, 1}
};

/**
 * Since normal values are pre-computed, we also have a limit on the number of
 * threads.
 */
const int MAX_THREADS_NORMAL_SPLIT = 12;

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
 * \param t A thread object, which itself includes the permutation, as well as
 * how many objectives to optimise.
 *
 * Note that objectives that are not optimised will still be calculated, they
 * just won't be optimised in any way.
 */
int solve(Env & e, Problem & p, int * result, double * rhs, Thread * t);

/**
 * Sets up the optimisation based on splitting the value of the final
 * objective.
 * \param e The relevant Env object
 * \param p The relevant Problem object
 * \param s A list in which to store solutions.
 * \param obj The objective to be optimised
 */
void split_setup(Env &e, Problem &p, std::list<int *> &s, int nObj);

/**
 * Finds out the limit (either maximum or minimum) of one objective in the
 * problem. This is only used to set up the splitting algorithm correctly.
 */
void get_limit(Env & e, Problem & p, int obj, double * rhs, int * result,
    const Sense sense);

/**
 * Optimise!
 * This function runs through the selected algorithm of Pettersson and
 * Ozlen[1]. It is designed to be run multi-threaded, hence using locking
 * mechanisms.
 *
 * \pram pFilename The file name of the file holding the problem description
 * \param all A Solutions object into which all solutions will be placed
 * \param infeasibles A Solutions object containing all infeasible problems
 * \param t A pointer to a Thread object that describes how this particular
 * thread should approach the problem.
 **/

/* Note that this function is templated. There's no evidence that this does
 * increase running time, but it's not like it decreases it either. */
template<Sense sense>
void optimise(const char * pFilename, Solutions & all, Solutions & infeasibles,
              Thread *t);

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
    ("split-normal,",
     po::bool_switch(&split_normal),
     "If splitting, assume a normal distribution for objective values\n"
     "Optional, defaults to False.")
    ("spread,s",
     po::bool_switch(&spread)->default_value(true),
     "Spread threads out over various subgroups of the symmetries (as opposed to clustering inside subgroups).\n"
     "Optional, defaults to True")
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

  if (split_normal && (num_threads > MAX_THREADS_NORMAL_SPLIT)) {
    std::cerr << "Error: split_normal can only handle at most "
      << MAX_THREADS_NORMAL_SPLIT << " threads." << std::endl;
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

  if (p.objcnt >= maxObjCount) {
    std::cerr << "Error: This version of moip_aira has been compiled to support at most " << maxObjCount << " objectives." << std::endl;
    return -ERR_AIRA;
  }

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

  outFile << std::endl << "Using improved algorithm at " << HASH << std::endl;

  /* Start the timer */
  starttime = clock();
  timespec start;
  clock_gettime(CLOCK_MONOTONIC, &start);
  startelapsed = start.tv_sec + start.tv_nsec/1e9;


  if (num_threads > S[p.objcnt].size())
    num_threads = S[p.objcnt].size();

  Solutions all(p.objcnt);

  std::list<Locking_Vars*> locking_var_list;

  std::list<Thread*> threads;
  if (split) {
    std::list<int *> sols;
    // Note that split_setup will farm out jobs to threads and everything.
    split_setup(e, p, sols, p.objcnt);
    double * lp = new double[p.objcnt];
    for(auto s : sols) {
      all.insert(lp, s, false /* infeasible */);
    }
  } else {
    // Not splitting.
    int * ordering = new int[p.objcnt];
    for (int c = 0; c < p.objcnt; ++c) {
      ordering[c] = c;
    }
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
  Solutions infeasibles(p.objcnt);
  std::list<std::thread> threadList;
  if (p.objsen == MIN) {
    for(Thread* thread: threads) {
      threadList.emplace_back(optimise<MIN>, pFilename.c_str(),
          std::ref(all), std::ref(infeasibles), thread);
    }
  } else {
    for(Thread* & thread: threads) {
      threadList.emplace_back(optimise<MAX>, pFilename.c_str(),
          std::ref(all), std::ref(infeasibles), thread);
    }
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
  for (auto& thread: threadList)
    thread.join();

  for (auto lv: locking_var_list)
    delete lv;

  /* Stop the clock. Sort and print results.*/
  endtime = clock();
  cpu_time_used=((double) (endtime - starttime)) / CLOCKS_PER_SEC;
  clock_gettime(CLOCK_MONOTONIC, &start);
  elapsedtime = (start.tv_sec + start.tv_nsec/1e9 - startelapsed);

  solcount = 0;
  all.sort_unique();

  for (const Result* soln: all) {
    if (soln->infeasible)
      continue;
    for (int i = 0; i < p.objcnt; i++) {
      outFile << soln->result[i] << "\t";
    }
    outFile << std::endl;
    solcount++;
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


void get_limit(Env & e, Problem & p, int obj, double * rhs, int * result, const Sense sense) {

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
    return;
  }
  status = CPXgetobjval (e.env, e.lp, &objval);
  if ( status ) {
    std::cerr << "Failed to obtain objective value." << std::endl;
    exit(0);
  }
  if ( abs(objval) > 1/p.mip_tolerance ) {
    while (abs(objval) > 1/p.mip_tolerance) {
      p.mip_tolerance /= 10;
    }
    CPXsetdblparam(e.env, CPXPARAM_MIP_Tolerances_MIPGap, p.mip_tolerance);
    status = CPXmipopt (e.env, e.lp);
    ipcount++;
    solnstat = CPXgetstat (e.env, e.lp);
    if ((solnstat == CPXMIP_INFEASIBLE) || (solnstat == CPXMIP_INForUNBD)) {
      return;
    }
    status = CPXgetobjval (e.env, e.lp, &objval);
    if ( status ) {
      std::cerr << "Failed to obtain objective value." << std::endl;
      exit(0);
    }
  }
  delete[] srhs;

  // Get the solution vector
  double * soln = new double[cur_numcols];
  CPXgetx(e.env, e.lp, soln, 0, cur_numcols - 1);
  // Now run through the rest of the objectives.
  for (int j = 0; j < p.objcnt; j++) {
    double res = 0;
    for(int i = 0; i < cur_numcols; ++i) {
      res += p.objcoef[j][i] * soln[i];
    }
    result[j] = round(res);
  }
  delete[] soln;
  return;
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
  for (int j_preimage = 0; j_preimage < t->nObj(); j_preimage++) {
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
    if ( abs(objval) > 1/p.mip_tolerance ) {
      while (abs(objval) > 1/p.mip_tolerance) {
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
  // Get the solution vector
  double * soln = new double[cur_numcols];
  CPXgetx(e.env, e.lp, soln, 0, cur_numcols - 1);
  // Now run through the rest of the objectives.
  for (int j_pre = t->nObj(); j_pre < p.objcnt; j_pre++) {
    int j = t->perm(j_pre);
    double res = 0;
    for(int i = 0; i < cur_numcols; ++i) {
      res += p.objcoef[j][i] * soln[i];
    }
    result[j] = round(res);
  }
  delete[] soln;

  delete[] srhs;

  return solnstat;
}

template<Sense sense>
void optimise(const char * pFilename, Solutions & all, Solutions & infeasibles,
    Thread *t) {
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
 #ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " splitting on ";
          std::cout << t->perm(t->nObj()-1) << std::endl;
          debug_mutex.unlock();
#endif
    rhs[t->perm(t->nObj()-1)] = t->split_start;
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
  if (solnstat == CPXMIP_INFEASIBLE) {
    std::cout << "infeasible";
  } else {
    for(int i = 0; i < p.objcnt; ++i) {
      std::cout << result[i] << ",";
    }
  }
  std::cout << std::endl;
  debug_mutex.unlock();
#endif

  /* Need to add a result to the list here*/
  if (solnstat == CPXMIP_INFEASIBLE) {
    infeasibles.insert(rhs, result, true);
  } else {
    if (split)
      all.insert(rhs, result, solnstat == CPXMIP_INFEASIBLE);
    else
      s.insert(rhs, result, solnstat == CPXMIP_INFEASIBLE);
  }
  // Note that if we are splitting, we aren't sharing.
  if (split) {
    if (p.objsen == MIN)
      t->split_stop--;
    else
      t->split_stop++;
  } else {
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
  for (int objective_counter = 1; objective_counter < t->nObj(); objective_counter++) {
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
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << j << "] to ±inf" << std::endl;
          debug_mutex.unlock();
#endif
      } else {
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << j << "] to " << *t->share_from[j] << std::endl;
          debug_mutex.unlock();
#endif
        rhs[j] = *t->share_from[j];
      }
    }
    if (split) {
        rhs[t->perm(t->nObj()-1)] = t->split_start;
    }
    /* Set rhs of current depth */
    if (sense == MIN) {
      rhs[objective] = max[objective]-1;
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << objective << "] to " << (max[objective]-1) << std::endl;
          debug_mutex.unlock();
#endif
    } else {
      rhs[objective] = min[objective]+1;
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << objective << "] to " << (min[objective]+1) << std::endl;
          debug_mutex.unlock();
#endif
    }
    if (split) {
      // check if we cross midpoint
      if (sense == MIN) {
        if (rhs[t->nObj()-1] < t->split_stop) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " reached split_stop";
          std::cout << " " << t->split_stop << std::endl;
          debug_mutex.unlock();
#endif
          break;
        }
      } else {
        if (rhs[t->nObj()-1] > t->split_stop) {
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " reached split_stop";
          std::cout << " " << t->split_stop << std::endl;
          debug_mutex.unlock();
#endif
          break;
        }
      }
    }
    max[objective] = (int) -CPX_INFBOUND;
    min[objective] = (int) CPX_INFBOUND;
    while (infcnt < objective_counter) {
      bool relaxed;
      int solnstat;
      /* Look for possible relaxations to the current problem*/
      const Result *relaxation;

#ifdef DEBUG_SOLUTION_SEARCH
      debug_mutex.lock();
      std::cout << "Thread " << t->id;
      debug_mutex.unlock();
#endif
      // First check if it's infeasible
      relaxation = infeasibles.find(rhs, p.objsen);
      if (relaxation == nullptr) {
        if (split) {
          relaxation = all.find(rhs, p.objsen);
        } else {
          relaxation = s.find(rhs, p.objsen);
        }
      }
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
        if (infeasible) {
          infeasibles.insert(rhs, result, true);
        } else {
          if (split) {
            all.insert(rhs, result, infeasible);
          } else {
            s.insert(rhs, result, infeasible);
          }
        }
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
        std::cout << "infeasible ";
      } else {
        for(int i = 0; i < p.objcnt; ++i) {
          std::cout << result[i] << ",";
        }
      }
      std::cout << " on line " << __LINE__ << std::endl;
      debug_mutex.unlock();
#endif
      if (split) {
        if (!infeasible) {
          // check if we cross midpoint
          if (infcnt == t->nObj() - 2) {
            if (sense == MIN) {
              if (rhs[t->nObj()-1] < t->split_stop) {
                infeasible = true;
#ifdef DEBUG
                debug_mutex.lock();
                std::cout << "Thread " << t->id << " reached split_stop";
                std::cout << " " << t->split_stop << std::endl;
                debug_mutex.unlock();
#endif
              }
            } else {
              if (rhs[t->nObj()-1] > t->split_stop) {
                infeasible = true;
#ifdef DEBUG
                debug_mutex.lock();
                std::cout << "Thread " << t->id << " reached split_stop";
                std::cout << " " << t->split_stop << std::endl;
                debug_mutex.unlock();
#endif
              }
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
      } else if (sharing && t->locks && t->locks[t->perm(infcnt+1)]) {
        std::unique_lock<std::mutex> lk(t->locks[t->perm(infcnt+1)]->status_mutex);
#ifdef DEBUG
          auto lock_num = t->perm(infcnt+1);
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << " got lock on " << lock_num << std::endl;
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
          std::cout << "apparently partner found something ";
          std::cout << "on obj " << obj << std::endl;
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
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << " should drop lock on " << lock_num << std::endl;
          debug_mutex.unlock();
#endif
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
        bool wait = false;
        Locking_Vars *lv = nullptr;
        if (t->locks)
          lv = t->locks[updated_objective];
        if ( lv != nullptr && !lv->any_complete()) {
          std::unique_lock<std::mutex> lk(lv->status_mutex);
          t->state = ThreadState::Waiting;
          if ( ! lv->all_done()) {
#ifdef DEBUG_SYNC
            debug_mutex.lock();
            std::cout << "Thread " << t->id <<  " done, at " << __LINE__ << std::endl;
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
            if (!lv->any_complete()) {
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "waiting on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
              lv->cv_subproblem_complete.wait(lk); // Wait for all other threads to update limits.
              t->state = ThreadState::Running;
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "released on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
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
            // This thread is last to finish this subproblem.
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
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting max[" << i << "] to " << *t->share_bounds[i] << " from " << max[i] << std::endl;
          debug_mutex.unlock();
#endif
                    max[i] = *t->share_bounds[i];
                  }
                } else {
                  if (*t->share_bounds[i] > min[i]) {
                    *t->share_bounds[i] = min[i];
                  } else {
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting min[" << i << "] to " << min[i] << std::endl;
          debug_mutex.unlock();
#endif
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
#ifdef DEBUG_SHARES
              debug_mutex.lock();
              std::cout << "Thread " << t->id << " ";
              std::cout << "set limit[" << updated_objective << "] to ";
              std::cout << *limit << std::endl;
              debug_mutex.unlock();
#endif
            }
#ifdef DEBUG_SHARES
            debug_mutex.lock();
            std::cout << "Thread " << t->id << ", ";
            std::cout << "max is " << max[updated_objective] << ", ";
            if (t->share_bounds[updated_objective] != nullptr)
              std::cout << "share_to is " << *t->share_to[updated_objective] << std::endl;
            debug_mutex.unlock();
#endif
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "releasing on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
            lk.unlock();
            t->state = ThreadState::Running;
            lv->cv_subproblem_complete.notify_all();
          }
        }
        if (lv != nullptr && !lv->any_complete()) {
          for(int i = 0; i <= infcnt; ++i) {
            if (sense == MIN) {
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "setting max[" << t->perm(i) << "] to -inf" << std::endl;
          debug_mutex.unlock();
#endif
              max[t->perm(i)] = (int)-CPX_INFBOUND;
              if (t->share_to[t->perm(i)] != nullptr) {
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "setting share_to[" << t->perm(i) << "] to inf" << std::endl;
          debug_mutex.unlock();
#endif
                *t->share_to[t->perm(i)] = (int)CPX_INFBOUND;
              }
              if (t->share_limit[t->perm(i)] != nullptr) {
#ifdef DEBUG_SHARES
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
              t->state = ThreadState::Waiting;
              if (! lv->all_done()) {
                if (lv->any_complete())
                  break;
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "waiting on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
                lv->cv.wait(lk);
                t->state = ThreadState::Running;
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "released on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
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
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "releasing on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
                t->state = ThreadState::Running;
                lk.unlock();
                lv->cv.notify_all();
              }
            }
          do {
            {
              std::unique_lock<std::mutex> lk(lv->status_mutex);
              t->state = ThreadState::Waiting;
              if (! lv->all_done()) {
                if (lv->any_complete())
                  break;
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "waiting on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
                lv->cv.wait(lk);
                t->state = ThreadState::Running;
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "released on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
              } else {
                // Last in
                lv->changed = false;
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "releasing on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
                lk.unlock();
                t->state = ThreadState::Running;
                lv->cv.notify_all();
              }
            }
            for (int i = 0; i <= infcnt; ++i) {
              int obj = t->perm(i);
              if (sense == MIN) {
                if (t->share_from[obj] != nullptr) {
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << " share_from[" << obj << "] is " << *t->share_from[obj] << std::endl;
          debug_mutex.unlock();
#endif
                  if (t->share_limit[obj] != nullptr) {
                    if (*t->share_limit[obj] > *t->share_from[obj]) {
                      lv->changed = true;
                      *t->share_limit[obj] = *t->share_from[obj];
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "setting share_limit[" << obj << "] to " << max[obj] << std::endl;
          debug_mutex.unlock();
#endif
                    }
                  }
                  if (t->share_to[obj] != nullptr) {
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << " share_to[" << obj << "] is " << *t->share_to[obj] << std::endl;
          debug_mutex.unlock();
#endif
                    if (*t->share_to[obj] > *t->share_from[obj]) {
                      lv->changed = true;
                      *t->share_to[obj] = *t->share_from[obj];
#ifdef DEBUG_SHARES
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
#ifdef DEBUG_SHARES
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
#ifdef DEBUG_SHARES
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
              t->state = ThreadState::Waiting;
              if (! lv->all_done()) {
                if (lv->any_complete())
                  break;
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "waiting on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
                lv->cv.wait(lk);
                t->state = ThreadState::Running;
#ifdef DEBUG_SYNC
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " ";
          std::cout << "released on lock " << lv << " at " << __LINE__ << std::endl;
          debug_mutex.unlock();
#endif
              } else {
                // Last in
                t->state = ThreadState::Running;
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
#ifdef DEBUG_SHARES
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
#ifdef DEBUG_SHARES
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
        if (sharing && (p.objcnt > 2) && (objective_counter == p.objcnt - 1)) {
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
#ifdef DEBUG_SHARES
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
#ifdef DEBUG_SHARES
              debug_mutex.lock();
              std::cout << "Thread " << t->id << " setting rhs[" << j;
              std::cout << "] to " << (*share_from - 1) << std::endl;
              debug_mutex.unlock();
#endif
              rhs[j] = *share_from - 1;
            } else {
#ifdef DEBUG_SHARES
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
#ifdef DEBUG_SHARES
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " setting share_to[" << j;
            std::cout << "] to " << *t->share_to[j] << std::endl;
            debug_mutex.unlock();
#endif
                }
              } else {
                if (*t->share_to[j] < *share_from) {
#ifdef DEBUG_SHARES
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
        /* If splitting, use split_start and not +/- infinity */
        if (split) {
          rhs[t->nObj()-1] = t->split_start;
        }
        /* In the case of a minimisation problem
         * set current level to max objective function value  -1 else set
         * current level to min objective function value  +1 */
        if (sense == MIN) {
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << objective << "] to " << max[objective]-1 <<std::endl;
          debug_mutex.unlock();
#endif
          rhs[objective] = max[objective]-1;
          max[objective] = (int) -CPX_INFBOUND;
        } else {
#ifdef DEBUG_SHARES
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
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "increasing depth from " << depth_level;
          std::cout << " when infcnt is " << infcnt << " and objective_counter is " << objective_counter << std::endl;
          debug_mutex.unlock();
#endif
        if (sharing) {
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
#ifdef DEBUG_SHARES
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
            std::cout << "setting rhs[" << depth << "] to " << (*share_from + 1) << std::endl;
            debug_mutex.unlock();
#endif
              rhs[depth] = *share_from + 1;
            } else {
#ifdef DEBUG_SHARES
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
            std::cout << "setting rhs[" << depth << "] to " << -CPX_INFBOUND << std::endl;
            debug_mutex.unlock();
#endif
              rhs[depth] = -CPX_INFBOUND;
            }
          }
        } else {
          if (sense == MIN) {
            rhs[depth] = CPX_INFBOUND;
          } else {
            rhs[depth] = -CPX_INFBOUND;
          }
        }

        depth_level++;
        depth = t->perm(depth_level);
#ifdef DEBUG_SHARES
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
            std::cout << "setting depth to " << depth << " and depth_level to " << depth_level << std::endl;
            debug_mutex.unlock();
#endif
        if (sense == MIN) {
          // Use the shared limit if the shared limit is smaller than max, or
          // if max is infinity, which happens if this thread found nothing
          // last round but another thread did find something.
          if (sharing && t->share_limit[depth] &&
              ((*t->share_limit[depth] < max[depth]) || (max[depth] == (int) -CPX_INFBOUND))) {
#ifdef DEBUG_SHARES
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
            std::cout << "setting rhs[" << depth << "] to " << (*t->share_limit[depth]+1) << std::endl;
            debug_mutex.unlock();
#endif
            rhs[depth] = *t->share_limit[depth] - 1;
          } else {
#ifdef DEBUG_SHARES
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
            std::cout << "setting rhs[" << depth << "] to " << (max[depth]+1) << std::endl;
            debug_mutex.unlock();
#endif
            rhs[depth] = max[depth]-1;
          }
          max[depth] = (int) -CPX_INFBOUND;
        } else {
          if (sharing && t->share_limit[depth] &&
              ((*t->share_limit[depth] > min[depth]) || (min[depth] == (int)CPX_INFBOUND))) {
#ifdef DEBUG_SHARES
            debug_mutex.lock();
            std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
            std::cout << "setting rhs[" << depth << "] to " << (*t->share_limit[depth]+1) << std::endl;
            debug_mutex.unlock();
#endif
            rhs[depth] = *t->share_limit[depth] + 1;
          } else {
#ifdef DEBUG_SHARES
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
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << " when infcnt is " << infcnt << " and objective_counter is " << objective_counter << std::endl;
          debug_mutex.unlock();
#endif
        if (sense == MIN) {
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << depth << "] to " << (max[depth]-1) << std::endl;
          debug_mutex.unlock();
#endif
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        } else {
#ifdef DEBUG_SHARES
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
#ifdef DEBUG_SHARES
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "setting rhs[" << depth << "] to " << -CPX_INFBOUND << std::endl;
          debug_mutex.unlock();
#endif
          rhs[depth] = max[depth]-1;
          max[depth] = (int) -CPX_INFBOUND;
        } else {
#ifdef DEBUG_SHARES
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
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << " just before end of loop with objective_counter = ";
          std::cout << objective_counter << " and infcnt = " << infcnt << std::endl;
          debug_mutex.unlock();
#endif
    }
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << " at end of loop with objective_counter = ";
          std::cout << objective_counter << " and infcnt = " << infcnt << std::endl;
          debug_mutex.unlock();
#endif

  }
#ifdef DEBUG
          debug_mutex.lock();
          std::cout << "Thread " << t->id << " at " << __LINE__ << " ";
          std::cout << "done!" << std::endl;
          debug_mutex.unlock();
#endif
  if (sharing) {
    t->state = ThreadState::Complete;
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
  s.sort_unique();

  all.merge(s);
  delete[] resultStore;
  delete[] rhs;
  delete[] min;
  delete[] max;
}

void split_optimise(Problem &p, std::list<int *>& s, int nObj, int max, int min) {
  double start_point, stop_point;
  if (p.objsen == MIN) {
    start_point = max;
    stop_point = min;
  } else {
    start_point = min;
    stop_point = max;
  }
  double split_start = start_point;
  double split_stop;
  double step_size = (stop_point - start_point)/num_threads;
  std::list<Thread *> threads;
  for (int t = 0; t < num_threads; ++t) {
    if (split_normal) {
      double start, stop;
      if (p.objsen == MIN) {
        double gap = (start_point - stop_point);
        stop = normal_values[num_threads][t]*gap + stop_point;
        start = normal_values[num_threads][t+1]*gap + stop_point;
      } else {
        double gap = (stop_point - start_point);
        start = normal_values[num_threads][t]*gap + start_point;
        stop = normal_values[num_threads][t+1]*gap + start_point;
      }
      threads.push_back(new Thread(t, nObj, p.objcnt, start, stop));
    } else {
      split_stop = split_start + step_size;
      threads.push_back(new Thread(t, nObj, p.objcnt, split_start, split_stop));
      split_start = split_stop;
    }
  }
  Solutions here(p.objcnt);
  Solutions infeasibles(p.objcnt);
  std::list<std::thread> threadList;
  if (p.objsen == MIN) {
    for(Thread* thread: threads) {
      threadList.emplace_back(optimise<MIN>, p.filename(),
          std::ref(here), std::ref(infeasibles), thread);
    }
  } else {
    for(Thread* & thread: threads) {
      threadList.emplace_back(optimise<MAX>, p.filename(),
          std::ref(here), std::ref(infeasibles), thread);
    }
  }
  for (auto& thread: threadList)
    thread.join();
  for(auto sol: here) {
    if (sol->infeasible)
      continue;
    int * n = new int[p.objcnt];
    for(int i = 0; i < p.objcnt; ++i) {
      n[i] = sol->result[i];
    }
    s.push_back(n);
  }
}

void split_setup(Env &e, Problem &p, std::list<int *> &s, int nObj) {
  if (nObj == 1) {
    int *min = new int[p.objcnt];
    get_limit(e, p, nObj - 1, p.rhs, min, p.objsen);
    s.push_back(min);
  } else {
#ifdef DEBUG
    std::cout << "Splitting on " << nObj << " objectives ... ";
#endif
    std::list<int *> sols;
    split_setup(e, p, sols, nObj - 1);
    int * res = new int[p.objcnt];
    int biggest, smallest;
    if (p.objsen == MIN) {
      get_limit(e, p, nObj-1, p.rhs, res, p.objsen);
      smallest = res[nObj - 1];
      biggest = (int)-CPX_INFBOUND;
      for(int * sol: sols) {
        if (sol[nObj-1] > biggest) {
          biggest = sol[nObj-1];
        }
      }
      if (biggest == smallest) {
        biggest = (int)CPX_INFBOUND;
      }
    } else {
      get_limit(e, p, nObj-1, p.rhs, res, p.objsen);
      biggest = res[nObj - 1];
      smallest = (int)CPX_INFBOUND;
      for(int * sol: sols) {
        if (sol[nObj-1] < smallest) {
          smallest = sol[nObj-1];
        }
      }
      if (biggest == smallest) {
        smallest = (int)-CPX_INFBOUND;
      }
    }
    delete[] res;
#ifdef DEBUG
    std::cout << " found range [" << smallest << ", " << biggest << "]";
    std::cout << " for objective " << (nObj-1) << std::endl;
#endif
    split_optimise(p, s, nObj, biggest, smallest);
  }
}
